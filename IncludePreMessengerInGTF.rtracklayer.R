library("rtracklayer")
library("Biostrings")

#' @param RessourceDIR Location of the raw ressources files directory
#' @param MMGTF name of the Mus Musculus GTF, supposed to be in ressource directory
#' @param MMfasta name of the Mus Musculus Fasta, supposed to be in ressource directory
#' @param DmelGTF Location of the Spike-In GTF, supposed to be in ressource directory
#' @param DmelFasta Location of the Spike-In Fasta , supposed to be in ressource directory
#' @param NewRefDIR Location of the Combined ressources files directory
RessourcesDIR<-"./TestInput/"
MMGTF<-"gencode.vM25.annotation.gtf"
MMfasta<-"GRCm38.primary_assembly.genome.fa"
DmelGTF<-"dmel-all-r6.60.gtf"
DmelFasta<-"dmel-all-chromosome-r6.60.fasta"
NewRefDIR<-"./TestOutput/"


## GTF Part

GTF.GR <- rtracklayer::import(paste0(RessourcesDIR,MMGTF))

# Cosmetic discarding
GTF.GR <- GTF.GR[GTF.GR$type!="CDS",]
GTF.GR <- GTF.GR[GTF.GR$type!="start_codon",]
GTF.GR <- GTF.GR[GTF.GR$type!="stop_codon",]

# supposed to be https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/gtf/dmel-all-r6.60.gtf.gz
FlyGTF.GR<-rtracklayer::import(paste0(RessourcesDIR,DmelGTF))

# Some genes have no strand
FlyGTF.GR<-FlyGTF.GR[strand(FlyGTF.GR)!="*"] # Strand inconnu

# Some genes are located on both strands
Strand_Gene<-unique(data.frame(FlyGTF.GR)[,c("strand","gene_id")])
Aberrant_Genes<-Strand_Gene$gene_id[duplicated(Strand_Gene$gene_id)]

FlyGTF.GR<-FlyGTF.GR[!mcols(FlyGTF.GR)$gene_id %in% Aberrant_Genes]
#SpikeInDmelChr3R       FlyBase mRNA    21360390        21377399        .       .       .       gene_id "SpikeInDmelFBgn0002781"; transcript_id "FBtr0084079"; gene_symbol "mod(mdg4)"; transcript_symbol "mod(mdg4)-RT";
#Error Message: Strand is neither '+' nor '-'!

# modification of chromosome (seqlevels) to add a flag as prefix
seqlevels(FlyGTF.GR)<-paste0("SpikeInDmelChr",seqlevels( FlyGTF.GR))
# modification of gene names (gene_id) to add a flag as prefix
FlyGTF.GR$gene_id<-paste0("SpikeInDmel",FlyGTF.GR$gene_id)
# modification of transcript names(transcript_id) to add a flag as prefix
FlyGTF.GR$transcript_id[!is.na(FlyGTF.GR$transcript_id)]<-paste0("SpikeInDmel",FlyGTF.GR$transcript_id[!is.na(FlyGTF.GR$transcript_id)])

FlyGTF.GR<-sort(FlyGTF.GR)

#Combine both GTF
NewGTF<-c(GTF.GR,FlyGTF.GR)


## FASTA PART

MouseGenome<-Biostrings::readDNAStringSet(paste0(RessourcesDIR,MMfasta))
DmelGenome<-Biostrings::readDNAStringSet(paste0(RessourcesDIR,DmelFasta))
# Add a flag on chromosome on the chromosome
names(DmelGenome)<-paste0("SpikeInDmelChr",names(DmelGenome))
MousePlusDmel<-c(MouseGenome,DmelGenome)

# simple test that gene can be extracted from genome
MousePlusDmelForTest<-MousePlusDmel
# Biostrings require exact correspondance between the seqnames and the names in the fasta file
names(MousePlusDmelForTest)<-gsub(" .*","",names(MousePlusDmelForTest))
# Just to avoid crash in R we will restrict on common chromosome between fasta file and GTF
CommonCHRs<-GenomicRanges::intersect(unique(as.vector(seqnames(NewGTF))),unique(names(MousePlusDmelForTest)))

NewGTF<-NewGTF[seqnames(NewGTF) %in% CommonCHRs]
MousePlusDmelForTest<-MousePlusDmelForTest[names(MousePlusDmelForTest) %in% CommonCHRs]

# The test per se, if it works it should have a biostring output equivalent to each entry of the GTF
MousePlusDmelForTest[NewGTF] 


# export the file
NewGTFNAME<-paste0(NewRefDIR,gsub(".gtf$","",MMGTF),"_",gsub(".gtf$","",DmelGTF),".SpikeFlag.gtf")
export(NewGTF,NewGTFNAME)
NewFastaNAME<-paste0(NewRefDIR,gsub(".fa$","",MMfasta),"_",gsub(".fasta$","",DmelFasta),".SpikeFlag.fasta")
writeXStringSet(MousePlusDmel,NewFastaNAME)

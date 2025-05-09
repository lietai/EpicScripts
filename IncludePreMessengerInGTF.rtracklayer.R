# add pre messenger for mus musculus
# add dmel
# tag dmel as spikein

library("rtracklayer")
library("Biostrings")
RessourcesDIR<-"../Ressources/RAW/"
#MMGTF<-"gencode.vM25.annotation.gtf"
MMGTF<-"mm10-rDNA_v1.0.gtf"
MMfasta<-"mm10-rDNA_v1.0.fa"
DmelGTF<-"dmel-all-r6.60.gtf"
DmelFasta<-"dmel-all-chromosome-r6.60.fasta"
NewRefDIR<-"../Ressources/BulkRef/"

GTF.GR <- rtracklayer::import(paste0(RessourcesDIR,MMGTF))

#::Begening of the mm10-rDNA gtf specific part

#Issues with UCSC like GTF:
# - No Gene type
# - issues with gene_id naming (they basically use the same name as the transcriptID)
# - issues with transcript_id naming for the multi transcripts gene they add _[0-9]*
#for the second and + transcript of the gene but not for the first
# - issues with ribosomal chromosome coordinate

#Remarks (should not be important)
# - exon id is transcript base meaning that two identical GenomicRanges belonging 
#to two distinct transcript will have two exon_id identifier

#Removing the protein information, too complicated for now
GTF.GR <- GTF.GR[GTF.GR$type!="CDS",]
GTF.GR <- GTF.GR[GTF.GR$type!="start_codon",]
GTF.GR <- GTF.GR[GTF.GR$type!="stop_codon",]

# Issues with coordonates of the ribosomal chromosome
GTF.GR[seqnames(GTF.GR)=="chrR" & GTF.GR$type=="exon"]<-resize(GTF.GR[seqnames(GTF.GR)=="chrR" & GTF.GR$type=="exon"],width(GTF.GR[seqnames(GTF.GR)=="chrR" & GTF.GR$type=="exon"])-1)

#issues with naming of the genes in ribosomal pseudo chromosome
mcols(GTF.GR[seqnames(GTF.GR)!="chrR"])$gene_id<-gsub("_.*$","",mcols(GTF.GR[seqnames(GTF.GR)!="chrR"])$gene_id)

#issues with genes that have only one transcript or the first transcript of the several transcript
#they are names like the gene while other have _[0-9]* at the end
#this will cause and issue for the interpretation of the results
UnfinishedNamedTranscript<-grep("_[0-9]*",GTF.GR$transcript_id,inver=TRUE)
GTF.GR$transcript_id[UnfinishedNamedTranscript]<-paste0(GTF.GR$transcript_id[UnfinishedNamedTranscript],"_1")

#Renaming for clarity
GTF.GR$transcript_id<-gsub("_","_Transcript",GTF.GR$transcript_id)

#Autonomous part of the sexual chromosome, Discared for now
chrXYGenes<-intersect(unique(mcols(GTF.GR[seqnames(GTF.GR)=="chrX"])$gene_id),unique(mcols(GTF.GR[seqnames(GTF.GR)=="chrY"])$gene_id))
GTF.GR<-GTF.GR[!GTF.GR$gene_id %in% chrXYGenes]



#Create a gene feature for the gtf
#gene will englobe the most extreme coordinates of the transcripts that rely on it
Transcripts.GR<-GTF.GR[GTF.GR$type=="transcript"]
Transcripts_by_gene.GR<-split(Transcripts.GR,Transcripts.GR$gene_id)
gene.GR<-lapply(Transcripts_by_gene.GR,reduce)
gene.GR<-unlist(as(gene.GR,"GRangesList"))
gene.GR$gene_id<-names(gene.GR)
gene.GR$type<-"gene"
names(gene.GR)<-NULL

#Reinject the Gene into the GTF
GTF.GR<-c(GTF.GR,gene.GR)
GTF.GR<-sort(GTF.GR)

#:: end of the mm10-rDNA gtf specific part



#Export to GTF
#https://support.bioconductor.org/p/70054/#70146

# Get MonoExonic Gene
EXONS.GR<-GTF.GR[GTF.GR$type=="exon"]
ExonsXGenes<-unique(mcols(EXONS.GR[,c("gene_id","exon_id")]))
NExonsXGenes<-table(ExonsXGenes$gene_id)
MonoExonicGenes<-names(NExonsXGenes)[NExonsXGenes==1]
#write.table(MonoExonicGenes,"MonoExonicGenes.txt",quote=FALSE,row.names=FALSE)

# Get the non Monoexonic Gene create a transcript flagged PreMRNA
# The preMRNA is based on the gene coordinate and not the transcript
# Its just a way to capture the intronic reads
NoMonoExonic.Genes.GR <- GTF.GR[GTF.GR$type=="gene" & !GTF.GR$gene_id %in% MonoExonicGenes]
NoMonoExonic.PreMessager.transcript<- NoMonoExonic.Genes.GR
NoMonoExonic.PreMessager.transcript$type<-"transcript"
NoMonoExonic.PreMessager.transcript$transcript_id<-paste0("PreMRNA_",NoMonoExonic.PreMessager.transcript$gene_id)

# Create artificial exon that match the PreMRNA
NoMonoExonic.PreMessager.exon<- NoMonoExonic.Genes.GR
NoMonoExonic.PreMessager.exon$type<-"exon"
NoMonoExonic.PreMessager.exon$transcript_id<-paste0("PreMRNA_",NoMonoExonic.PreMessager.transcript$gene_id)
NoMonoExonic.PreMessager.exon$exon_id<-paste0("PreMRNA_",NoMonoExonic.PreMessager.transcript$gene_id)
NoMonoExonic.PreMessager.exon$exon_number<-1

#Combine everything
NewGTF<-c(GTF.GR,NoMonoExonic.PreMessager.exon,NoMonoExonic.PreMessager.transcript)
NewGTF<-sort(NewGTF)
#export(NewGTF,"test.gtf")


FlyGTF.GR<-rtracklayer::import(paste0(RessourcesDIR,DmelGTF))
FlyGTF.GR<-FlyGTF.GR[strand(FlyGTF.GR)!="*"] # Strand inconnu
Strand_Gene<-unique(data.frame(FlyGTF.GR)[,c("strand","gene_id")])
Aberrant_Genes<-Strand_Gene$gene_id[duplicated(Strand_Gene$gene_id)]
FlyGTF.GR<-FlyGTF.GR[!mcols(FlyGTF.GR)$gene_id %in% Aberrant_Genes]
#SpikeInDmelChr3R       FlyBase mRNA    21360390        21377399        .       .       .       gene_id "SpikeInDmelFBgn0002781"; transcript_id "FBtr0084079"; gene_symbol "mod(mdg4)"; transcript_symbol "mod(mdg4)-RT";
#Error Message: Strand is neither '+' nor '-'!
seqlevels(FlyGTF.GR)<-paste0("SpikeInDmelChr",seqlevels( FlyGTF.GR))
FlyGTF.GR$gene_id<-paste0("SpikeInDmel",FlyGTF.GR$gene_id)
FlyGTF.GR$transcript_id[!is.na(FlyGTF.GR$transcript_id)]<-paste0("SpikeInDmel",FlyGTF.GR$transcript_id[!is.na(FlyGTF.GR$transcript_id)])

FlyGTF.EXONS.GR<-FlyGTF.GR[FlyGTF.GR$type=="exon"]
FlyGTF.EXONS.GR$exon_id<-paste(seqnames(FlyGTF.EXONS.GR),start(FlyGTF.EXONS.GR),end(FlyGTF.EXONS.GR),strand(FlyGTF.EXONS.GR),sep="_")
FlyGTF.ExonsXGenes<-unique(mcols(FlyGTF.EXONS.GR[,c("gene_id","exon_id")]))
FlyGTF.NExonsXGenes<-table(FlyGTF.ExonsXGenes$gene_id)
FlyGTF.MonoExonicGenes<-names(FlyGTF.NExonsXGenes)[FlyGTF.NExonsXGenes==1]

FlyGTF.NoMonoExonic.Genes.GR <- FlyGTF.GR[FlyGTF.GR$type=="gene" & !FlyGTF.GR$gene_id %in% FlyGTF.MonoExonicGenes]
FlyGTF.NoMonoExonic.PreMessager.transcript<- FlyGTF.NoMonoExonic.Genes.GR
FlyGTF.NoMonoExonic.PreMessager.transcript$type<-"transcript"
FlyGTF.NoMonoExonic.PreMessager.transcript$transcript_id<-paste0("PreMRNA_",FlyGTF.NoMonoExonic.PreMessager.transcript$gene_id)

FlyGTF.NoMonoExonic.PreMessager.exon<- FlyGTF.NoMonoExonic.Genes.GR
FlyGTF.NoMonoExonic.PreMessager.exon$type<-"exon"
FlyGTF.NoMonoExonic.PreMessager.exon$transcript_id<-paste0("PreMRNA_",FlyGTF.NoMonoExonic.PreMessager.transcript$gene_id)
FlyGTF.NoMonoExonic.PreMessager.exon$exon_id<-paste0("PreMRNA_",FlyGTF.NoMonoExonic.PreMessager.transcript$gene_id)
FlyGTF.NoMonoExonic.PreMessager.exon$exon_number<-1

FlyGTF.GR<-c(FlyGTF.GR,FlyGTF.NoMonoExonic.PreMessager.exon,FlyGTF.NoMonoExonic.PreMessager.transcript)
FlyGTF.GR<-sort(FlyGTF.GR)


NewGTF<-c(NewGTF,FlyGTF.GR)




MouseGenome<-readDNAStringSet(paste0(RessourcesDIR,MMfasta))
DmelGenome<-readDNAStringSet(paste0(RessourcesDIR,DmelFasta))
names(DmelGenome)<-paste0("SpikeInDmelChr",names(DmelGenome))
MousePlusDmel<-c(MouseGenome,DmelGenome)


#simple test that gene can be extracted from genome
MousePlusDmelForTest<-MousePlusDmel
names(MousePlusDmelForTest)<-gsub(" .*","",names(MousePlusDmelForTest))

CommonCHRs<-intersect(unique(as.vector(seqnames(NewGTF))),unique(names(MousePlusDmelForTest)))
NewGTF<-NewGTF[seqnames(NewGTF) %in% CommonCHRs]
MousePlusDmelForTest<-MousePlusDmelForTest[names(MousePlusDmelForTest) %in% CommonCHRs]
MousePlusDmelForTest[NewGTF] 


export(NewGTF,paste0(NewRefDIR,gsub(".gtf$","",MMGTF),".preMRNA",gsub(".gtf$","",DmelGTF),".SpikeFlag.gtf"))
writeXStringSet(MousePlusDmel,paste0(NewRefDIR,gsub(".fa$","",MMfasta),".",gsub(".fasta$","",DmelFasta),".SpikeFlag.fasta"))

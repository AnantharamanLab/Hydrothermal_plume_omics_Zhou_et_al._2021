# ----------------------------------------------------------------------#
# This R programme is for DEGseq analysis, need five parameters
# the 1st parameter: the input read count file
# the 2nd parameter: the control sample,seperated by ','
# the 3rd parameter: the treat sample, seperated by ','
# the 4th parameter: the output directory with prefix
# the 5th parameter: the the pair information for each sample, seperated by ','
# joeyxu@
# date: 2015-4-28
# 
#-----------------------------------------------------------------------#
#options <- commandArgs(trailingOnly = T)
infile = "Functional_analysis_summary.mdf.txt";
control = "CymD.D.LB,CymS.D.B";
treat = "LBAb.D.NBB,LBTa.D.APB,LBMa.D.BPB";
output = "Functional_analysis_summary.B.Cym-v-Lau";
#pair="1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8";
outfile =paste(output,".diff.xls",sep="");

library("DESeq2");

data=read.table(infile,sep="\t",header=T,row.names=1);
data=data[rowSums(data) != 0,];
#colnames(data)<-gsub("X","",colnames(data))
samples=names(data)

controls=unlist(strsplit(control,","))
treats=unlist(strsplit(treat,","));
#pairs=unlist(strsplit(pair,","));

control.index=match(controls,samples);
treat.index=match(treats,samples);

data.new=data[,c(control.index,treat.index)];
data.new=as.matrix(data.new);

group=c(rep("control",length(control.index)),rep("treat",length(treat.index)));
colData<-data.frame(row.names=c(controls,treats),sample=factor(group));
dds<-DESeqDataSetFromMatrix(countData=data.new, colData=colData, design=~sample);
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)
#dds <- DESeq(dds,fitType="local");
res <- results(dds);
sig<-res[which(res$padj<0.05),]
head=(sig)
print(head)
res<-res[order(res$padj),]
res<-res[,c(2,5,6)]
write.table(res,outfile,sep = "\t",quote=F);

#-----------------------------------------------------------------------#
#options <- commandArgs(trailingOnly = T)
infile = "Functional_analysis_summary.mdf.txt";
control = "CymD.D.LRP,CymD.D.MRP,CymS.D.HRP.1,CymS.D.LRP.1";
treat = "LBAb.D.HRP,LBAb.D.LRP,LBKM.D.HRP,LBKM.D.LRP,LBKM.D.MRP,LBKM.D.NBP,LBTa.D.RP,LBMa.D.RP,LBTu.D.RP";
output = "Functional_analysis_summary.P.Cym-v-Lau";
#pair="1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8";
outfile =paste(output,".diff.xls",sep="");

library("DESeq2");

data=read.table(infile,sep="\t",header=T,row.names=1);
data=data[rowSums(data) != 0,];
#colnames(data)<-gsub("X","",colnames(data))
samples=names(data)

controls=unlist(strsplit(control,","))
treats=unlist(strsplit(treat,","));
#pairs=unlist(strsplit(pair,","));

control.index=match(controls,samples);
treat.index=match(treats,samples);

data.new=data[,c(control.index,treat.index)];
data.new=as.matrix(data.new);

group=c(rep("control",length(control.index)),rep("treat",length(treat.index)));
colData<-data.frame(row.names=c(controls,treats),sample=factor(group));
dds<-DESeqDataSetFromMatrix(countData=data.new, colData=colData, design=~sample);
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)
#dds <- DESeq(dds,fitType="local");
res <- results(dds);
sig<-res[which(res$padj<0.05),]
head=(sig)
print(head)
res<-res[order(res$padj),]
res<-res[,c(2,5,6)]
write.table(res,outfile,sep = "\t",quote=F);

#-----------------------------------------------------------------------#
#options <- commandArgs(trailingOnly = T)
infile = "Functional_analysis_summary.mdf.txt";
control = "CymD.D.LRP,CymD.D.MRP,CymS.D.HRP.1,CymS.D.LRP.1";
treat = "GyBn.D.NBP";
output = "Functional_analysis_summary.P.Cym-v-GyBn";
#pair="1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8";
outfile =paste(output,".diff.xls",sep="");

library("DESeq2");

data=read.table(infile,sep="\t",header=T,row.names=1);
data=data[rowSums(data) != 0,];
#colnames(data)<-gsub("X","",colnames(data))
samples=names(data)

controls=unlist(strsplit(control,","))
treats=unlist(strsplit(treat,","));
#pairs=unlist(strsplit(pair,","));

control.index=match(controls,samples);
treat.index=match(treats,samples);

data.new=data[,c(control.index,treat.index)];
data.new=as.matrix(data.new);

group=c(rep("control",length(control.index)),rep("treat",length(treat.index)));
colData<-data.frame(row.names=c(controls,treats),sample=factor(group));
dds<-DESeqDataSetFromMatrix(countData=data.new, colData=colData, design=~sample);
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)
#dds <- DESeq(dds,fitType="local");
res <- results(dds);
sig<-res[which(res$padj<0.05),]
head=(sig)
print(head)
res<-res[order(res$padj),]
res<-res[,c(2,5,6)]
write.table(res,outfile,sep = "\t",quote=F);

#-----------------------------------------------------------------------#
#options <- commandArgs(trailingOnly = T)
infile = "Functional_analysis_summary.mdf.txt";
control = "GyBn.D.NBP";
treat = "LBAb.D.HRP,LBAb.D.LRP,LBKM.D.HRP,LBKM.D.LRP,LBKM.D.MRP,LBKM.D.NBP,LBTa.D.RP,LBMa.D.RP,LBTu.D.RP";
output = "Functional_analysis_summary.P.GyBn-v-Lau";
#pair="1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8";
outfile =paste(output,".diff.xls",sep="");

library("DESeq2");

data=read.table(infile,sep="\t",header=T,row.names=1);
data=data[rowSums(data) != 0,];
#colnames(data)<-gsub("X","",colnames(data))
samples=names(data)

controls=unlist(strsplit(control,","))
treats=unlist(strsplit(treat,","));
#pairs=unlist(strsplit(pair,","));

control.index=match(controls,samples);
treat.index=match(treats,samples);

data.new=data[,c(control.index,treat.index)];
data.new=as.matrix(data.new);

group=c(rep("control",length(control.index)),rep("treat",length(treat.index)));
colData<-data.frame(row.names=c(controls,treats),sample=factor(group));
dds<-DESeqDataSetFromMatrix(countData=data.new, colData=colData, design=~sample);
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)
#dds <- DESeq(dds,fitType="local");
res <- results(dds);
sig<-res[which(res$padj<0.05),]
head=(sig)
print(head)
res<-res[order(res$padj),]
res<-res[,c(2,5,6)]
write.table(res,outfile,sep = "\t",quote=F);

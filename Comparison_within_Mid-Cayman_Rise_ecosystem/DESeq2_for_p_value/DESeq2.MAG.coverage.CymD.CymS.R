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
infile = "CymD.10times.transpose.MAG.MetaG.coverage.txt";
control = "CymD.D.LB";
treat = "CymD.D.LRP,CymD.D.MRP";
output = "Background-v-Plume.CymD.MAG.MetaG.HydrothermalPlume.coverage";
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


#------------------------------------------------
infile = "CymS.10times.transpose.MAG.MetaG.coverage.txt";
control = "CymS.D.B";
treat = "CymS.D.HRP.1,CymS.D.LRP.1";
output = "Background-v-Plume.CymS.MAG.MetaG.HydrothermalPlume.coverage";
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

#------------------------------------------------

infile = "CymD.10times.transpose.MAG.MetaG.coverage.txt";
control = "CymD.D.LB";
treat = "CymD.D.PD.FS851,CymD.D.PD.FS852,CymD.D.PD.FS854,CymD.D.PD.FS856";
output = "Background-v-Plume.CymD.MAG.MetaG.VentFluid.coverage";
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


#----------------------------------------------------------
infile = "CymS.10times.transpose.MAG.MetaG.coverage.txt";
control = "CymS.D.B";
treat = "CymS.D.VD.FS841,CymS.D.VD.FS842,CymS.D.VD.FS844,CymS.D.VD.FS848,CymS.D.VD.FS849,CymS.D.VD.FS866,CymS.D.VD.FS872,CymS.D.VD.FS874,CymS.D.VD.FS877,CymS.D.VD.FS879,CymS.D.VD.FS881";
output = "Background-v-Plume.CymS.MAG.MetaG.VentFluid.coverage";
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
infile = "CymD.10times.transpose.MAG.MetaG.coverage.txt";
control = "CymD.D.LRP,CymD.D.MRP";
treat = "CymD.D.PD.FS851,CymD.D.PD.FS852,CymD.D.PD.FS854,CymD.D.PD.FS856";
output = "Background-v-Plume.CymD.MAG.MetaG.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymS.10times.transpose.MAG.MetaG.coverage.txt";
control = "CymS.D.HRP.1,CymS.D.LRP.1";
treat = "CymS.D.VD.FS841,CymS.D.VD.FS842,CymS.D.VD.FS844,CymS.D.VD.FS848,CymS.D.VD.FS849,CymS.D.VD.FS866,CymS.D.VD.FS872,CymS.D.VD.FS874,CymS.D.VD.FS877,CymS.D.VD.FS879,CymS.D.VD.FS881";
output = "Background-v-Plume.CymS.MAG.MetaG.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymD.100times.MAG.MetaT.coverage.txt";
control = "CymD.C.LB,CymD.C.HB";
treat = "CymD.C.LRP,CymD.C.MRP,CymD.C.HRP";
output = "Background-v-Plume.CymD.MAG.MetaT.HydrothermalPlume.coverage";
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
infile = "CymD.100times.MAG.MetaT.coverage.txt";
control = "CymD.C.LB,CymD.C.HB";
treat = "CymD.C.PD.FS851,CymD.C.PD.FS852,CymD.C.PD.FS854,CymD.C.PD.FS856";
output = "Background-v-Plume.CymD.MAG.MetaT.VentFluid.coverage";
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
infile = "CymD.100times.MAG.MetaT.coverage.txt";
control = "CymD.C.LRP,CymD.C.MRP,CymD.C.HRP";
treat = "CymD.C.PD.FS851,CymD.C.PD.FS852,CymD.C.PD.FS854,CymD.C.PD.FS856";
output = "Background-v-Plume.CymD.MAG.MetaT.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymS.100times.MAG.MetaT.coverage.txt";
control = "CymS.C.B,CymS.C.NBB.1,CymS.C.NBB.2";
treat = "CymS.C.LRP.1,CymS.C.LRP.2,CymS.C.HRP.1,CymS.C.HRP.2";
output = "Background-v-Plume.CymS.MAG.MetaT.HydrothermalPlume.coverage";
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
infile = "CymS.100times.MAG.MetaT.coverage.txt";
control = "CymS.C.B,CymS.C.NBB.1,CymS.C.NBB.2";
treat = "CymS.C.VD.FS841,CymS.C.VD.FS844,CymS.C.VD.FS848,CymS.C.VD.FS872,CymS.C.VD.FS879,CymS.C.VD.FS881";
output = "Background-v-Plume.CymS.MAG.MetaT.VentFluid.coverage";
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
infile = "CymS.100times.MAG.MetaT.coverage.txt";
control = "CymS.C.LRP.1,CymS.C.LRP.2,CymS.C.HRP.1,CymS.C.HRP.2";
treat = "CymS.C.VD.FS841,CymS.C.VD.FS844,CymS.C.VD.FS848,CymS.C.VD.FS872,CymS.C.VD.FS879,CymS.C.VD.FS881";
output = "Background-v-Plume.CymS.MAG.MetaT.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymD.KEGGlevel3.MetaG.abundance.txt";
control = "CymD.D.LB";
treat = "CymD.D.LRP,CymD.D.MRP";
output = "Background-v-Plume.CymD.KEGGlevel3.MetaG.HydrothermalPlume.coverage";
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


#------------------------------------------------
infile = "CymS.KEGGlevel3.MetaG.abundance.txt";
control = "CymS.D.B";
treat = "CymS.D.HRP.1,CymS.D.LRP.1";
output = "Background-v-Plume.CymS.KEGGlevel3.MetaG.HydrothermalPlume.coverage";
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

#------------------------------------------------

#------------------------------------------------

infile = "CymD.KEGGlevel3.MetaG.abundance.txt";
control = "CymD.D.LB";
treat = "CymD.D.PD.FS851,CymD.D.PD.FS852,CymD.D.PD.FS854,CymD.D.PD.FS856";
output = "Background-v-Plume.CymD.KEGGlevel3.MetaG.VentFluid.coverage";
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


#----------------------------------------------------------
infile = "CymS.KEGGlevel3.MetaG.abundance.txt";
control = "CymS.D.B";
treat = "CymS.D.VD.FS841,CymS.D.VD.FS842,CymS.D.VD.FS844,CymS.D.VD.FS848,CymS.D.VD.FS849,CymS.D.VD.FS866,CymS.D.VD.FS872,CymS.D.VD.FS874,CymS.D.VD.FS877,CymS.D.VD.FS879,CymS.D.VD.FS881";
output = "Background-v-Plume.CymS.KEGGlevel3.MetaG.VentFluid.coverage";
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
infile = "CymD.KEGGlevel3.MetaG.abundance.txt";
control = "CymD.D.LRP,CymD.D.MRP";
treat = "CymD.D.PD.FS851,CymD.D.PD.FS852,CymD.D.PD.FS854,CymD.D.PD.FS856";
output = "Background-v-Plume.CymD.KEGGlevel3.MetaG.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymS.KEGGlevel3.MetaG.abundance.txt";
control = "CymS.D.HRP.1,CymS.D.LRP.1";
treat = "CymS.D.VD.FS841,CymS.D.VD.FS842,CymS.D.VD.FS844,CymS.D.VD.FS848,CymS.D.VD.FS849,CymS.D.VD.FS866,CymS.D.VD.FS872,CymS.D.VD.FS874,CymS.D.VD.FS877,CymS.D.VD.FS879,CymS.D.VD.FS881";
output = "Background-v-Plume.CymS.KEGGlevel3.MetaG.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymD.Fun2MetaG.txt";
control = "CymD.D.LB";
treat = "CymD.D.LRP,CymD.D.MRP";
output = "Background-v-Plume.CymD.Fun.MetaG.HydrothermalPlume.coverage";
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


#------------------------------------------------
infile = "CymS.Fun2MetaG.txt";
control = "CymS.D.B";
treat = "CymS.D.HRP.1,CymS.D.LRP.1";
output = "Background-v-Plume.CymS.Fun.MetaG.HydrothermalPlume.coverage";
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

#------------------------------------------------

infile = "CymD.Fun2MetaG.txt";
control = "CymD.D.LB";
treat = "CymD.D.PD.FS851,CymD.D.PD.FS852,CymD.D.PD.FS854,CymD.D.PD.FS856";
output = "Background-v-Plume.CymD.Fun.MetaG.VentFluid.coverage";
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


#----------------------------------------------------------
infile = "CymS.Fun2MetaG.txt";
control = "CymS.D.B";
treat = "CymS.D.VD.FS841,CymS.D.VD.FS842,CymS.D.VD.FS844,CymS.D.VD.FS848,CymS.D.VD.FS849,CymS.D.VD.FS866,CymS.D.VD.FS872,CymS.D.VD.FS874,CymS.D.VD.FS877,CymS.D.VD.FS879,CymS.D.VD.FS881";
output = "Background-v-Plume.CymS.Fun.MetaG.VentFluid.coverage";
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
infile = "CymD.Fun2MetaG.txt";
control = "CymD.D.LRP,CymD.D.MRP";
treat = "CymD.D.PD.FS851,CymD.D.PD.FS852,CymD.D.PD.FS854,CymD.D.PD.FS856";
output = "Background-v-Plume.CymD.Fun.MetaG.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymS.Fun2MetaG.txt";
control = "CymS.D.HRP.1,CymS.D.LRP.1";
treat = "CymS.D.VD.FS841,CymS.D.VD.FS842,CymS.D.VD.FS844,CymS.D.VD.FS848,CymS.D.VD.FS849,CymS.D.VD.FS866,CymS.D.VD.FS872,CymS.D.VD.FS874,CymS.D.VD.FS877,CymS.D.VD.FS879,CymS.D.VD.FS881";
output = "Background-v-Plume.CymS.Fun.MetaG.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymD.KEGGlevel3.MetaT.abundance.txt";
control = "CymD.C.LB,CymD.C.HB";
treat = "CymD.C.LRP,CymD.C.MRP,CymD.C.HRP";
output = "Background-v-Plume.CymD.KEGGlevel3.MetaT.HydrothermalPlume.coverage";
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


#------------------------------------------------
infile = "CymS.KEGGlevel3.MetaT.abundance.txt";
control = "CymS.C.B,CymS.C.NBB.1,CymS.C.NBB.2";
treat = "CymS.C.LRP.1,CymS.C.LRP.2,CymS.C.HRP.1,CymS.C.HRP.2";
output = "Background-v-Plume.CymS.KEGGlevel3.MetaT.HydrothermalPlume.coverage";
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

#------------------------------------------------

#------------------------------------------------

infile = "CymD.KEGGlevel3.MetaT.abundance.txt";
control = "CymD.C.LB,CymD.C.HB";
treat = "CymD.C.PD.FS851,CymD.C.PD.FS852,CymD.C.PD.FS854,CymD.C.PD.FS856";
output = "Background-v-Plume.CymD.KEGGlevel3.MetaT.VentFluid.coverage";
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


#----------------------------------------------------------
infile = "CymS.KEGGlevel3.MetaT.abundance.txt";
control = "CymS.C.B,CymS.C.NBB.1,CymS.C.NBB.2";
treat = "CymS.C.VD.FS841,CymS.C.VD.FS844,CymS.C.VD.FS848,CymS.C.VD.FS872,CymS.C.VD.FS879,CymS.C.VD.FS881";
output = "Background-v-Plume.CymS.KEGGlevel3.MetaT.VentFluid.coverage";
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
infile = "CymD.KEGGlevel3.MetaT.abundance.txt";
control = "CymD.C.LRP,CymD.C.MRP,CymD.C.HRP";
treat = "CymD.C.PD.FS851,CymD.C.PD.FS852,CymD.C.PD.FS854,CymD.C.PD.FS856";
output = "Background-v-Plume.CymD.KEGGlevel3.MetaT.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymS.KEGGlevel3.MetaT.abundance.txt";
control = "CymS.C.LRP.1,CymS.C.LRP.2,CymS.C.HRP.1,CymS.C.HRP.2";
treat = "CymS.C.VD.FS841,CymS.C.VD.FS844,CymS.C.VD.FS848,CymS.C.VD.FS872,CymS.C.VD.FS879,CymS.C.VD.FS881";
output = "Background-v-Plume.CymS.KEGGlevel3.MetaT.HydrothermalPlume.vs.VentFluid.coverage";
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
infile = "CymD.Fun2MetaT.txt";
control = "CymD.C.LB,CymD.C.HB";
treat = "CymD.C.LRP,CymD.C.MRP,CymD.C.HRP";
output = "Background-v-Plume.CymD.Fun.MetaT.HydrothermalPlume.coverage";
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


#------------------------------------------------
infile = "CymS.Fun2MetaT.txt";
control = "CymS.C.B,CymS.C.NBB.1,CymS.C.NBB.2";
treat = "CymS.C.LRP.1,CymS.C.LRP.2,CymS.C.HRP.1,CymS.C.HRP.2";
output = "Background-v-Plume.CymS.Fun.MetaT.HydrothermalPlume.coverage";
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

#------------------------------------------------
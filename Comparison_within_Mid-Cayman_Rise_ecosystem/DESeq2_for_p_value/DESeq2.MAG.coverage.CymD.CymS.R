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
write.table(res,outfile,sep = "\t",quote=F)
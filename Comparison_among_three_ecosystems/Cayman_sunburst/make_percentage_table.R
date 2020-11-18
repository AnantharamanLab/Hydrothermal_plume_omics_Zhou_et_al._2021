data1<-read.table("MAG_average_coverage.Group2Row_Mean_MetaG.xls",sep = "\t",head=T, row.names = 1)
data2<-as.matrix(data1)
data3<-prop.table(data2,2)
write.table(data3,"MAG_average_coverage.Group2Row_Mean_MetaG.percentage.xls",sep = "\t",quote=F);
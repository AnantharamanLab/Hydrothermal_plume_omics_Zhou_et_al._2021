library("pheatmap")

x<-read.table("CymD.MAG.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt",sep="\t",head=T,row.names=1)
pdf(file="CymD.MAG.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.pdf",width=25,height=25)
pheatmap(x,scale="row",color =colorRampPalette(c("blue", "yellow"))(12),cluster_cols = F,show_colnames = T,cluster_row=F,border_color="NA", fontsize = 25)
dev.off()

library("pheatmap")

x<-read.table("CymS.MAG.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt",sep="\t",head=T,row.names=1)
pdf(file="CymS.MAG.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.pdf",width=25,height=25)
pheatmap(x,scale="row",color =colorRampPalette(c("blue", "yellow"))(12),cluster_cols = F,show_colnames = T,cluster_row=F,border_color="NA", fontsize = 25)
dev.off()


library("pheatmap")

x<-read.table("CymD.MAG.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt",sep="\t",head=T,row.names=1)
pdf(file="CymD.MAG.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.pdf",width=25,height=25)
pheatmap(x,scale="row",color =colorRampPalette(c("blue", "yellow"))(12),cluster_cols = F,show_colnames = T,cluster_row=F,border_color="NA", fontsize = 25)
dev.off()


library("pheatmap")

x<-read.table("CymS.MAG.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt",sep="\t",head=T,row.names=1)
pdf(file="CymS.MAG.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.pdf",width=25,height=25)
pheatmap(x,scale="row",color =colorRampPalette(c("blue", "yellow"))(12),cluster_cols = F,show_colnames = T,cluster_row=F,border_color="NA", fontsize = 25)
dev.off()

library("pheatmap")

x<-read.table("CymD.Fun.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt",sep="\t",head=T,row.names=1)
pdf(file="CymD.Fun.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.pdf",width=25,height=25)
pheatmap(x,scale="row",color =colorRampPalette(c("blue", "yellow"))(12),cluster_cols = F,show_colnames = T,cluster_row=F,border_color="NA", fontsize = 25)
dev.off()

library("pheatmap")

x<-read.table("CymS.Fun.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt",sep="\t",head=T,row.names=1)
pdf(file="CymS.Fun.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.pdf",width=25,height=25)
pheatmap(x,scale="row",color =colorRampPalette(c("blue", "yellow"))(12),cluster_cols = F,show_colnames = T,cluster_row=F,border_color="NA", fontsize = 25)
dev.off()

library("pheatmap")

x<-read.table("CymD.Fun.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt",sep="\t",head=T,row.names=1)
pdf(file="CymD.Fun.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.pdf",width=25,height=25)
pheatmap(x,scale="row",color =colorRampPalette(c("blue", "yellow"))(12),cluster_cols = F,show_colnames = T,cluster_row=F,border_color="NA", fontsize = 25)
dev.off()

library("pheatmap")

x<-read.table("CymS.Fun.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt",sep="\t",head=T,row.names=1)
pdf(file="CymS.Fun.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.pdf",width=25,height=25)
pheatmap(x,scale="row",color =colorRampPalette(c("blue", "yellow"))(12),cluster_cols = F,show_colnames = T,cluster_row=F,border_color="NA", fontsize = 25)
dev.off()



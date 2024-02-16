setwd("~/projects/neoags/")
load("preProcess.RData")
load("xysWithAnnotations")
options(stringsAsFactors = F)
library(pheatmap)
library(imcRtools)
library(Rtsne)

l.markers.orig <- c("X141Pr_CD14.ome","X142Nd_FoxP3.ome","X143Nd_epcam.ome","X145Nd_CD4.ome","X146Nd_CD8a.ome",
                    "X147Sm_CollagenI.ome","X154Sm_CD163.ome","X159Tb_CD68.ome","X160Gd_CD11b.ome","X162Dy_CD11c.ome",
                    "X161Dy_CD20.ome","X170Er_CD3e.ome","X174Yb_CD57.ome","X176Yb_Pan.cytokeratin.ome","X194Pt_aSMA.ome",
                    "X198Pt_Vimentin.ome","X89Y_CD45.ome")
name.tbl <- data.frame(ref=colnames(xys)[5:43],
                       alt=c("IL6","CD14","Foxp3","EPCAM","HLA-I","CD4","CD8a","CollagenI","B2M",
                             "HLA-II","TIM3","PSMB8","PSME1","TIGIT","CD163","IDO1","PD-L1","LAG3",
                             "CD68","CD11b","CD20","CD11c","CD40","GranzymeB","PD-1","Ki67","CTSL&CTSB",
                             "PSMD7","TAP1&2","CD3e","TNFα","STAT3","SNAIL&SLUG","CD57","SMAD2",
                             "Pan-Cytokeratin","αSMA","Vimentin","CD45"))
l.markers <- c("CD14","Foxp3","EPCAM","CD4","CD8a","CollagenI","CD163","CD68","CD11b","CD11c",
               "CD20","CD3e","CD57","Pan-Cytokeratin","αSMA","Vimentin","CD45")

exprs.m <- aggregate(xys[,l.markers.orig], by=list(xys$Phenograph), mean)
rownames(exprs.m) <- exprs.m$Group.1
exprs.m$Group.1 <- NULL
colnames(exprs.m) <- name.tbl[match(colnames(exprs.m), name.tbl$ref),"alt"]

grDevices::cairo_pdf("figures/heatmap.pdf", height=5, width=10, onefile = F)
pheatmap(exprs.m, border_color = "white", scale = "none",color = paletteer_c("grDevices::Geyser", 30))
dev.off()

roi.total <- dplyr::count(xys, ROI)
ct.prop <- dplyr::count(xys, ROI, Phenograph)
ct.prop$total <- roi.total[match(ct.prop$ROI, roi.total$ROI),"n"] 
ct.prop$freq <- ct.prop$n/ct.prop$total
ct.prop.w <- reshape(ct.prop[,c(1,2,5)], idvar = "ROI", timevar = "Phenograph", direction = "wide")
ct.prop.w[is.na(ct.prop.w)] <- 0
ct.prop <- reshape2::melt(ct.prop.w)
ct.prop$variable <- gsub("freq\\.","",ct.prop$variable)

hm <- pheatmap(exprs.m,border_color = "white", scale = "none")
hm.orders <- rownames(exprs.m)[hm$tree_row$order]
ct.prop.draw <- ct.prop
ct.prop.draw$variable <- factor(ct.prop.draw$variable, levels = hm.orders)
ggplot(ct.prop.draw, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill=variable),color="white") + 
  geom_jitter(aes(fill=variable),color="white", pch=21,stroke=0.1) + theme_classic()+
  scale_fill_manual(values = ct.cols)+
  theme(legend.position = "none",axis.text.x = element_blank())+
  ylab("Frequency")+xlab("")
ggsave("figures/ct-boxplot-all.pdf", width = 6, height = 3)

set.seed(42)
sample.idx <- data.frame(id=seq(1,nrow(xys)), Phenograph=xys$Phenograph) %>%
  group_by(Phenograph) %>% dplyr::slice_sample(prop=.1)
tsne_out <- Rtsne(xys[sample.idx$id,l.markers.orig], check_duplicates=FALSE,set.seed=42)

ct.order <- c("Epithelials","DedifferentiatedEpithelials","Epithelials/Fibroblasts","Fibroblasts",
              "ResidentMφ","InfiltratingMφ","Epithelials/ResidentMφ/Treg","NK",
              "CD4+TCells","CD8+TCells",'VIM+TCells',"BCells","Lineage-")
ct.cols <- setNames(c(paletteer_d("ggthemes::Classic_Green_Orange_12"),
                      paletteer_d("ggthemes::Classic_Blue_Red_12")[1]),ct.order)
ct.colsWithNum <- ct.cols
names(ct.colsWithNum) <- paste0(seq(1,13),":",names(ct.cols))

dat <- data.frame(tsne_out$Y)
dat$Phenograph <- xys$Phenograph[sample.idx$id]
colnames(dat)[1:2] <- c("TSNE1","TSNE2")
dat$Phenograph <- paste0(match(dat$Phenograph,ct.order),":",dat$Phenograph)
dat$Phenograph <- factor(dat$Phenograph,levels = names(ct.colsWithNum))
pos.df <- aggregate(dat[,c("TSNE1","TSNE2")], median, by=list(dat$Phenograph))
p <- ggplot(data = dat, aes(TSNE1,TSNE2))+geom_point(aes(color=Phenograph))+ 
  scale_color_manual(values = ct.colsWithNum, name="")+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  geom_text(data=pos.df, aes(x=TSNE1, y=TSNE2, label=gsub("(.*):.*","\\1",Group.1)),size=10)+
  theme_void()+theme(legend.text=element_text(size=12, face = "bold"))
p
grDevices::cairo_pdf("figures/IMCTsne-202312.pdf", height = 8, width = 12)
print(p)
dev.off()

cbind(xys[sample.idx$id,5:43],dat) %>% select(-Phenograph) %>% 
  melt(.,id.vars=c("TSNE1","TSNE2")) -> dat.draw
p.list <- lapply(colnames(xys)[5:43], function(i) 
{ggplot(data = dat.draw[dat.draw$variable==i, ], aes(TSNE1,TSNE2))+geom_point(aes(color=value))+ 
    scale_color_paletteer_c("ggthemes::Classic Orange-Blue","",direction = -1)+
    # geom_text(data = pos.df, aes(x=TSNE1,y=TSNE2,label=Group.1))+
    theme_void()})
lapply(seq(5,43), function(i){
  ggsave(paste0("imc-figures/tsne-expr/", colnames(xys)[i],".tiff"), p.list[[i-4]], height = 5, width = 6.5)
})

files <- unique(xys$ROI)
spaList = foreach(file=files,.packages = c("imcRtools","doParallel"),.errorhandling="pass") %dopar% {
  dat <- xys[xys$ROI==file,]
  cd <- DataFrame(Pos_X = dat$X_position, Pos_Y = dat$Y_position, z=rownames(dat), 
                  ImageNb=file,CellType=dat$Phenograph)
  spe <- SpatialExperiment(assay = t(dat[,5:43]), colData = cd, sample_id=file,
                           spatialDataNames = c("z"), spatialCoordsNames = c("Pos_X", "Pos_Y"))
  spe[["Pos_X"]] <- dat$X_position
  spe[["Pos_Y"]] <- dat$Y_position
  
  spe <- buildSpatialGraph(spe, img_id = "ImageNb",type = "knn",k = 10)
  nbrs <- data.frame(colPair(spe, "knn_interaction_graph"))
  nbrs <- apply(nbrs, 2, function(x)rownames(dat)[x])
  
  spe <- buildSpatialGraph(spe, img_id = "ImageNb",type = "expansion",threshold = 20)
  cts <- unique(xys$Phenograph)
  patchList <- foreach(ct=cts,.packages = "imcRtools",.errorhandling="remove", n = cts, .combine = c) %dopar% {
    tryCatch({
      rl <- list()
      rl[[n]] <- patchDetection(spe, patch_cells = spe$CellType == ct,
                                colPairName = "expansion_interaction_graph",
                                expand_by = 20, min_patch_size=50,
                                img_id = "ImageNb")
      rl
    })
  }
  return(list(nbrs,patchList))
}

load("trait")
sample.info <- read.table("sample-info.csv", header = F, sep=",")
sample.info <- dplyr::left_join(sample.info,trait,by=c("V2"="sample"))
sample.info$tumor_stage.diagnoses[sample.info$tumor_stage.diagnoses=="I"] <- "I/II"
sample.info$tumor_stage.diagnoses[sample.info$tumor_stage.diagnoses=="II"] <- "I/II"
sample.info$tumor_stage.diagnoses[sample.info$tumor_stage.diagnoses=="III"] <- "III/IV"
sample.info$tumor_stage.diagnoses[sample.info$tumor_stage.diagnoses=="IV"] <- "III/IV"
colnames(sample.info)
dplyr::count(sample.info,tumor_stage.diagnoses,microsatellite_instability)
sample.info$microsatellite_instability <- factor(sample.info$microsatellite_instability, levels = c("MSI-H","MSS"))
sample.info$tumor_stage.diagnoses <- factor(sample.info$tumor_stage.diagnoses, levels = c("I/II","III/IV"))

library(paletteer)
library(RColorBrewer)
library(ggsci)
library(cowplot)

col=list(msi=setNames(pal_d3()(2),c("MSS","MSI-H")),
         cms=setNames(pal_npg()(4),c("CMS1","CMS2","CMS3","CMS4")),
         stage=setNames(as.character(paletteer_d("ggsci::deep_purple_material"))[c(2,7)],c("I/II","III/IV")))

ct.d <- inner_join(sample.info, ct.prop.w, by=c("V1"="ROI"))
colnames(ct.d) <- gsub("freq\\.","",colnames(ct.d))
roi.order <- ct.d %>% arrange(microsatellite_instability, tumor_stage.diagnoses,desc(Epithelials)) %>% select(V1) 
ct.prop.draw$ROI <- factor(ct.prop.draw$ROI, levels = roi.order$V1)
ct.prop.draw$variable <- factor(ct.prop.draw$variable, levels = ct.order)
p1 <- ggplot(ct.prop.draw, aes(x=ROI, y=value, fill=variable))+geom_bar(position = "fill",stat = "identity")+
  theme_classic()+ylab("Frequency")+xlab("")+
  scale_fill_manual(values = ct.cols)+
  theme(axis.text.x = element_blank())
p1
sample.info$V1 <- factor(sample.info$V1, levels = roi.order$V1)
ms.p <- ggplot(sample.info[,c("V1","microsatellite_instability")], aes(x=V1,y=1))+
  geom_tile(aes(fill=microsatellite_instability), color="white")+
  scale_fill_manual(values = col$msi)+
  theme_classic()+xlab("")+ylab("")+theme(axis.text.x = element_blank(), 
                                          axis.ticks.x = element_blank(),
                                          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                          axis.line.x = element_blank(), axis.line.y = element_blank())

ss.p <- ggplot(sample.info[,c("V1","tumor_stage.diagnoses")], aes(x=V1,y=1))+
  geom_tile(aes(fill=tumor_stage.diagnoses), color="white")+
  scale_fill_manual(values = col$stage)+
  theme_classic()+xlab("")+ylab("")+theme(axis.text.x = element_blank(), 
                                          axis.ticks.x = element_blank(),
                                          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                          axis.line.x = element_blank(), axis.line.y = element_blank())
cms.p <- ggplot(sample.info[,c("V1","cms")], aes(x=V1,y=1))+
  geom_tile(aes(fill=cms), color="white")+
  scale_fill_manual(values = col$cms)+
  theme_classic()+xlab("")+ylab("")+theme(axis.text.x = element_blank(), 
                                          axis.ticks.x = element_blank(),
                                          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                          axis.line.x = element_blank(), axis.line.y = element_blank())

ggsave("figures/ctprops.pdf", height = 5, width = 8)
plot_grid(p1,ms.p, ss.p, cms.p,ncol=1, rel_heights = c(5,1,1,1))
dev.off()

sampleih.mhc1 <- read.table("~/projects/neoantigens/gastrointestinal/CRC-4.0/data/in-house/MHCI_results.gz", header = F)
new <- read.delim("~/projects/neoantigens/gastrointestinal/CRC-2020/pna-new/MHCI_results.tsv", header = T, sep="\t", fill=T)
colnames(sampleih.mhc1) <- colnames(new)
mhc1 <- rbind(sampleih.mhc1, new)
mhc1$id <- paste(mhc1$chromosome, mhc1$start, mhc1$reference, mhc1$alteration, sep=":")
mhc1 <- mhc1[mhc1$sample %in% trait$sample,]
mhc1.count <- dplyr::count(unique(mhc1[,c("sample","id")]),sample)

sampleih.mhc2 <- read.table("~/projects/neoantigens/gastrointestinal/CRC-4.0/data/in-house/MHCII_results.gz", header = F)
new <- read.delim("~/projects/neoantigens/gastrointestinal/CRC-2020/pna-new/MHCII_results.tsv", header = T, sep="\t", fill=T)
colnames(sampleih.mhc2) <- colnames(new)
mhc2 <- rbind(sampleih.mhc2, new)
mhc2$id <- paste(mhc2$chromosome, mhc2$start, mhc2$reference, mhc2$alteration, sep=":")
mhc2 <- mhc2[mhc2$sample %in% trait$sample,]
mhc2.count <- dplyr::count(unique(mhc2[,c("sample","id")]),sample)

library(reshape2)
ct.prop$V2 <- sample.info[match(ct.prop$ROI, sample.info$V1),"V2"]
ct.prop$mhc1.count <- mhc1.count[match(ct.prop$V2, mhc1.count$sample),"n"]
ct.prop$mhc2.count <- mhc2.count[match(ct.prop$V2, mhc2.count$sample),"n"]
mhc1.cor <- as.data.frame(do.call(rbind,lapply(ct.order[1:12], function(ct){
  unlist(cor.test(ct.prop[ct.prop$variable==ct,"value"],ct.prop[ct.prop$variable==ct,"mhc1.count"])[3:4])
}))) 
mhc1.cor$ct <- ct.order[1:12]
mhc1.cor$sig <- ifelse(mhc1.cor$p.value<.05,"*","")
mhc2.cor <- as.data.frame(do.call(rbind,lapply(ct.order[1:12], function(ct){
  unlist(cor.test(ct.prop[ct.prop$variable==ct,"value"],ct.prop[ct.prop$variable==ct,"mhc2.count"])[3:4])
}))) 
mhc2.cor$ct <- ct.order[1:12]
mhc2.cor$sig <- ifelse(mhc2.cor$p.value<.05,"*","")

cor.df <- inner_join(mhc1.cor, mhc2.cor, by="ct")
colnames(cor.df)[c(2,6)] <- c("MHCI","MHCII")
cor.df.m <- melt(cor.df[c(3,2,6)])
cor.df.m$ct <- factor(cor.df.m$ct, levels = ct.order[1:12])
p <- ggplot() + 
  geom_tile(data=cor.df.m,aes(ct,variable), color = "grey90",lwd = .5,linetype = 1, fill=NA)+
  geom_point(data=cor.df.m,aes(ct,variable,size=abs(value), color=value))+
  
  geom_text(data=cor.df[cor.df$sig.x=="*",], aes(x=ct,y=1,label=sig.x),size=5,vjust=.75)+
  geom_text(data=cor.df[cor.df$sig.y=="*",], aes(x=ct,y=2,label=sig.y),size=5)+
  scale_color_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-.31,.31),na.value = NA,
                          name="Correlation\nCoefficients")+
  scale_size(limits=c(0,.31),guide="none")+
    theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.margin=margin(t=80))+
  xlab("")+ylab("")

grDevices::cairo_pdf("figures/NeoAgsWithCTProportion.pdf", width = 8, height = 2.5)
print(p)
dev.off()

f.markers <- setdiff(colnames(xys), l.markers.orig)[6:22]
f.order <- c("X144Nd_HLA.I.ome","X149Sm_HLA.II.ome",
             "X148Nd_beta2m.ome","X151Eu_psmb8.ome","X168Er_PSMDs.ome","X152Sm_psme1.ome",
             "X167Er_CTSL_CTSB.ome","X169Tm_TAP1_2.ome",
             "X150Nd_HAVCR2.ome","X153Eu_TIGIT.ome","X155Gd_IDO1.ome","X156Gd_PD.L1.ome","X158Gd_lag3.ome","X163Dy_CD40.ome","X165Ho_PD.1.ome",
             "X166Er_Ki.67.ome","X164Dy_GranzymeB.ome")
f.order <- name.tbl[match(f.order, name.tbl$ref),"alt"]
exprs <- xys
exprs <- dplyr::left_join(exprs, sample.info, by=c("ROI"="V1"))
exprs <- exprs[exprs$Phenograph!="Lineage-",]

stats <- exprs[,c(f.markers,"ROI","Phenograph")] %>% group_by(ROI,Phenograph) %>% dplyr::select(-ROI,-Phenograph) %>%
  summarise_all(mean) %>% melt
stats$Phenograph <- factor(stats$Phenograph, levels = ct.order[1:12])
stats$variable <- name.tbl[match(stats$variable,name.tbl$ref),"alt"]
stats$variable <- factor(stats$variable, levels = f.order)
stats <- stats %>% filter(!variable %in% c("Ki67","GranzymB"))
stats <- stats %>% filter(variable != "GranzymeB")

stats.l <- exprs[exprs$tumor_stage.diagnoses=="I/II",c(f.markers,"ROI","Phenograph")] 
stats.h <- exprs[exprs$tumor_stage.diagnoses=="III/IV",c(f.markers,"ROI","Phenograph")]

s.stats <- data.frame(ct=as.character(),marker=as.character(),fc=as.numeric(), pval=as.numeric())
for (ct in ct.order[1:12]){
  for (marker in f.markers){
    test <- t.test(stats.h[stats.h$Phenograph==ct,marker], stats.l[stats.l$Phenograph==ct,marker])
    pval <- test$p.value
    fc <- test$estimate[1]/test$estimate[2]
    s.stats[nrow(s.stats)+1,] <- c(ct, marker, fc, pval)
  }
}

s.stats$fc <- log2(as.numeric(s.stats$fc))
s.stats$label <- ifelse(as.numeric(s.stats$pval)<0.05&abs(s.stats$fc)>1, "p<0.05&|log2(FoldChange)|>1",
                        "NotSignificant")
s.stats$ct <- factor(s.stats$ct, levels = ct.order[1:12])
s.stats$marker <- name.tbl[match(s.stats$marker,name.tbl$ref),"alt"]
s.stats$marker <- factor(s.stats$marker, levels = f.order)
s.stats <- s.stats %>% filter(!marker %in% c("Ki67","GranzymB"))
s.stats <- s.stats %>% filter(marker != "GranzymeB")

stats.mss <- exprs[exprs$microsatellite_instability=="MSS",c(f.markers,"ROI","Phenograph")]
stats.msi <- exprs[exprs$microsatellite_instability=="MSI-H",c(f.markers,"ROI","Phenograph")]

ms.stats <- data.frame(ct=as.character(),marker=as.character(),fc=as.numeric(), pval=as.numeric())
for (ct in ct.order[1:12]){
  for (marker in f.markers){
    test <- t.test(stats.msi[stats.msi$Phenograph==ct,marker], stats.mss[stats.mss$Phenograph==ct,marker])
    pval <- test$p.value
    fc <- test$estimate[1]/test$estimate[2]
    ms.stats[nrow(ms.stats)+1,] <- c(ct, marker, fc, pval)
  }
}
ms.stats$fc <- log2(as.numeric(ms.stats$fc))
ms.stats$label <- ifelse(as.numeric(ms.stats$pval)<0.05&abs(ms.stats$fc)>1, "p<0.05&|log2(FoldChange)|>1",
                        "NotSignificant")
ms.stats$ct <- factor(ms.stats$ct, levels = ct.order[1:12])
ms.stats$marker <- name.tbl[match(ms.stats$marker,name.tbl$ref),"alt"]
ms.stats$marker <- factor(ms.stats$marker, levels = f.order)
ms.stats <- ms.stats %>% filter(!marker %in% c("Ki67","GranzymB"))
ms.stats <- ms.stats %>% filter(marker != "GranzymeB")


p <- ggplot(data=stats, aes(x=Phenograph,y=variable))+geom_tile(aes(fill=value), color="grey80")+
  scale_fill_paletteer_c("grDevices::Purples", "Expression", direction = -1,limits=c(0,1))+
  geom_point(data=s.stats, aes(x=ct, y=marker,color=fc,size=label),position = position_nudge(x=-.2))+
  geom_point(data=ms.stats, aes(x=ct, y=marker,color=fc,size=label),position = position_nudge(x=.2))+
  scale_size_manual(values = c(1,5))+
  scale_color_gradientn(colours=as.character(rev(paletteer_d("RColorBrewer::RdYlBu"))),
                        values =c(0, seq(.3,.7, length.out=9), 1),
                        name="log2(FoldChange)", limits=c(-4.5,4.5))+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),
                                     legend.margin=margin(t=-20))
grDevices::cairo_pdf("figures/ExprsComparisons.pdf",height=5, width=9)
print(p)
dev.off()



exprs$mhc1 <- mhc1.count[match(exprs$V2, mhc1.count$sample),"n"]
exprs$mhc2 <- mhc2.count[match(exprs$V2, mhc2.count$sample),"n"]

mhc1.cor <- as.data.frame(do.call(rbind,lapply(f.markers, function(marker){
  setNames(as.data.frame(do.call(rbind,lapply(split(exprs, factor(exprs$Phenograph)), 
                                              function(x) {
                                                avrg <- x[,c("ROI",marker,"mhc1")] %>% group_by(ROI) %>% 
                                                  summarise_at(c(eval(marker),"mhc1"),mean) %>% as.data.frame()
                                                unlist(cor.test(avrg[,2],avrg[,3])[3:4])
                                              }))),
           c("p.value","estimate")) %>% rownames_to_column()
})))

mhc1.cor$marker <- rep(name.tbl[match(f.markers,name.tbl$ref),"alt"],each=12)

mhc2.cor <- as.data.frame(do.call(rbind,lapply(f.markers, function(marker){
  setNames(as.data.frame(do.call(rbind,lapply(split(exprs, factor(exprs$Phenograph)), 
                                              function(x) {
                                                avrg <- x[,c("ROI",marker,"mhc2")] %>% group_by(ROI) %>% 
                                                  summarise_at(c(eval(marker),"mhc2"),mean) %>% as.data.frame()
                                                unlist(cor.test(avrg[,2],avrg[,3])[3:4])
                                              }))),
           c("p.value","estimate")) %>% rownames_to_column()
})))
mhc2.cor$marker <- rep(name.tbl[match(f.markers,name.tbl$ref),"alt"],each=12)
count.cor <- inner_join(mhc1.cor, mhc2.cor,by=c("rowname","marker"))
count.cor$rowname <- factor(count.cor$rowname, levels = ct.order[1:12])
count.cor$marker <- factor(count.cor$marker, levels = f.order)
count.cor <- count.cor %>% filter(!marker %in% c("GranzymeB","Ki67"))

p <- ggplot(count.cor,aes(x = rowname,y = marker)) +
  geom_jjtriangle(aes(fill = estimate.x),type = 'ul',color="white",size=.5) +
  scale_fill_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),
                       limits=c(-1,1),name="Correlation\nCoefficients") +
  new_scale_fill() +
  geom_jjtriangle(aes(fill = estimate.y),type = 'br', color="white",size=.5) +
  scale_fill_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-1,1),guide="none") +
  geom_text(data=count.cor[count.cor$p.value.x<.05,],
            aes(x = rowname,y = marker),label="*",position = position_nudge(x=-.2,y=.1),size=5)+
  geom_text(data=count.cor[count.cor$p.value.y<.05,],
            aes(x = rowname,y = marker),label="*",position = position_nudge(x=.2,y=-.3),size=5)+
  theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab("")+coord_fixed()
grDevices::cairo_pdf("figures/NeoAgsWithExpression.pdf", height=7, width=6)
print(p)
dev.off()


nbrs <- as.data.frame(do.call(rbind, lapply(spaList, function(x)x[[1]])))
nbrs$Phenograph <- xys[match(nbrs$to, rownames(xys)),"Phenograph"]
nbrs.df <- dplyr::count(nbrs, from, Phenograph)
nbrs.df <- reshape(nbrs.df, idvar = "from", timevar = "Phenograph", direction = "wide")
rownames(nbrs.df) <- nbrs.df$from
nbrs.df$from <- NULL
nbrs.df[is.na(nbrs.df)] <-0
nbrs.df$ROI <- gsub("(.*)-.*","\\1", rownames(nbrs.df))
nbrs.df$Phenograph <- xys[match(rownames(nbrs.df), rownames(xys)),"Phenograph"]
colnames(nbrs.df) <- gsub("n\\.", "",colnames(nbrs.df))

nbrs.df <- dplyr::left_join(nbrs.df, sample.info, by=c("ROI"="V1"))
nbrs.df <- nbrs.df[nbrs.df$Phenograph!="Lineage-",]
nbrs.df$`Lineage-` <- NULL

library(jjPlot)
library(dplyr)

stats <- nbrs.df[,c(1:12,14)] %>%
  group_by(Phenograph) %>% dplyr::select(-Phenograph) %>%
  summarise_all(mean) %>% melt
stats$Phenograph <- factor(stats$Phenograph, levels = ct.order[1:12])
stats$variable <- factor(stats$variable, levels = ct.order[1:12])

stats.l <- nbrs.df[nbrs.df$tumor_stage.diagnoses=="I/II",c(1:12,13,14)] %>%
  melt(id.vars=c("ROI","Phenograph"))
stats.h <- nbrs.df[nbrs.df$tumor_stage.diagnoses=="III/IV",c(1:12,13,14)] %>%
  melt(id.vars=c("ROI","Phenograph"))

s.stats <- data.frame(center=as.character(),neighbour=as.character(),fc=as.numeric(), pval=as.numeric())
for (ct in ct.order[1:12]){
  for (nb in ct.order[1:12]){
    test <- t.test(stats.h[stats.h$Phenograph==ct&stats.h$variable==nb,"value"], 
                   stats.l[stats.l$Phenograph==ct&stats.l$variable==nb,"value"])
    pval <- test$p.value
    fc <- test$estimate[1]/test$estimate[2]
    s.stats[nrow(s.stats)+1,] <- c(ct, nb, fc, pval)
  }
}
s.stats$fc <- as.numeric(s.stats$fc)
s.stats$fc <- log2(s.stats$fc)
s.stats$label <- ifelse(as.numeric(s.stats$pval)<0.05&abs(fc)>1, 
                        "p<0.05&|log2(FoldChange)|>1","NotSignificant")
s.stats$center <- factor(s.stats$center, levels = ct.order[1:12])
s.stats$neighbour <- factor(s.stats$neighbour, levels = ct.order[1:12])

stats.mss <- nbrs.df[nbrs.df$microsatellite_instability=="MSS",c(1:12,13,14)] %>% 
  melt(id.vars=c("ROI","Phenograph"))
stats.msi <- nbrs.df[nbrs.df$microsatellite_instability=="MSI-H",c(1:12,13,14)] %>%
  melt(id.vars=c("ROI","Phenograph"))

ms.stats <- data.frame(center=as.character(),neighbour=as.character(),fc=as.numeric(), pval=as.numeric())
for (ct in ct.order[1:12]){
  for (nb in ct.order[1:12]){
    test <- t.test(stats.msi[stats.msi$Phenograph==ct&stats.msi$variable==nb,"value"], 
                   stats.mss[stats.mss$Phenograph==ct&stats.mss$variable==nb,"value"])
    pval <- test$p.value
    fc <- test$estimate[1]/test$estimate[2]
    ms.stats[nrow(ms.stats)+1,] <- c(ct, nb, fc, pval)
  }
}
ms.stats$fc <- as.numeric(ms.stats$fc)
ms.stats$fc <- log2(ms.stats$fc)
ms.stats$label <- ifelse(as.numeric(ms.stats$pval)<0.05&abs(fc)>1, 
                        "p<0.05&|log2(FoldChange)|>1","NotSignificant")
ms.stats$center <- factor(ms.stats$center, levels = ct.order[1:12])
ms.stats$neighbour <- factor(ms.stats$neighbour, levels = ct.order[1:12])



p <- ggplot(data=stats, aes(x=variable,y=Phenograph))+geom_tile(aes(fill=value), color="grey80")+
  scale_fill_paletteer_c("grDevices::Purples", "CellCounts", direction = -1,limits=c(0,10))+
  geom_point(data=s.stats, aes(x=neighbour, y=center,color=fc,size=label),position = position_nudge(x=-.2))+
  geom_point(data=ms.stats, aes(x=neighbour, y=center,color=fc,size=label),position = position_nudge(x=.2))+
  scale_size_manual(values = c(1,5))+
  scale_color_gradientn(colours=as.character(rev(paletteer_d("RColorBrewer::RdYlBu"))),
                        values =c(0, seq(.3,.7, length.out=9), 1),limits=c(-6, 6),
                        name="log(FoldChange)")+
 theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-20))
grDevices::cairo_pdf("figures/CCIComparisons.pdf",height=5, width=10)
print(p)
dev.off()

library(jjPlot)
library(ggnewscale)

nbrs.df$mhc1 <- mhc1.count[match(nbrs.df$V2, mhc1.count$sample),"n"]
nbrs.df$mhc2 <- mhc2.count[match(nbrs.df$V2, mhc2.count$sample),"n"]

mhc1.cor <- as.data.frame(do.call(rbind,lapply(ct.order[1:12], function(ct){
  setNames(as.data.frame(do.call(rbind,lapply(split(nbrs.df, factor(nbrs.df$Phenograph)), 
                                              function(x) {
                                                avrg <- x[,c("ROI",ct,"mhc1")] %>% group_by(ROI) %>% 
                                                  summarise_at(c(eval(ct),"mhc1"),mean) %>% as.data.frame()
                                                unlist(cor.test(avrg[,2],avrg[,3])[3:4])
                                              }))),
           c("p.value","estimate")) %>% rownames_to_column()
})))

mhc1.cor$center <- rep(ct.order[1:12],each=12)
mhc2.cor <- as.data.frame(do.call(rbind,lapply(ct.order[1:12], function(ct){
  setNames(as.data.frame(do.call(rbind,lapply(split(nbrs.df, factor(nbrs.df$Phenograph)), 
                                              function(x) {
                                                avrg <- x[,c("ROI",ct,"mhc2")] %>% group_by(ROI) %>% 
                                                  summarise_at(c(eval(ct),"mhc2"),mean) %>% as.data.frame()
                                                unlist(cor.test(avrg[,2],avrg[,3])[3:4])
                                              }))),
           c("p.value","estimate")) %>% rownames_to_column()
})))
mhc2.cor$center <- rep(ct.order[1:12],each=12)
count.cor <- inner_join(mhc1.cor, mhc2.cor,by=c("rowname","center"))
count.cor$rowname <- factor(count.cor$rowname, levels = ct.order[1:12])
count.cor$center <- factor(count.cor$center, levels = ct.order[1:12])

p <- ggplot(count.cor,aes(x = rowname,y = center)) +
  geom_jjtriangle(aes(fill = estimate.x),type = 'ul',color="white",size=.5) +
  scale_fill_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),
                       limits=c(-1,1),name="Correlation\nCoefficients") +
  new_scale_fill() +
  geom_jjtriangle(aes(fill = estimate.y),type = 'br', color="white",size=.5) +
  scale_fill_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-1,1),guide="none") +
  geom_text(data=count.cor[count.cor$p.value.x<.05,],
            aes(x = rowname,y = center),label="*",
            position = position_nudge(x=-.2,y=.1),size=5)+
  geom_text(data=count.cor[count.cor$p.value.y<.05,],
            aes(x = rowname,y = center),label="*",
            position = position_nudge(x=.2,y=-.3),size=5)+
  theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("")+ylab("")+coord_fixed()

grDevices::cairo_pdf("figures/NeoAgsWithCCI.pdf", height=6, width=6)
print(p)
dev.off()

apps <- f.order[1:8]
icbs <- f.order[9:15]

app.exprs <- exprs[,c(name.tbl[match(apps,name.tbl$alt),"ref"],"ROI","Phenograph")] %>% group_by(ROI,Phenograph) %>%
  summarise_at(name.tbl[match(apps,name.tbl$alt),"ref"], mean)
icb.exprs <- exprs[,c(name.tbl[match(icbs,name.tbl$alt),"ref"],"ROI","Phenograph")] %>% group_by(ROI,Phenograph) %>%
  summarise_at(name.tbl[match(icbs,name.tbl$alt),"ref"], mean)

mean.exprs <- exprs %>%
  dplyr::select(c(name.tbl[match(f.order,name.tbl$alt),"ref"],"ROI","Phenograph")) %>% group_by(ROI,Phenograph) %>%
  summarise_at(name.tbl[match(f.order,name.tbl$alt),"ref"], mean)

cor.test <- data.frame(ct1=as.character(), ct2=as.character(),app=as.character(),icb=as.character(),
                       cor=as.numeric(),p=as.numeric())
for (ct1 in ct.order[c(1,2,4,5,6,8,9,10,11,12)]){
  for(ct2 in ct.order[c(1,2,4,5,6,8,9,10,11,12)]){
    for(app in apps){
      for(icb in icbs){
        dat1 <- mean.exprs[mean.exprs$Phenograph==ct1,c(name.tbl[match(app,name.tbl$alt),"ref"],"ROI")]
        dat2 <- mean.exprs[mean.exprs$Phenograph==ct2,c(name.tbl[match(icb,name.tbl$alt),"ref"],"ROI")]
        dat <- inner_join(dat1, dat2,by="ROI")
        test <- cor.test(dat[,1]%>%pull,dat[,3]%>%pull)
        cor.test[nrow(cor.test)+1,] <- c(ct1, ct2, app, icb, test$estimate, test$p.value)
      }
    }
  }
}
cor.test <- cor.test %>% mutate(cor=as.numeric(cor)) %>% mutate(p=as.numeric(p))

mean.exprs <- exprs %>% 
  dplyr::select(c(name.tbl[match(f.order,name.tbl$alt),"ref"],"ROI","Phenograph","X115In_IL6.ome")) %>% 
  group_by(ROI,Phenograph) %>%
  summarise_at(c(name.tbl[match(f.order,name.tbl$alt),"ref"],"X115In_IL6.ome"), mean)

cor.il6 <- data.frame(ct1=as.character(), ct2=as.character(),gene=as.character(),
                      cor=as.numeric(),p=as.numeric())
for (ct1 in ct.order[c(1,2,4,5,6,8,9,10,11,12)]){
  for (ct2 in ct.order[c(1,2,4,5,6,8,9,10,11,12)]){
    for(gene in c(apps,icbs)){
      dat1 <- mean.exprs[mean.exprs$Phenograph==ct1,c(name.tbl[match(gene,name.tbl$alt),"ref"],"ROI")]
      dat2 <- mean.exprs[mean.exprs$Phenograph==ct2,c("X115In_IL6.ome","ROI")]
      dat <- inner_join(dat1, dat2,by="ROI")
      test <- cor.test(dat[,1]%>%pull,dat[,3]%>%pull)
      cor.il6[nrow(cor.il6)+1,] <- c(ct1, ct2, gene, test$estimate, test$p.value)
    }
  }
}
cor.il6 <- cor.il6 %>% mutate(cor=as.numeric(cor)) %>% mutate(p=as.numeric(p))
cor.il6$label <- paste0(cor.il6$ct1,"(", cor.il6$gene, ")")
View(cor.il6 %>% filter(p<.05))


app.order <- paste0(rep(ct.order[c(1,2,4,5,6,8,9,10,11,12)], each=8), "(", rep(apps,10),")")
icb.order <- paste0(rep(ct.order[c(1,2,4,5,6,8,9,10,11,12)], each=7), "(", rep(icbs,10),")")
cor.test$label_y <- paste0(cor.test$ct1,"(", cor.test$app, ")")
cor.test$label_x <- paste0(cor.test$ct2,"(", cor.test$icb, ")")
cor.test <- cor.test%>% mutate(label_x=factor(label_x, levels = icb.order)) %>%
  mutate(label_y=factor(label_y, levels = app.order))

p <- ggplot(data=cor.test %>% filter(p<0.05) %>% filter(!is.na(label_x)) %>% filter(!is.na(label_y)), aes(x=label_x,y=label_y))+
  geom_tile(aes(fill=cor), color="grey80")+
  scale_fill_paletteer_c("grDevices::Fall", "Correlation", direction = 1, limits=c(-1,1))+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-5))+
  coord_equal()+theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))

grDevices::cairo_pdf("figures/APP&ICBRegulation.pdf", height=12, width=10)
print(p)
dev.off()

p1 <- ggplot(cor.il6 %>% filter(label %in% app.order) %>% 
  mutate(label=factor(label, levels = app.order)) %>%
    mutate(ct2=factor(ct2,levels=ct.order[c(1,2,4,5,6,8,9,10,11,12)])),aes(x=label,y=ct2))+
  geom_point(aes(color=cor,size=ifelse(p<0.05,"p<0.05","p≥0.05")))+
  scale_size_manual(values = c(5,1),name="p-value")+
  scale_color_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-1,1),na.value = NA,
                        name="Correlation\nCoefficients")+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-5))+
  coord_equal()+theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))

p2 <- ggplot(cor.il6 %>% filter(label %in% icb.order) %>% 
         mutate(label=factor(label, levels = icb.order)) %>%
                  mutate(ct2=factor(ct2,levels=ct.order[c(1,2,4,5,6,8,9,10,11,12)])),aes(x=label,y=ct2))+
  geom_point(aes(color=cor,size=ifelse(p<0.05,"p<0.05","p≥0.05")))+
  scale_size_manual(values = c(5,1),name="p-value")+
  scale_color_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-1,1),na.value = NA,
                        name="Correlation\nCoefficients")+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-5))+
  coord_equal()+theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))


grDevices::cairo_pdf("figures/IL6RegulationAPP.pdf", height=15, width=15)
print(p1)
dev.off()
grDevices::cairo_pdf("figures/IL6RegulationICB.pdf", height=15, width=15)
print(p2)
dev.off()

mean.exprs <-
  dplyr::select(c(name.tbl[match(f.order,name.tbl$alt),"ref"],"ROI","Phenograph","X171Yb_TNFa.ome")) %>% group_by(ROI,Phenograph) %>%
  summarise_at(c(name.tbl[match(f.order,name.tbl$alt),"ref"],"X171Yb_TNFa.ome"), mean)

cor.tnfa <- data.frame(ct1=as.character(), ct2=as.character(),gene=as.character(),
                      cor=as.numeric(),p=as.numeric())
for (ct1 in ct.order[c(1,2,4,5,6,8,9,10,11,12)]){
  for (ct2 in ct.order[c(1,2,4,5,6,8,9,10,11,12)]){
    for(gene in c(apps,icbs)){
      dat1 <- mean.exprs[mean.exprs$Phenograph==ct1,c(name.tbl[match(gene,name.tbl$alt),"ref"],"ROI")]
      dat2 <- mean.exprs[mean.exprs$Phenograph==ct2,c("X171Yb_TNFa.ome","ROI")]
      dat <- inner_join(dat1, dat2,by="ROI")
      test <- cor.test(dat[,1]%>%pull,dat[,3]%>%pull)
      cor.tnfa[nrow(cor.tnfa)+1,] <- c(ct1, ct2, gene, test$estimate, test$p.value)
    }
  }
}
cor.tnfa <- cor.tnfa %>% mutate(cor=as.numeric(cor)) %>% mutate(p=as.numeric(p))
cor.tnfa$label <- paste0(cor.tnfa$ct1,"(", cor.tnfa$gene, ")")

p1 <- ggplot(cor.tnfa %>% filter(label %in% app.order) %>% 
               mutate(label=factor(label, levels = app.order)) %>%
               mutate(ct2=factor(ct2,levels=ct.order[c(1,2,4,5,6,8,9,10,11,12)])),aes(x=label,y=ct2))+
  geom_point(aes(color=cor,size=ifelse(p<0.05,"p<0.05","p≥0.05")))+
  scale_size_manual(values = c(5,1),name="p-value")+
  scale_color_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-1,1),na.value = NA,
                        name="Correlation\nCoefficients")+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-5))+
  coord_equal()+theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))

p2 <- ggplot(cor.tnfa %>% filter(label %in% icb.order) %>% 
               mutate(label=factor(label, levels = icb.order)) %>%
               mutate(ct2=factor(ct2,levels=ct.order[c(1,2,4,5,6,8,9,10,11,12)])),aes(x=label,y=ct2))+
  geom_point(aes(color=cor,size=ifelse(p<0.05,"p<0.05","p≥0.05")))+
  scale_size_manual(values = c(5,1),name="p-value")+
  scale_color_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-1,1),na.value = NA,
                        name="Correlation\nCoefficients")+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-5))+
  coord_equal()+theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))


grDevices::cairo_pdf("figures/TnfRegulationAPP.pdf", height=15, width=15)
print(p1)
dev.off()
grDevices::cairo_pdf("figures/TnfRegulationICB.pdf", height=15, width=15)
print(p2)
dev.off()

mean.exprs <- exprs %>% filter(microsatellite_instability=="MSI-H") %>%
  dplyr::select(c(name.tbl[match(f.order,name.tbl$alt),"ref"],"ROI","Phenograph","X172Yb_stat3.ome")) %>% 
  group_by(ROI,Phenograph) %>%
  summarise_at(c(name.tbl[match(f.order,name.tbl$alt),"ref"],"X172Yb_stat3.ome"), mean)

cor.stat3 <- data.frame(ct1=as.character(), ct2=as.character(),gene=as.character(),
                       cor=as.numeric(),p=as.numeric())
for (ct1 in ct.order[c(1,2,4,5,6,8,9,10,11,12)]){
  for (ct2 in ct.order[c(1,2,4,5,6,8,9,10,11,12)]){
    for(gene in c(apps,icbs)){
      dat1 <- mean.exprs[mean.exprs$Phenograph==ct1,c(name.tbl[match(gene,name.tbl$alt),"ref"],"ROI")]
      dat2 <- mean.exprs[mean.exprs$Phenograph==ct2,c("X172Yb_stat3.ome","ROI")]
      dat <- inner_join(dat1, dat2,by="ROI")
      test <- cor.test(dat[,1]%>%pull,dat[,3]%>%pull)
      cor.stat3[nrow(cor.stat3)+1,] <- c(ct1, ct2, gene, test$estimate, test$p.value)
    }
  }
}
cor.stat3 <- cor.stat3 %>% mutate(cor=as.numeric(cor)) %>% mutate(p=as.numeric(p))
cor.stat3$label <- paste0(cor.stat3$ct1,"(", cor.stat3$gene, ")")

p1 <- ggplot(cor.stat3 %>% filter(label %in% app.order) %>% 
               mutate(label=factor(label, levels = app.order)) %>%
               mutate(ct2=factor(ct2,levels=ct.order[c(1,2,4,5,6,8,9,10,11,12)])),aes(x=label,y=ct2))+
  geom_point(aes(color=cor,size=ifelse(p<0.05,"p<0.05","p≥0.05")))+
  scale_size_manual(values = c(5,1),name="p-value")+
  scale_color_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-1,1),na.value = NA,
                        name="Correlation\nCoefficients")+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-5))+
  coord_equal()+theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))

p2 <- ggplot(cor.stat3 %>% filter(label %in% icb.order) %>% 
               mutate(label=factor(label, levels = icb.order)) %>%
               mutate(ct2=factor(ct2,levels=ct.order[c(1,2,4,5,6,8,9,10,11,12)])),aes(x=label,y=ct2))+
  geom_point(aes(color=cor,size=ifelse(p<0.05,"p<0.05","p≥0.05")))+
  scale_size_manual(values = c(5,1),name="p-value")+
  scale_color_gradientn(colors=rev(paletteer_d("RColorBrewer::RdYlBu")),limits=c(-1,1),na.value = NA,
                        name="Correlation\nCoefficients")+
  theme_bw()+xlab("")+ylab("")+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.margin=margin(t=-5))+
  coord_equal()+theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))


grDevices::cairo_pdf("figures/Stat3RegulationAPPMSI.pdf", height=15, width=15)
print(p1)
dev.off()
grDevices::cairo_pdf("figures/Stat3RegulationICBMSI.pdf", height=15, width=15)
print(p2)
dev.off()


stroke<-Read10X("/home/linux/GY/data_ann/gwas_aging/mydata/stroke_ad/GSE167593_stroke/")
stroke<-CreateSeuratObject(stroke,project = "stroke")
stroke$"condition"<-"HS"
stroke@meta.data[colnames(stroke)[grep("-1",colnames(stroke))],"condition"]<-"WT_sc"


stroke@meta.data[colnames(stroke)[grep("-2",colnames(stroke))],"condition"]<-"IS"

wt_sn<-Read10X("/home/linux/GY/data_ann/gwas_aging/mydata/stroke_ad/GSM4160645_WT_cor/")

wt_sn<-CreateSeuratObject(wt_sn,project = "WT_sn")

ad_cor<-Read10X("/home/linux/GY/data_ann/gwas_aging/mydata/stroke_ad/GSM4160645_WT_5XFAD/")
ad_cor<-CreateSeuratObject(ad_cor,project = "ad_cor")


adstroke<-merge(stroke,y=c(wt_sn,ad_cor))

head(adstroke@meta.data)
adstroke$condition<-recode(adstroke$condition,
                           "HS1"="HS")

selected_c<-WhichCells(adstroke,expression = nFeature_RNA>200)
selected_f<-rownames(adstroke)[Matrix::rowSums(adstroke@assays$RNA@counts>0)>3]
adstroke<-subset(adstroke, features=selected_f, cells=selected_c)
adstroke<-NormalizeData(adstroke,normalization.method = "LogNormalize",scale.factor = 10000)


adstroke<-FindVariableFeatures(adstroke,selection.method = "vst",nfeatures = rownames(adstroke))


all.genes<-rownames(adstroke)
adstroke<-ScaleData(adstroke,features = all.genes)


adstroke<-SCTransform(adstroke,vars.to.regress = c("nFeature_RNA"))%>%
  RunPCA(npcs = 50)

adstroke<-RunPCA(adstroke,features = VariableFeatures(object = adstroke))



adstroke=adstroke %>% RunHarmony("study", plot_convergence = TRUE,
                                 reduction="pca",project.dim=F)

adstroke<-adstroke %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) 


adstroke<-FindNeighbors(adstroke, reduction = "pca",dims = 1:30)%>%
  FindClusters(resolution = 0.4)  
## Sonia Leach
## 11/25/2026

## Code to reproduce results for Kleist et al, 2027 manuscript 
## (see object details and sessionInfo() at end)

## conda activate R443

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratExtend)    ## for theme_umap_arrows()
library(SeuratWrappers)  ## for imputation by RunALRA()
library(RColorBrewer)    ## for colorRampPalette()

set.seed(1234)
DATADIR <- 'sample_filtered_feature_bc_matrix'

DATADIR = '/dartfs-hpc/rc/home/w/f007t2w/Projects/Rosato/GEO_submission/250626_10x__mouseOT1'
##============================
## RNA expression
##============================
data <- Read10X(data.dir = DATADIR)
data[["Antibody Capture"]]@Dimnames[[1]][1] <- 'Not.Used'
data[["Antibody Capture"]]@Dimnames[[1]][2] <- 'TrmOnlyTumor'
data[["Antibody Capture"]]@Dimnames[[1]][3] <- 'TcircmOnlyTumor'
data[["Antibody Capture"]]@Dimnames[[1]][4] <- 'IntactSpleen'

obj <- CreateSeuratObject(counts = data[["Gene Expression"]],
                          min.cells = 3, min.features = 200,
                          project = 'MouseOT1')
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^mt-")
obj <- subset(obj,
              subset = (nFeature_RNA > 500   &
                        nFeature_RNA < 8520  &  ## MAD=5
                        nCount_RNA   < 73148 &  ## MAD=3
                        percent.mt   < 4        ## MAD=3
              ))
dim(obj) ## [1] 23087  7390

##============================
## HTO hashes - Samples
##============================
obj[['sample']] <- CreateAssayObject(
    counts = data$`Antibody Capture`[1:4,colnames(obj)])
obj <- NormalizeData(obj, assay = "sample", normalization.method = "CLR")
obj <- HTODemux(obj, assay = "sample", positive.quantile = 0.99)
## Cutoff for Not.Used : 9 reads
## Cutoff for TrmOnlyTumor : 1064 reads <-- note high threshold
## Cutoff for TcircmOnlyTumor : 20 reads
## Cutoff for IntactSpleen : 9 reads
Idents(obj) <- obj$hash.ID

## Establish where TrmOnlyTumor highest in Doublet/Negative calls
vals <- paste0(obj$hash.ID, ':', obj$sample_maxID)
where <- which(obj$hash.ID == obj$sample_maxID)
## Keep all proper singlet calls
vals[where] <- as.character(obj$hash.ID[where])  
## Simplify names for Negative/Doublet calls
vals <- gsub('Negative', 'Neg',
             gsub('Doublet', 'Dbl',
                  gsub('OnlyTumor', '',
                       gsub('Intact', '',vals)))) 
vals <- factor(vals)
obj$hashMax <- vals

## Manually recover more TrmOnly candidates based on normalized read
## count for doublet/negatives [thresholds determined by extensive
## investigation]
trmOnlydata <- obj[["sample"]]@data['TrmOnlyTumor',]
myAssgn <- obj$hash.ID
where <- ((obj$hashMax=='Dbl:Trm' & trmOnlydata>4) |
         ( obj$hashMax=='Neg:Trm' & trmOnlydata>1))
myAssgn[where] <- 'TrmOnlyTumor'
obj$hash.IDmanual <- myAssgn

## Keep only singlets and change sample name in orig.ident
obj <- subset(obj,
              subset = hash.IDmanual %in% c('Doublet', 'Negative'),
              invert=TRUE)
obj$orig.ident <- obj$hash.IDmanual
table(obj$orig.ident)
##TcircmOnlyTumor     Not.Used    TrmOnlyTumor    IntactSpleen
##           1461         1370            1634            1741
Idents(obj) <- 'orig.ident'
dim(obj) ##  [1] 23087  6206

## Remove 'Not.Used' cells
obj <- subset(obj, subset=orig.ident=='Not.Used', invert=TRUE)
obj$orig.ident <- factor(obj$orig.ident,
                        levels=c('TrmOnlyTumor',
                                 'TcircmOnlyTumor',
                                 'IntactSpleen'))

table(obj$orig.ident)
##   TrmOnlyTumor TcircmOnlyTumor    IntactSpleen
##           1634            1461            1741


##============================
## HTO hashes - Antibodies
##============================
obj[['HTO']] <- CreateAssayObject(
    counts = data$`Antibody Capture`[5:6,colnames(obj)])
obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
## Cutoff for CD103 : 30 reads
## Cutoff for CD69 : 46 reads

##============================
## Do standard SCT analysis
##============================
options(future.globals.maxSize = 1000 * 1024^2)  ## xMB x 1024^2 [x=500 default]
set.seed(143)
obj <- SCTransform(obj, assay="RNA", method = "glmGamPoi",
                   vars.to.regress=c("percent.mt"),
                   conserve.memory = TRUE,
                   return.only.var.genes = FALSE)
markers.remove <- c(grep(pattern = "^Trav|^Trbv|^Trgv|^Trdv|^Trbc|^Trac|^Trdc",
                         rownames(obj), value = TRUE), ## TCR genes
                    c("Gm42418",'AY036118','Lars2')) ## lnRNA.Junk
VariableFeatures(obj) <- setdiff(VariableFeatures(obj), markers.remove)
obj <- PrepSCTFindMarkers(obj) ## Do before DE analysis
obj <- RunPCA(obj, npcs = 30)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- RunUMAP(obj, reduction = 'pca', dims = 1:30)
obj <- FindClusters(obj, resolution=0.6, verbose=FALSE)

## DimPlot should look like a screaming baby chicken with tongue (cluster 9)
## and a little cluster in front of chick (cluster 6)
DimPlot(obj, label=TRUE, label.size=8) +
    geom_vline(xintercept=5) +
    theme(aspect.ratio=1)

## Remove junk clusters and random cells off main cell mass 
obj <- subset(obj,
             subset = (umap_1 > 5) | (SCT_snn_res.0.6 %in% c('6', '9')),
             invert=TRUE)
table(obj$orig.ident)
##   TrmOnlyTumor TcircmOnlyTumor    IntactSpleen
##           1540            1400            1737

## Visualize again
DimPlot(obj, label=TRUE, label.size=8) +
    geom_vline(xintercept=5) +
    theme(aspect.ratio=1)

##============================
## Impute with ALRA
##============================

start_time <- Sys.time()
## ALRA uses DefaultAssay(obj)$data slot
obj <- RunALRA(obj, use.mkl = FALSE, mkl.seed = 143)  
print(difftime(Sys.time(), start_time, units = "mins"))
obj <- RenameAssays(obj, alra = 'ALRA')
## Rank k = 45

## 0.03% of the values became negative in the scaling process and were
## set to zero
## The matrix went from 20.45% nonzero to 50.67% nonzero
## Time difference of 0.4858825 mins
## Renaming default assay from alra to ALRA

##============================
## Score Residency and Circulating
##============================

## Lists taken from PMID:35537413
## Jaiswal_SuppInfo2_ViralTcellSigs_1-s2.0-S1535610822001672-mmc2.xlsx
## columns N/O Non-exclusive TRMup TRMdn [transferred from human to
## mouse via homologene v1.4.68.19.3.27 R library]

## Jaiswal Non-exclusive TRM-up, mouse homologs
Residency <- c("AA467197", "Abcb1a", "Abr", "Adam19", "Aim1", "Aldoa",
              "Aplp2", "Arf6", "Arrdc3", "Arsb", "Atf4", "Atf6",
              "Atp2b4", "B4galnt4", "Bag3", "Baiap2", "Bcl2a1d",
              "Bhlhe40", "Btg2", "Ccdc107", "Ccl5", "Ccr2", "Cd244",
              "Cd28", "Cdh1", "Cdkn1a", "Chn2", "Cish", "Crbn",
              "Crem", "Csrnp1", "Cst7", "Ctla4", "Ctsw", "Cxcr3",
              "Cxcr6", "D16Ertd472e", "Dapk2", "Ddit3", "Ddx3x",
              "Dennd4a", "Dennd5a", "Dnaja1", "Dnajb4", "Dnajb6",
              "Dock5", "Dusp1", "Ell2", "Erbb2", "Fam204a", "Fasl",
              "Fgl2", "Fnbp1", "Fosl2", "Fth1", "Gadd45b", "Gla",
              "Glmn", "Gpr171", "Gpr34", "Gpr65", "Gpr68", "Gzma",
              "Gzmk", "H2-M3", "H2-T23", "Hist1h1c", "Hopx", "Hpgds",
              "Hsp90ab1", "Hspa5", "Hspa8", "Hsph1", "Hyou1", "Icos",
              "Id2", "Ifng", "Igf2r", "Il18rap", "Il2ra", "Il2rb",
              "Il4ra", "Inpp4b", "Inpp5d", "Irf4", "Isg20", "Itga1",
              "Itgae", "Itgax", "Jun", "Junb", "Klf6", "Litaf",
              "Mapk6", "Mapkapk3", "Mtmr1", "Nabp1", "Ndfip1",
              "Neurl3", "Nfil3", "Nkg7", "Notch4", "Nr4a1", "Nr4a2",
              "Nr4a3", "Nt5e", "Odc1", "Osbpl3", "P2ry10", "Pdia3",
              "Pelo", "Per1", "Pik3ap1", "Plek", "Plekha5", "Plk3",
              "Ppp2r2d", "Ppp3cc", "Ptger4", "Ptges3", "Ptk2b",
              "Ptp4a1", "Ptpn22", "Ptprj", "Pygl", "Qpct", "Rab8b",
              "Rgs1", "Rgs2", "Rnf149", "Rpl30", "Runx2", "S100a6",
              "Samd3", "Sat1", "Serpinb9", "Sh2d2a", "Skil", "Slamf7",
              "Smim3", "Spty2d1", "Srsf3", "St3gal6", "Suco", "Sympk",
              "Sytl2", "Sytl3", "Tbx21", "Tesc", "Tigit", "Tinf2",
              "Tiparp", "Tm2d2", "Tmed5", "Tmem2", "Tnfaip3",
              "Tnfrsf1b", "Tob2", "Ttc39b", "Ubc", "Vps37b", "Ybx3",
              "Zeb2", "Zfp36l2")
## Jaiswal Non-exclusive TRM-dn, mouse homologs
Circulating <- c("Actn1", "Adcy6", "Add1", "Add3", "Aff3", "Ampd1",
                "Arhgap15", "Arhgef18", "Arid4a", "Art4", "Atic",
                "Atp1b1", "Atp1b3", "Atp8a1", "Bach2", "Baz2b",
                "Ccr7", "Cd2ap", "Daf2", "Chia1", "Cnga1", "Dapl1",
                "Dgka", "Dnah8", "Dnm2", "Dock11", "Dsg1b", "Emc9",
                "Eng", "Epha6", "Ephx1", "Fam53b", "Fam63b", "Fam65b",
                "Galnt10", "Gin1", "Slc28a2", "Gramd4", "Gria3",
                "H2-Oa", "H2-Ob", "Hectd2", "Icam2", "Ift80", "Ikbke",
                "Il6ra", "Il6st", "Kbtbd11", "Klhdc1", "Lef1",
                "Mab21l2", "Map4k2", "Mboat1", "Mogat1", "Ms4a4a",
                "N4bp2", "Numb", "Oma1", "Parp9", "Pdk1", "Pik3ip1",
                "Pole2", "Polr3b", "Pot1a", "Prkcq", "Rapgef4",
                "Rasa3", "Rasgrp2", "Rassf2", "Rbks", "Rnf144a",
                "Rpf2", "S1pr1", "Scml4", "Sele", "Sell", "Setx",
                "Sh3bp5", "Sh3pxd2a", "Sidt1", "Slamf6", "Slfn5",
                "Smad5", "Smc4", "Sptbn1", "Srpk1", "Stk38", "Tcf4",
                "Tet1", "Tlr1", "Traf3ip3", "Treml2", "Tspan32",
                "Ttc13", "Usp28", "Usp33", "Wdr37", "Xkrx", "Ylpm1",
                "Zdhhc13", "Zfc3h1", "Zfp82")

GeneLists <- list(Residency=Residency, Circulating=Circulating)
## Subset lists to genes measured
GeneLists <- lapply(GeneLists, intersect, rownames(obj)) 

DefaultAssay(obj) <- 'ALRA'  ## use imputed data for scoring
for (genelist in names(GeneLists)) {
    cat('Scoring list', genelist, '\n')
    obj <- AddModuleScore(object = obj,
                          features = GeneLists[genelist],
                          name = paste0('ModuleScore_', genelist))
}
colnames(obj@meta.data) <- gsub('^ModuleScore\\_(.+)1$', 'Score.\\1',
                               colnames(obj@meta.data))

##=========================
## Establish color scheme
##=========================
sampleColors <- setNames(c(gold =rgb(228,169,70, maxColorValue=255),
                           slate=rgb(97,141,183, maxColorValue=255),
                           teal =rgb(139,201,184, maxColorValue=255)),
                         c("TcircmOnlyTumor","TrmOnlyTumor","IntactSpleen"))
fig.samplenames <- setNames(c(expression(T[CircM]*"-only Tumor"),
                         expression(T[RM]*"-only Tumor"),
                         expression('Spleen')),
                       c("TcircmOnlyTumor","TrmOnlyTumor","IntactSpleen"))

col20 <- c(darkblue="#1F77B4FF",
           lightblue="#AEC7E8FF",
           darkorange="#FF7F0EFF",
           lightorange="#FFBB78FF",
           darkgreen="#2CA02CFF",
           lightgreen="#98DF8AFF",
           darkred="#D62728FF",
           lightred="#FF9896FF",
           darkpurple="#9467BDFF",
           lightpurple = "#C5B0D5FF",
           darkbrown="#8C564BFF",
           lightbrown="#C49C94FF",
           darkpink = "#E377C2FF",
           lightpink = "#F7B6D2FF",
           darkgrey = "#7F7F7FFF",
           lightgrey = "#C7C7C7FF",
           darkyellow = "#BCBD22FF",
           lightyellow = "#DBDB8DFF",
           darkaqua = "#17BECFFF",
           lightaqua = "#9EDAE5FF")

## Reorder cluster names and numbers based on identities
obj$ClusterNums <- paste0('C', as.character(obj$seurat_clusters))
Idents(obj) <- 'ClusterNums'
table(Idents(obj))
##  C2   C0   C3   C4   C1   C5   C8   C7
## 775 1153  693  650 1137  142   57   70
obj <- RenameIdents(obj,
                   'C3' = 0,
                   'C2' = 1,
                   'C0' = 2,
                   'C5' = 3,
                   'C8' = 4,
                   'C7' = 5,
                   'C1' = 6,
                   'C4' = 7)
obj$ClusterNums <- Idents(obj)
figlegend.samplenames <-c(
    '0' = expression("0 TIL CD103+"~T[RM]),
    '1' = expression("1 TIL ISG"~T[RM]),
    '2' = expression("2 TIL"~T[CM]),
    '3' = expression("3 TIL"~Tox^{hi}),
    '4' = expression("4 TIL Proliferating"),
    '5' = expression("5 Transitioning"),
    '6' = expression("6 Spleen"~T[EM]),
    '7' = expression("7 Spleen"~T[CM])
)
fig.colors <- c(
    "0" = unname(col20["darkred"]),    ## TIL CD103+ Trm
    "1" = unname(col20["darkorange"]), ## TIL ISG Trm
    "2" = unname(col20["darkblue"]),   ## TIL Tcm
    "3" = unname(col20["lightaqua"]),  ## TIL Tox.hi
    "4" = unname(col20["lightgreen"]), ## TIL Proliferating
    "5" = unname(col20["darkgreen"]),  ## Transitioning
    "6" = unname(col20["darkpurple"]), ## Spleen Tem
    "7" = unname(col20["darkpink"])    ## Spleen Tcm
)

## Visualize color settings, etc
DimPlot(obj, group.by='ClusterNums', split.by='orig.ident',
        pt.size=0.1, cols=fig.colors, label=TRUE,label.size=8) +
    scale_color_manual(values = fig.colors,
                       labels = figlegend.samplenames) +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size = 5)))  +
    theme_umap_arrows() + ggtitle('')


##===================================================================
##    Visualizations
##===================================================================

##===================================================================
## (.) UMAP of transcriptional states annotated by cluster name.

## (.) UMAP of OT-I T cells annotated by experimental group/tissue.
##===================================================================
pt.size=0.6
obj$ClusterNums <- factor(obj$ClusterNums, levels=as.character(0:7))
p.UMAP <- DimPlot(obj, group.by='ClusterNums',
                 pt.size=pt.size, cols=fig.colors,
                 label=TRUE,label.size=8) +
    scale_color_manual(values = fig.colors,
                       labels = figlegend.samplenames) +  
    theme(aspect.ratio=1,
          plot.margin = margin(10, 10, 20, 20, "pt"),
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 14)   
          ) +
    guides(
        color = guide_legend(override.aes = list(size = 5)) 
    )  +
    theme_umap_arrows(anchor_x = unit(6, "mm"),
                      anchor_y = unit(6, "mm")
                      ) +
    ggtitle('')

p.UMAP.samp <- DimPlot(obj, group.by='orig.ident',
                   pt.size=pt.size, cols=sampleColors) +
    scale_color_manual(values = sampleColors,
                       labels = fig.samplenames) +  
    theme(aspect.ratio=1,
          plot.margin = margin(10, 10, 20, 20, "pt"),
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 14)   
          ) +
    guides(
        color = guide_legend(override.aes = list(size = 5)) # Increase Shape Size
    )  +
    theme_umap_arrows() +
    ggtitle('')


ggsave(filename="PLOT_UMAPclusters.png", plot=p.UMAP,
       width=8, height=8, units="in")
ggsave(filename="PLOT_UMAPclusters.eps", plot=p.UMAP,
       width=8, height=8, units="in", dpi=300)

ggsave(filename="PLOT_UMAPsamples.png", plot=p.UMAP.samp,
       width=8, height=8, units="in")
ggsave(filename="PLOT_UMAPsamples.eps", plot=p.UMAP.samp,
       width=8, height=8, units="in", dpi=300)

##===================================================================
## (.) Relative frequencies of OT-I T cell clusters in each
## experimental group/tissue.

## (.) Relative frequencies of each experimental group/tissue in each
## cluster.
##===================================================================
obj$ClusterNums <- factor(obj$ClusterNums, levels=as.character(0:7))
Idents(obj) <- 'ClusterNums'

## Cluster order changed to emphasize differences between samples
tempobj <- data.frame(Treatment=obj@meta.data[,'orig.ident'],
                     CellType=obj@meta.data[,'ClusterNums'])
tempobj$CellType <- factor(tempobj$CellType,
                          levels = c("2","0","1","3","4","5","7","6"))
p.SampleCompBars <-
    ggplot(tempobj, aes(x=Treatment, fill=CellType)) +
    geom_bar(position='fill', color='black') +  ## scaled to 100%
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_flip() + theme_minimal() +
    ylab('% of sample composition') + xlab('') +
    theme(axis.text=element_text(size=20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold"),
          legend.title=element_text(size=18),
          legend.text=element_text(size=16)
          ) +
    scale_fill_manual(values=fig.colors, name='Cluster') +
    scale_x_discrete(labels = fig.samplenames) +
    ggtitle('Sample Composition')

## Cluster order changed to emphasize differences between samples
tempobj$CellType <- factor(tempobj$CellType,
                          levels = c("0","2","4","1","3","5","7","6"))
p.ClusterCompBars <-
    ggplot(tempobj, aes(x=CellType, fill=Treatment)) +
    geom_bar(position = position_fill(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_flip() + theme_minimal() +
    ylab('% of cluster composition ') + xlab('') +
    theme(axis.text=element_text(size=20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold"),
          legend.title=element_text(size=18),
          legend.text=element_text(size=16)
          ) +
    scale_fill_manual(values=sampleColors, name='',
                      labels=fig.samplenames) +
    ggtitle('Cluster Composition')


ggsave('PLOT_SampleCompositionBarPlots.png', p.SampleCompBars,
       height=4, width=9, units='in')
ggsave('PLOT_SampleCompositionBarPlots.eps', p.SampleCompBars,
       height=4, width=9, units='in', dpi=300)

ggsave('PLOT_ClusterCompositionBarPlots.png', p.ClusterCompBars,
       height=4, width=9, units='in')
ggsave('PLOT_ClusterCompositionBarPlots.eps', p.ClusterCompBars,
       height=4, width=9, units='in', dpi=300)


##===================================================================
## (.) Dot plot of select circulating and resident memory T cell
## transcripts across clusters.
##===================================================================

features.ResCirc = c(
    ##~~~Residency
    'Zfp683', ## Zfp683 == 'Hobit'
    'Itgae',
    'Zap70',
    'Cxcr6',
    'Prdm1', ## Prdm1=='BLIMP1'
    'Runx3',
    'Dusp5',
    'Bhlhe40',
    'Cd69',
    ##~~~Circulating
    'Zeb2',
    'Klrg1',
    'Cx3cr1',
    'S1pr5',
    'S1pr4',
    'S1pr1',
    'Klf2',
    'Il7r',
    'Sell',
    'Ccr7'
)
lineIdx <- grep('Cd69', features.ResCirc) + 0.5

## Cluster order changed to emphasize differences between Res/Circ
dotplotOrder <- as.character(c( 7,2,5,6,  1,3,4,0))  ## left to right
obj$ClusterNums <- factor(obj$ClusterNums, levels=dotplotOrder)

dp.ResCirc <-
    DotPlot(obj, assay = 'ALRA',
                      group.by = 'ClusterNums',
                      features = features.ResCirc,
                      cluster.idents = FALSE, cols = 'RdBu',
                      col.max = 2, col.min = -2) +
    theme(axis.text.x = element_text(vjust = 0.5,
                                     hjust = 0.7,
                                     size = 14))+  ## cluster ids
    theme(axis.text.y = element_text(face = "italic",
                                     size = 14)) +  ## gene names
    coord_flip() +
    geom_vline(xintercept=lineIdx, col=rgb(36,122,190,maxColorValue=255),
               linetype='dotted', linewidth=0.7) +
    ggtitle('Memory T cell subsets')


ggsave(filename="PLOT_DotplotResCirc.png", plot=dp.ResCirc,
       width=6, height=5.2, units="in")
ggsave(filename="PLOT_DotplotResCirc.eps", plot=dp.ResCirc,
       width=6, height=5.2, units="in", dpi=300)

##===================================================================
## (.) Violin plots (left) and feature plots (right) 40 depicting
## circulating and residency module scores across clusters. Module
## scores generated from Jaiswal PMID:35537413
##===================================================================

## Violin Plots
scoreMedians <- obj@meta.data %>% select(contains('Score.')) %>%
    summarise(across(everything(), \(x) median(x, na.rm=TRUE))) %>%
    unlist()

vp.Circ <- VlnPlot(obj, features=c('Score.Circulating'),
                  group.by='ClusterNums',
             sort=TRUE, cols=fig.colors) +
        ggtitle('Circulating') +
    geom_hline(yintercept=scoreMedians['Score.Circulating'],
               linetype='dashed') +
    NoLegend() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 1),
          axis.text = element_text(size = rel(1.5)),
          plot.title = element_text(size = rel(1.5)))

vp.Res <- VlnPlot(obj, features=c('Score.Residency'), group.by='ClusterNums',
                 sort=TRUE, cols=fig.colors) +
    ggtitle('Residency') +
    geom_hline(yintercept=scoreMedians['Score.Residency'],
               linetype='dashed') +
    NoLegend() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 1),
          axis.text = element_text(size = rel(1.5)),
          plot.title = element_text(size = rel(1.5)))


ggsave(filename="PLOT_VlnPlotCirculating.png", vp.Circ,
       width=8, height=3, units='in')
ggsave(filename="PLOT_VlnPlotCirculating.eps", vp.Circ,
       width=8, height=3, units='in', dpi=300)

ggsave(filename="PLOT_VlnPlotResidency.png", vp.Res,
       width=8, height=3, units='in')
ggsave(filename="PLOT_VlnPlotResidency.eps", vp.Res,
       width=8, height=3, units='in', dpi=300)


## Feature Plots
rwb <- colorRampPalette(colors = c("red", "white", "blue"))
colors.RWB   <- rev(rwb(7))

fp.Circ <- FeaturePlot(obj,features='Score.Circulating',
                      order=TRUE, pt.size=0.5) +
    scale_colour_gradientn(colours = colors.RWB, name='Score') +
    ggtitle('Circulating') + 
    theme(aspect.ratio=1,
          plot.title = element_text(size = 20)) +
    theme_umap_arrows( anchor_x = unit(3.2, "mm"),
                      anchor_y = unit(3.2, "mm"),
                      line_length = unit(14, "mm"),
                      text_offset_x = unit(1, "mm"),
                      text_offset_y = unit(1, "mm"),
                      text_size = 7)


fp.Res <- FeaturePlot(obj,features='Score.Residency',
                     order=TRUE, pt.size=0.5) +
    scale_colour_gradientn(colours = colors.RWB, name='Score') +
    ggtitle('Residency') + 
    theme(aspect.ratio=1,
          plot.title = element_text(size = 20)) + 
    theme_umap_arrows( anchor_x = unit(3.2, "mm"),
                      anchor_y = unit(3.2, "mm"),
                      line_length = unit(14, "mm"),
                      text_offset_x = unit(1, "mm"),
                      text_offset_y = unit(1, "mm"),
                      text_size = 7)

ggsave(filename="PLOT_UMAPCirculating.png", fp.Circ,
       width=4, height=4, units='in')
ggsave(filename="PLOT_UMAPCirculating.eps", fp.Circ,
       width=4, height=4, units='in', dpi=300)

ggsave(filename="PLOT_UMAPResidency.png", fp.Res,
       width=4, height=4, units='in')
ggsave(filename="PLOT_UMAPResidency.eps", fp.Res,
       width=4, height=4, units='in', dpi=300)

##===================================================================
## () Dot plot of select cluster-defining genes.
##===================================================================

## Cluster defining genes dotplot
features.ClusterDefn = list(
    Type1INF =
        c('Ifitm1','Ifitm2','Ifitm3',
          'Irf7','Isg15','Ifit3','Ifit1'),
    CytokineChemokine =
        c('Cxcl10','Ccl3','Ccl4','Ccl5',
          'Ifng',
          'Il2ra',  ## 
          'Tgfb1','Tnf'),
    Cytotoxic =
        c('Gzmk','Gzmb','Gzma',
          'Prf1', 'Klrg1', 'Nkg7', 'Fasl',
          'Klrk1', 'Klrd1'),
    InhibExhaust =
        c('Lag3','Ctla4','Pdcd1','Tigit',
          'Entpd1','Havcr2','Cd244a',
          'Slamf6'),
    TxFactors =
        c('Stat1','Stat3', 'Fos', 'Tbx21',
          'Foxo1','Tcf7','Tox', 'Hif1a' 
          ),
    CellCycle=
        c("Mki67", "Pcna", "Ccna2", "Top2a"
          )
)


## Cluster order changed to emphasize differences between Res/Circ
dotplotOrder <- as.character(c( 7,2,5,6,  1,3,4,0))  ## bottom to top
obj$ClusterNums <- factor(obj$ClusterNums, levels=dotplotOrder)

dp.ClustDefn <-
    DotPlot(obj, assay = 'ALRA',
                    group.by= 'ClusterNums',
                    features = features.ClusterDefn,
                    cols = 'RdBu',col.max = 1.5,
                    col.min=-1.5, cluster.idents=FALSE) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)) +
    theme(axis.text.x = element_text(size = 12,face = "italic"))


ggsave(filename="PLOT_DotplotClusterDefn.png",
       dp.ClustDefn, width=16, height=4, units="in")
ggsave(filename="PLOT_DotplotClusterDefn.eps",
       dp.ClustDefn, width=16, height=4, units="in", dpi=300)


##===================================================================
## (.) Volcano plot of DEGs between Cluster 0 OT-I T cells from
## TCircM-only v.s. TRM-only tumors.
##===================================================================

vo.0 <-
    VolcanoPlot(subset(obj, subset= ClusterNums=='0'),
                group.by='orig.ident',
                assay='SCT',
                ident.1 = "TrmOnlyTumor",
                ident.2 = "TcircmOnlyTumor",
                log.base = "2"    # Use log2 instead of natural log
                ) +
    labs(subtitle =figlegend.samplenames['0']) + 
    theme(
        plot.subtitle = element_text(hjust = 0.5)
    )


ggsave('PLOT_VolcanoPlotsTrmvsTcirc_CN0.png', vo.0,
       width=6, height=5, units='in')
ggsave('PLOT_VolcanoPlotsTrmvsTcirc_CN0.eps', vo.0,
       width=6, height=5, units='in', dpi=300)


##===================================================================
## Save processed data for GEO supplement
##===================================================================
obj
## An object of class Seurat 
## 58497 features across 4677 samples within 5 assays 
## Active assay: ALRA (17702 features, 0 variable features)
##  1 layer present: data
##  4 other assays present: RNA, sample, HTO, SCT
##  2 dimensional reductions calculated: pca, umap

colnames(obj@meta.data)
##  [1] "orig.ident"                   "nCount_RNA"                  
##  [3] "nFeature_RNA"                 "percent.mt"                  
##  [5] "nCount_sample"                "nFeature_sample"             
##  [7] "sample_maxID"                 "sample_secondID"             
##  [9] "sample_margin"                "sample_classification"       
## [11] "sample_classification.global" "hash.ID"                     
## [13] "hashMax"                      "hash.IDmanual"               
## [15] "nCount_HTO"                   "nFeature_HTO"                
## [17] "HTO_maxID"                    "HTO_secondID"                
## [19] "HTO_margin"                   "HTO_classification"          
## [21] "HTO_classification.global"    "nCount_SCT"                  
## [23] "nFeature_SCT"                 "SCT_snn_res.0.6"             
## [25] "seurat_clusters"              "Score.Residency"             
## [27] "Score.Circulating"            "ClusterNums"    

saveRDS(obj, 'seurat_obj.mouseOT1.Kleist2026.RDS')
write.csv(obj@meta.data, 'seurat_obj_metadata.mouseOT1.Kleist2026.csv')

##===================================================================
sessionInfo()
## R version 4.4.3 (2025-02-28)
## Platform: x86_64-conda-linux-gnu
## Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

## Matrix products: default
## BLAS/LAPACK: /dartfs-hpc/rc/home/w/f007t2w/.conda/envs/R443/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

## time zone: America/New_York
## tzcode source: system (glibc)

## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     

## other attached packages:
##  [1] ggrepel_0.9.6          tidyr_1.3.1            rlist_0.4.6.2         
##  [4] SeuratWrappers_0.4.0   future_1.67.0          RColorBrewer_1.1-3    
##  [7] cowplot_1.2.0          SeuratExtend_1.2.4     SeuratExtendData_0.3.0
## [10] dplyr_1.1.4            patchwork_1.3.2        ggplot2_4.0.0         
## [13] Seurat_5.3.1.9999      SeuratObject_5.2.0     sp_2.2-0              

## loaded via a namespace (and not attached):
##   [1] jsonlite_2.0.0              magrittr_2.0.4             
##   [3] ggbeeswarm_0.7.2            spatstat.utils_3.2-0       
##   [5] farver_2.1.2                zlibbioc_1.52.0            
##   [7] vctrs_0.6.5                 ROCR_1.0-11                
##   [9] DelayedMatrixStats_1.28.1   spatstat.explore_3.5-3     
##  [11] htmltools_0.5.8.1           S4Arrays_1.6.0             
##  [13] SparseArray_1.6.2           sctransform_0.4.2          
##  [15] parallelly_1.45.1           KernSmooth_2.23-26         
##  [17] htmlwidgets_1.6.4           ica_1.0-3                  
##  [19] plyr_1.8.9                  plotly_4.11.0              
##  [21] zoo_1.8-14                  igraph_2.2.1               
##  [23] mime_0.13                   lifecycle_1.0.4            
##  [25] pkgconfig_2.0.3             rsvd_1.0.5                 
##  [27] Matrix_1.7-4                R6_2.6.1                   
##  [29] fastmap_1.2.0               GenomeInfoDbData_1.2.13    
##  [31] MatrixGenerics_1.18.1       fitdistrplus_1.2-4         
##  [33] shiny_1.11.1                digest_0.6.37              
##  [35] S4Vectors_0.44.0            tensor_1.5.1               
##  [37] RSpectra_0.16-2             irlba_2.3.5.1              
##  [39] GenomicRanges_1.58.0        labeling_0.4.3             
##  [41] progressr_0.18.0            spatstat.sparse_3.1-0      
##  [43] httr_1.4.7                  polyclip_1.10-7            
##  [45] abind_1.4-8                 compiler_4.4.3             
##  [47] remotes_2.5.0               withr_3.0.2                
##  [49] S7_0.2.0                    fastDummies_1.7.5          
##  [51] R.utils_2.13.0              MASS_7.3-65                
##  [53] DelayedArray_0.32.0         tools_4.4.3                
##  [55] vipor_0.4.7                 lmtest_0.9-40              
##  [57] otel_0.2.0                  beeswarm_0.4.0             
##  [59] httpuv_1.6.16               future.apply_1.20.0        
##  [61] goftest_1.2-3               R.oo_1.27.1                
##  [63] glmGamPoi_1.18.0            glue_1.8.0                 
##  [65] nlme_3.1-168                promises_1.5.0             
##  [67] Rtsne_0.17                  cluster_2.1.8.1            
##  [69] reshape2_1.4.4              generics_0.1.4             
##  [71] gtable_0.3.6                spatstat.data_3.1-9        
##  [73] R.methodsS3_1.8.2           data.table_1.17.8          
##  [75] XVector_0.46.0              BiocGenerics_0.52.0        
##  [77] spatstat.geom_3.6-0         RcppAnnoy_0.0.22           
##  [79] RANN_2.6.2                  pillar_1.11.1              
##  [81] stringr_1.6.0               spam_2.11-1                
##  [83] RcppHNSW_0.6.0              later_1.4.4                
##  [85] splines_4.4.3               lattice_0.22-7             
##  [87] survival_3.8-3              deldir_2.0-4               
##  [89] tidyselect_1.2.1            miniUI_0.1.2               
##  [91] pbapply_1.7-4               gridExtra_2.3              
##  [93] IRanges_2.40.1              SummarizedExperiment_1.36.0
##  [95] scattermore_1.2             stats4_4.4.3               
##  [97] Biobase_2.66.0              matrixStats_1.5.0          
##  [99] stringi_1.8.7               UCSC.utils_1.2.0           
## [101] lazyeval_0.2.2              codetools_0.2-20           
## [103] tibble_3.3.0                BiocManager_1.30.26        
## [105] cli_3.6.5                   uwot_0.2.4                 
## [107] xtable_1.8-4                reticulate_1.44.0          
## [109] Rcpp_1.1.0                  GenomeInfoDb_1.42.3        
## [111] globals_0.18.0              spatstat.random_3.4-2      
## [113] png_0.1-8                   ggrastr_1.0.2              
## [115] spatstat.univar_3.1-4       parallel_4.4.3             
## [117] dotCall64_1.2               sparseMatrixStats_1.18.0   
## [119] listenv_0.10.0              viridisLite_0.4.2          
## [121] scales_1.4.0                ggridges_0.5.7             
## [123] crayon_1.5.3                purrr_1.2.0                
## [125] rlang_1.1.6                
 

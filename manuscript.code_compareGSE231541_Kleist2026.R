## Sonia Leach
## 11/25/2026

## Code to reproduce results for Kleist et al, 2027 manuscript 
## (see object details and sessionInfo() at end)

## conda activate R443


## Compares manuscript data to public data from GEO GSE231541:
## CITE-seq analysis of virus-specific memory P14 (α-viral) and
## tumor-specific OT-I (α-tumor) CD8+ T cells isolated from tumors,
## spleens, or MFPs of EO771.ova-bearing, LCMV-immune mice (n = 9).

##==============================================
## From GEO entry:

## Cell barcodes and unique molecular identifiers (UMIs) associated
## with the aligned reads were corrected and filtered. Filtered
## gene-barcode matrices containing only barcodes with UMI counts
## passing threshold for cell detection were imported to Seurat
## v4.0.

## Basic QC: Cells with a high percentage of UMIs from mitochondrial genes
## (> 5% of a cell’s total UMI count) and expressing < 1000 genes were
## removed. Cells with an RNA library size between 1,000 and 20,000
## UMI counts and a protein library size between 100 and 2,500 UMI
## counts were retained.


## Contaminating B cells (uniquely expressing
## high levels of Cd79a and Cd19) and macrophages (uniquely expressing
## high levels of Csf2ra) were removed. After contaminating cells
## (non-CD8+ T cells) were removed, T-cell receptor (TCR) genes were
## removed so they would not influence subsequent analyses.

## A list of ## murine TCR genes were downloaded from IMGT (The
## International ## ImMunoGeneTics Information website,
## http://www.imgt.org) and converted to gene symbols with SynGo.

## The filtered and integrated Seurat object was converted to an
## AnnData object for analyses with totalVI (Total Variational
## Inference) and scVI-tools (single- cell Variational Inference
## tools).


## To identify P14 (Thy1.1+), OT-I (CD45.1+), and endogenous
## polyclonal CD8+ (Thy1.1− and CD45.1−) T cells by protein expression
## of Thy1.1 and CD45.1, we used totalVI, which models the data to the
## calculated ‘denoised’ protein expression values that have been
## corrected for background noise (60). Before modeling with totalVI,
## an initial gene filter removed genes expressed in fewer than four
## cells, and the top 6000 highly variable genes (HVGs) were selected
## by the Seurat v3 method (56) as implemented by scVI (60, 61). We
## ran totalVI modeling with a 20-dimensional latent space and default
## parameters to produce denoised protein data for all proteins, and
## each 10× lane was treated as a batch. After modeling, the totalVI
## denoised expression values for Thy1.1 and CD45.1 were used to
## determine the specificity of the CD8+ T cells. Within the
## endogenous, polyclonal CD8+ T cell population, contaminating naïve
## CD8+ T cells were also removed. Specifically, endogenous polyclonal
## CD8+ T cells derived from the spleen that expressed both CD103 and
## CD62L were filtered out (log2 denoised expression value > 2 for
## both proteins). After determining cell specificity and filtering
## out endogenous naïve polyclonal CD8+ T cells, the raw expression
## values for Thy1.1 and CD45.1 were removed, and totalVI was used to
## model the data again and determine protein expression values for
## all other proteins (see table S1). The denoised protein expression
## values for each cell were saved as metadata observations for
## downstream analyses. After filtering, the CITE-seq dataset
## contained a total of 7754 OT-I; 8781 P14; and 7686 CD44hi
## polyclonal, endogenous CD8+ T cells derived from the spleen, tumor,
## and MFP.

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratExtend)
library(cowplot)
library(RColorBrewer)

mycolors <- c(AntiTumor=rgb(79,198,0, maxColorValue = 255), ## green
             AntiViral=rgb(90,1,7, maxColorValue = 255), ## dark brown
             MFP      =rgb(162,0,88, maxColorValue = 255), ## raspberry
             SPL      =rgb(183,151,96, maxColorValue = 255), ## coffee
             TUMOR    =rgb(59,92,255, maxColorValue = 255) ## blueberry
             )

getObject <- function(fname, pname) {

    set.seed(1234)
    
    ##============================
    ## RNA
    ##============================
    data <- Read10X(data.dir=fname)
    
    print(names(data))
    ##[1] "Gene Expression"  "Antibody Capture" "Custom"
    print(rownames(data$"Antibody Capture"))
    ## [1] "LY6C"   "CD62L"  "CD103"  "OTI"    "CD69"   "CD127"  "KLRG1"  "P14"
    ## [9] "CX3CR1" "LY108"  "PD1"    "CD39"   "CD49a"  "P2RX7"  "TIM3"   "CXCR6"
    ##[17] "TIGIT"  "CD38"
    print(rownames(data$"Custom"))
    ##[1] "SPL"   "MFP"   "TUMOR"
    
    ## P14 is anti-viral
    ## OT1 is anti-tumor
    
    
    obj <- CreateSeuratObject(counts = data[["Gene Expression"]],
                             min.cells = 3, min.features = 200,
                             project = pname)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern ="^mt-")
    
    ##============================
    ## HTO hashes - SAMPLE
    ##============================
    
    obj[['sample']] <- CreateAssayObject(counts = data[["Custom"]][,colnames(obj)])
    obj <- NormalizeData(obj, assay = "sample", normalization.method = "CLR")
    obj <- HTODemux(object = obj, assay = "sample",positive.quantile = 0.99)
    Idents(obj) <- obj$hash.ID
    obj <- subset(obj, subset = sample_classification.global=='Singlet')
    obj$Tissue <- obj$sample_classification
    
    ##============================
    ## HTO hashes - PROTEINS
    ##============================
    
    obj[['PROT']] <- CreateAssayObject(counts = data$`Antibody Capture`[,colnames(obj)])
    obj <- NormalizeData(obj, assay = "PROT", normalization.method = "CLR")
    obj <- HTODemux(object = obj, assay = "PROT", positive.quantile = 0.99)
    
    ##============================
    ## BASIC QC AS IN PAPER
    ##============================
    obj <- subset(obj, subset =
                           nFeature_RNA > 1000 &
                           nCount_RNA>1000 &
                           nCount_RNA<20000 &
                           nCount_PROT>100 &
                           nCount_PROT<2500 &
                       percent.mt < 5)
    
    return(obj)
}



##=================================================
## PROCESS GSE231541-A AND GSE231541-B
##=================================================
obj.list =
    lapply(c('A', 'B'),
           function(x) {
               getObject(fname = sprintf('data/GSE231541-%s',x), pname=sprintf('Lane %s',x))
           })
sapply(obj.list, dim)
##           [,1]  [,2]
## [1,] 19782 20183
## [2,] 13047 16240

obj = merge(x=obj.list[[1]],  y=obj.list[-1])
rm(obj.list); gc()
obj
## An object of class Seurat 
## 20780 features across 29287 samples within 3 assays 
## Active assay: RNA (20759 features, 0 variable features)
##  2 layers present: counts.Lane A, counts.Lane B
##  2 other assays present: sample, PROT


##=====================================
## Filtering as in original manuscript:
##=====================================
## Remove B cells (uniquely expressing high levels of Cd79a and Cd19)
## Remove macrophages (uniquely expressing high levels of Csf2ra)
## Remove non-CD8+ T cells
## Remove TCR genes so they would not influence subsequent analyses.
obj <- subset(obj, subset = prot_OTI > 1 | prot_P14 > 1) ## only anti-viral or anti-tumor
obj <- NormalizeData(obj)  ## need data slot for the subset= options
obj <- subset(obj, subset = Cd79a > 0 | Cd19 > 0, invert=TRUE) ## remove B cells
obj <- subset(obj, subset = Csf2ra > 0, invert=TRUE) ## remove macrophages
obj <- subset(obj, subset = Cd8b1 > 0 | Cd8a > 0) ## only Cd8b1
obj$Specificity.AntiTumor <- obj[["PROT"]]@counts["OTI",]>0
obj$Specificity.AntiVirus <- obj[["PROT"]]@counts["P14",]>0
table(Tum=obj$Specificity.AntiTumor, Vir=obj$Specificity.AntiVirus)
##        Vir
## Tum     FALSE TRUE
##   FALSE     0 4472
##   TRUE   4875 3139
## Take out double positives to be most conservative
obj <- subset(obj, subset= Specificity.AntiTumor & Specificity.AntiVirus, invert=TRUE) 
table(Tum=obj$Specificity.AntiTumor, Vir=obj$Specificity.AntiVirus)
##        Vir
## Tum     FALSE TRUE
##   FALSE     0 4472
##   TRUE   4875    0
obj$Specificity <- rep('anti-Tumor', ncol(obj))
obj$Specificity[obj$Specificity.AntiVirus] <- 'anti-Virus'

## 'Remove' TCR genes effectively by setting them equal to zero 
tcr.pattern  <- "^(?i)tr[abdg][avjc]"   # ignores case; drops TCR V/J genes (e.g., Trbv12-1)
tcr.genes <-  grep(tcr.pattern, rownames(obj), value=TRUE)
inA <- intersect(tcr.genes, rownames(obj[["RNA"]]$"counts.Lane A"))
inB <- intersect(tcr.genes, rownames(obj[["RNA"]]$"counts.Lane B"))
obj[["RNA"]]$"counts.Lane A"[inA,] <- 0
obj[["RNA"]]$"data.Lane A"[inA,] <- 0
obj[["RNA"]]$"counts.Lane B"[inB,] <- 0
obj[["RNA"]]$"data.Lane B"[inB,] <- 0

##=====================================
## Standard Pipeline to visualize as in original paper
##=====================================

obj <- NormalizeData(obj, assay='RNA')
obj <- FindVariableFeatures(obj,assay='RNA')
obj <- ScaleData(obj, assay='RNA', vars.to.regress=c("percent.mt"))
obj <- RunPCA(obj, npcs=30)
obj <- FindNeighbors(obj, assay='RNA', dims=1:30, reduction='pca')
obj <- FindClusters(obj, assay='RNA', resolution=0.8)
obj <- RunUMAP(obj, assay='RNA', dims=1:30, reduction='pca')

X <- DimPlot(obj, group.by='Tissue', cols=mycolors[3:5]) + 
    theme(aspect.ratio=1) + coord_flip() + scale_x_reverse() + scale_y_reverse() +
    ggtitle('Location') +
    scale_color_manual(labels = c("MFP" = "MFP", "SPL" = "Spleen", "TUMOR" = "Tumor"),
                       values=mycolors[3:5]) +
    theme(
        axis.title = element_blank(), 
        axis.text = element_blank(),  
        axis.ticks = element_blank(), 
        axis.line = element_blank(),  
        legend.text = element_text(size = 16),
        legend.title = element_blank()
        ) +
    guides(color = guide_legend(override.aes = list(size = 5)))

Y <- DimPlot(obj, group.by='Specificity', cols=unname(mycolors[1:2])) +
    theme(aspect.ratio=1) + coord_flip() + scale_x_reverse() + scale_y_reverse() +
    ggtitle('Specificity') +
    theme(
        axis.title = element_blank(), 
        axis.text = element_blank(),  
        axis.ticks = element_blank(), 
        axis.line = element_blank(),  
        legend.text = element_text(size = 16),
        legend.title = element_blank(), 
        ) +
    guides(color = guide_legend(override.aes = list(size = 5)))


X | Y

##====================================================
## Compare to Kleist OT-I data (X6)
##====================================================



##====================
## Rename as masopust data
masopust <- obj
table(masopust$Tissue, masopust$Specificity)
##       anti-Tumor anti-Virus
## MFP            8        577
## SPL         3153       2038
## TUMOR       1714       1857
masopust$FigClass <- paste(as.character(masopust$Specificity), as.character(masopust$Tissue))
data.frame(table(masopust$FigClass))
##               Var1 Freq
## 1   anti-Tumor MFP    8
## 2   anti-Tumor SPL 3153
## 3 anti-Tumor TUMOR 1714
## 4   anti-Virus MFP  577
## 5   anti-Virus SPL 2038
## 6 anti-Virus TUMOR 1857

## Create simplified merged object to facilitate side by side (add RNA
## normalization)
masopust2 <- masopust
masopust2@meta.data <- masopust2@meta.data[,c('FigClass'), drop=FALSE]
DefaultAssay(masopust2) <- 'RNA'
masopust2 <- DietSeurat(masopust2, assays='RNA', layers='counts')
gc()

##====================
## Read in Kleist data
obj <- readRDS('seurat_obj.mouseOT1.Kleist2026.RDS')
dotplotOrder <- as.character(c( 7,2,5,6,  1,3,4,0))  ## dots look best in Mem T-cell
obj$ClusterNums <- factor(obj$ClusterNums, levels=dotplotOrder)

## Create simplified merged object to facilitate side by side (add RNA
## normalization)
obj2 <- obj
obj2$FigClass <- as.character(obj2$orig.ident)
obj2@meta.data <- obj2@meta.data[,c('FigClass', 'ClusterNums'), drop=FALSE]
DefaultAssay(obj2) <- 'RNA'
obj2 <- DietSeurat(obj2, assays='RNA', layers='counts')
gc()

##====================
## Merge data and subset to just tissue types we want
obj2 <- merge(obj2, masopust2)
data.frame(table(obj2$FigClass))
##               Var1 Freq
## 1   anti-Tumor MFP    8
## 2   anti-Tumor SPL 3153
## 3 anti-Tumor TUMOR 1714
## 4   anti-Virus MFP  577
## 5   anti-Virus SPL 2038
## 6 anti-Virus TUMOR 1857
## 7     IntactSpleen 1737
## 8  TcircmOnlyTumor 1400
## 9     TrmOnlyTumor 1540

mylevels <- c('TcircmOnlyTumor', 'TrmOnlyTumor',
             'anti-Virus TUMOR', 'anti-Tumor TUMOR', 'anti-Virus MFP')
obj2 <- subset(obj2, subset=FigClass %in% mylevels)
obj2$FigClass <- factor(obj2$FigClass, levels=mylevels)


##====================
## Log1p normalization on RNA counts so normalizes w/in sample
DefaultAssay(obj2) <- 'RNA'
obj2 <- JoinLayers(obj2)  ## join SeuratProject and MouseOT1
obj2[["RNA"]] <- split(obj2[["RNA"]], f=obj2$FigClass)  ## split by figClass instead
obj2 <- NormalizeData(obj2)  ## only need 'data' slot for DotPlot, not scale.data
obj2 <- JoinLayers(obj2)


##====================
## Breakout Kleist data into subclusters in TumorTrm-rich clusters or Spleen-rich clustrs
obj3 <- obj2
obj3$CircOrNot <- as.character(obj3$FigClass)
obj3$CircOrNot[which(obj3$ClusterNums==2)] =
    paste0(obj3$CircOrNot[which(obj3$ClusterNums==2)], '_TumorTcm')
where <- which(obj3$ClusterNums %in% c(1,3,4,0))
obj3$CircOrNot[where] <- paste0(obj3$CircOrNot[where], '_TumorTrm')
where <- which(obj3$ClusterNums %in% c(7,5,6))
obj3$CircOrNot[where] <- paste0(obj3$CircOrNot[where], '_Spleen')

## Ignore Spleen-rich clusters
mylevs <- c(
    'TcircmOnlyTumor_TumorTcm',
    'TrmOnlyTumor_TumorTcm',
    'TcircmOnlyTumor_TumorTrm',
    'TrmOnlyTumor_TumorTrm',
    'anti-Virus TUMOR',
    'anti-Tumor TUMOR',
    'anti-Virus MFP'
)
obj3 <- subset(obj3, subset=CircOrNot %in% mylevs)
obj3$CircOrNot <- factor(obj3$CircOrNot, levels=mylevs)


##==================================================
## Show Residency and Circulating signature genes across samples
##==================================================

## Memory t-cell subsets dotplot
features.ResCirc <- c(  ## bottom to top
    ##~~~Exhaustion  
    'Lag3',
    'Ctla4',
    'Pdcd1',
    'Tigit',
    'Entpd1',
    'Havcr2',
    'Cd244a',
    'Slamf6',
    'Tox',
    ##~~~Residency
    'Zfp683', ## Zfp683 ==='Hobit',
    'Itgae',  
    'Zap70',
    'Cxcr6',
    'Il2ra',  
    'Prdm1', 
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

lineIdx <- c(grep('Cd69', features.ResCirc) + 0.5,
            grep('Tox', features.ResCirc) + 0.5)

dp.ResCircGSE =
    DotPlot(obj3, assay = 'RNA',
            group.by = 'CircOrNot',
            features = features.ResCirc,
            cluster.idents = FALSE, cols = 'RdBu',
            col.max = 1.5, col.min = -1.5
            ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.7, size = 14))+  ## identities vjust LR, hjust TD
    theme(axis.text.y = element_text(face = "italic", size = 14)) +  ## gene names
    coord_flip() +
    geom_vline(xintercept=lineIdx,
               col=rgb(36,122,190,maxColorValue=255),
               linetype='dotted', linewidth=0.7) +
    ggtitle('Memory T cell subsets') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 1, r = 1, b = 1, l = 3, unit = "cm")
          ) +
    geom_hline(yintercept=4.5,
               col=rgb(36,122,190,maxColorValue=255),
               linetype='solid', linewidth=0.7) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank()
          )



ggsave("PLOT_GSE231541_DotplotResCirc_breakoutTcmTrm.png", dp.ResCircGSE,
       width=7, height=10, units="in")



##==============================================
sessionInfo()

## R version 4.4.3 (2025-02-28)
## Platform: x86_64-conda-linux-gnu
## Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

## Matrix products: default
## BLAS/LAPACK: R443/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

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
## [1] stats     graphics  grDevices utils     datasets  methods   base     

## other attached packages:
##  [1] RColorBrewer_1.1-3     cowplot_1.2.0          SeuratExtend_1.2.4    
##  [4] SeuratExtendData_0.3.0 dplyr_1.1.4            patchwork_1.3.2       
##  [7] ggplot2_4.0.0          Seurat_5.3.1.9999      SeuratObject_5.2.0    
## [10] sp_2.2-0              

## loaded via a namespace (and not attached):
##  [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3         
##  [4] rlang_1.1.6            magrittr_2.0.4         RcppAnnoy_0.0.22      
##  [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7        
## [10] compiler_4.4.3         spatstat.geom_3.6-0    png_0.1-8             
## [13] vctrs_0.6.5            reshape2_1.4.4         stringr_1.6.0         
## [16] pkgconfig_2.0.3        fastmap_1.2.0          promises_1.5.0        
## [19] purrr_1.2.0            jsonlite_2.0.0         goftest_1.2-3         
## [22] later_1.4.4            spatstat.utils_3.2-0   irlba_2.3.5.1         
## [25] parallel_4.4.3         cluster_2.1.8.1        R6_2.6.1              
## [28] ica_1.0-3              stringi_1.8.7          spatstat.data_3.1-9   
## [31] reticulate_1.44.0      parallelly_1.45.1      spatstat.univar_3.1-4 
## [34] lmtest_0.9-40          scattermore_1.2        Rcpp_1.1.0            
## [37] tensor_1.5.1           future.apply_1.20.0    zoo_1.8-14            
## [40] sctransform_0.4.2      httpuv_1.6.16          Matrix_1.7-4          
## [43] splines_4.4.3          igraph_2.2.1           tidyselect_1.2.1      
## [46] abind_1.4-8            codetools_0.2-20       spatstat.random_3.4-2 
## [49] miniUI_0.1.2           spatstat.explore_3.5-3 listenv_0.10.0        
## [52] lattice_0.22-7         tibble_3.3.0           plyr_1.8.9            
## [55] withr_3.0.2            shiny_1.11.1           S7_0.2.0              
## [58] ROCR_1.0-11            Rtsne_0.17             future_1.67.0         
## [61] fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7       
## [64] fitdistrplus_1.2-4     pillar_1.11.1          KernSmooth_2.23-26    
## [67] plotly_4.11.0          generics_0.1.4         RcppHNSW_0.6.0        
## [70] scales_1.4.0           globals_0.18.0         xtable_1.8-4          
## [73] glue_1.8.0             lazyeval_0.2.2         tools_4.4.3           
## [76] data.table_1.17.8      RSpectra_0.16-2        RANN_2.6.2            
## [79] dotCall64_1.2          grid_4.4.3             tidyr_1.3.1           
## [82] nlme_3.1-168           cli_3.6.5              spatstat.sparse_3.1-0 
## [85] spam_2.11-1            viridisLite_0.4.2      uwot_0.2.4            
## [88] gtable_0.3.6           digest_0.6.37          progressr_0.18.0      
## [91] ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2          
## [94] htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7            
## [97] mime_0.13              MASS_7.3-65           
 

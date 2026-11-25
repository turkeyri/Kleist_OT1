## Sonia Leach
## 11/25/2026

## Code to reproduce results for Kleist et al, 2027 manuscript with
## human GBM data

## (see object details and sessionInfo() at end)

##===============================
## Load Libraries
##===============================

library(scRepertoire)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scales)
library(tidyr)
library(ragg)  ## cairo_ps and sans fonts
library(ggVennDiagram)
library(ggpubr)
library(rstatix)


## set.seed(1234)
DATADIR <- 'sample_filtered_feature_bc_matrix'


##====================
## RNA object
##====================
data <- Read10X(data.dir = DATADIR)
names(data)
## [1] "Gene Expression"  "Antibody Capture"
rownames(data[['Antibody Capture']])
## [1] "Blood"     "Tumor"     "Adj-Brian" "COVID"     "EBV"       "Flu"
## [7] "CMV"       "CD69.1"    "CD103"

data[["Antibody Capture"]]@Dimnames[[1]][3]<-"Adj-Brain" # fix typo

obj <- CreateSeuratObject(counts = data[["Gene Expression"]],
                          min.cells = 3, min.features = 200,
                          project = 'HumanGBM')
dim(obj) ## [1] 24636  6458

obj[["percent.MT"]] <- PercentageFeatureSet(obj, pattern ="^MT-")
obj <- subset(obj, subset =
                       nFeature_RNA > 500 &
                       nFeature_RNA < 4000 &
                       percent.MT <10 &
                       percent.MT > 0 & ## 12 have percent.MT=0
                       nCount_RNA < 10000)
dim(obj) ## [1] 24636  5487


##=============================
## Hash tags
##=============================
obj[['tissue']] <- CreateAssayObject(
    counts = data$`Antibody Capture`[1:3,colnames(obj)])
obj <- NormalizeData(obj, assay = "tissue", normalization.method = "CLR")
obj <- HTODemux(object = obj,  assay = 'tissue', positive.quantile = 0.95)
## Cutoff for Blood : 5 reads
## Cutoff for Tumor : 109 reads  ## v453 106
## Cutoff for Adj-Brain : 112 reads ## v453 104

obj <- subset(obj,
              subset = tissue_classification %in%
                  c('Blood', 'Tumor','Adj-Brain'))
table(obj$tissue_classification)
##Adj-Brain     Blood     Tumor
##     1965       671      2128
dim(obj) ## [1] 24636  4764
obj$orig.ident <- obj$tissue_classification
Idents(obj) <- 'orig.ident'

##=========================================
## RNA-seq analysis
##=========================================

obj <- SCTransform(obj, vst.flavor = "v2", vars.to.regress = c("percent.MT"),
                   return.only.var.genes=FALSE, verbose = TRUE)

## Since are only CD8 t-cells, remove TCR Genes from
## VariableFeatures(obj) so fluctuations not just because of clonotype
markers.remove <- grep(pattern = "^TRAV|^TRBV|^TRGV|^TRDV|^TRBC|^TRAC|^TRDC",
                       rownames(obj), value=TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj), markers.remove)
obj <- RunPCA(obj, assay = 'SCT')
obj <- RunUMAP(obj, dims = 1:30)
obj <- FindNeighbors(obj, dims = 1:30, assay = "SCT")
obj <- FindClusters(obj, resolution = 0.3)
Idents(obj) <- "SCT_snn_res.0.3"


(DimPlot(obj, label = FALSE, group.by = 'tissue_classification') +
 DimPlot(obj, label = TRUE) ) / 
 DimPlot(obj, label = TRUE, split.by = 'tissue_classification') 


table(Idents(obj))
##    0    1    2    3    4    5    6    7    8    9
## 1133  864  862  647  481  329  153  128  100   67


##=========================================
## TCR contig info
##=========================================

## Were there any collisions of barcodes (so had to change -1 to -2,
## etc?) Ans. no
table(suffix=gsub('.+-', '', rownames(obj[[]])),
      tissue=obj$tissue_classification)
## suffix   Adj-Brain Blood Tumor
##      1        1965   671  2128

## Subset contig info by cells kept in obj
contigs <- read.csv('filtered_contig_annotations.csv')
myBarcodes <-  split(colnames(obj), obj$tissue_classification) ## split by samp
contig.list <-  lapply(myBarcodes,
                         function(x) {
                             contigs[contigs$barcode %in% x,]
                         })
names(contig.list) <- names(myBarcodes)
sapply(contig.list, dim)
##              Adj-Brain Blood Tumor
## [1,]      1767   785  2247
## [2,]        31    31    31

combined <- combineTCR(contig.list,
                       samples = names(contig.list),  ## need for w/ expr
                       removeNA = FALSE,
                       removeMulti = FALSE,
                       filterMulti = FALSE,
                       filter.nonproductive = TRUE
                       )
## Save info about all clones
## Stack lists in combined [can separate by combined$sample]
flattened <- do.call('rbind',
                    combined)[, c('barcode', 'sample',
                                  'cdr3_nt2', 'cdr3_aa2')]
rownames(flattened) <- c()
flattened$sample <- factor(flattened$sample,
                          levels=c('Blood', 'Tumor', 'Adj-Brain'))

unique.nt.aa <- unique(flattened[,c("cdr3_nt2", "cdr3_aa2")])
unique.nt.aa <- unique.nt.aa[!is.na(unique.nt.aa[,1]),]
rownames(unique.nt.aa) <- unique.nt.aa[,1]


## Which clones are in which samples?
## verified this is what clonalcompare() does to compute proportion
tt <- table(flattened$cdr3_nt2, flattened$sample)
tt <- tt[order(tt[,1], tt[,2], tt[,3]),]
tt.p <- sweep(tt, MARGIN=2, STAT=apply(tt,2,sum), FUN='/')

## Annotate tally with nt and aa info
CloneInfo <- data.frame(as.data.frame.matrix(tt),
                    as.data.frame.matrix(tt.p),
                    unique.nt.aa[rownames(tt),]
                    )
colnames(CloneInfo) <- c(
    "Blood",                "Tumor",                "Adj.Brain",
    "Blood_proportion",     "Tumor_proportion", "Adj.Brain_proportion",
    "original.clones",      "original.clones.aa")
CloneInfo$WhereString <-
    gsub('\\.$', '',
         gsub('^\\.', '',
              gsub('\\.\\.', '.',
                   apply(cbind(ifelse(CloneInfo$Blood>0, 'Blood', ''),
                               ifelse(CloneInfo$Tumor>0, 'Tumor', ''),
                               ifelse(CloneInfo$Adj.Brain>0, 'Brain', '')),
                         1, paste0, collapse='.'))))
## Are there shared clones?
cbind(table(CloneInfo$WhereString))
##                   [,1]
## Blood              195
## Blood.Brain          6
## Blood.Tumor         18
## Blood.Tumor.Brain   36
## Brain              153
## Tumor              399
## Tumor.Brain         43

## What does this CloneInfo look like (minus the rownames)
x <- rbind(head(CloneInfo,2), tail(CloneInfo,2)) ; rownames(x) = c(); x
##   Blood Tumor Adj.Brain Blood_proportion Tumor_proportion Adj.Brain_proportion
## 1     0     0         1       0.00000000      0.000000000          0.001004016
## 2     0     0         1       0.00000000      0.000000000          0.001004016
## 3    31    15        24       0.07045455      0.013215859          0.024096386
## 4    47     7        13       0.10681818      0.006167401          0.013052209
##                                 original.clones original.clones.aa
## 1 TGCAGCGTCGCACCTAGCGGGAGCGGGGGAGATACGCAGTATTTT    CSVAPSGSGGDTQYF
## 2          TGCAGCGTGGTCACCAGCTCCTACGAGCAGTACTTC       CSVVTSSYEQYF
## 3       TGCGCCAGCAGCCAAGATCGTAACACTGAAGCTTTCTTT      CASSQDRNTEAFF
## 4 TGTGCCAGCAGTGAGACTAGCGGGGGATACAATGAGCAGTTCTTC    CASSETSGGYNEQFF
##         WhereString
## 1             Brain
## 2             Brain
## 3 Blood.Tumor.Brain
## 4 Blood.Tumor.Brain


## Before combine TCR with expression, strip barcode off each combined
combined2 <- lapply(combined,
                   function(x) {
                       x$barcode = gsub('.+\\_', '', x$barcode)
                       return(x)
                       })
obj <- combineExpression(combined2,
                         obj,
                         chain='TRB',
                         cloneCall="nt",
                         group.by = "sample")
table(obj$tissue_classification, !is.na(obj$clonalFrequency))
##           FALSE TRUE
## Adj-Brain   969  996
## Blood       231  440
## Tumor       993 1135

## Limit obj to just those that have TCR seq data
obj <- subset(obj, subset=!is.na(clonalFrequency))
table(obj$tissue_classification, !is.na(obj$clonalFrequency))
##           TRUE
## Adj-Brain  996
## Blood      440
## Tumor     1135

## Add clone info as metadata
obj$CDR3TCRB=gsub('.+\\_', '', obj$CTnt)
table(obj$CDR3TCRB %in% rownames(CloneInfo), useNA='ifany') ##  TRUE 2571
obj <- AddMetaData(obj, metadata=CloneInfo[obj$CDR3TCRB,])


##===================================================================
##    Visualizations
##===================================================================
Idents(obj) <- 'SCT_snn_res.0.3'

##===================================================================
## Alluvial plots depicting overlap of the top 10 CD8+ T cell
## clonotypes across Patient 018 blood, high grade glioma tumor, and
## peritumoral brain, with flow width proportional to clone frequency.
##===================================================================

Top10Table <- clonalCompare(combined, top.clones=10,
                 clone.call='nt',
                 chain = 'TRB',
                 relabel.clones=TRUE,
                 palette= 'viridis',
                 samples=c('Blood', 'Tumor', 'Adj-Brain'),
                 order.by=c('Blood', 'Tumor', 'Adj-Brain'),
                 graph = 'alluvial',
                 exportTable=TRUE
                 )

## Reorder clones so that in alluvial plot go in an order that
## emphasizes those common to all tissues, then Tumor shared with 1
## tissue, then only Tumor

## clonalCompare will reorder based on appearence in tallest bar
## (Adj.Brain) then next bar, then next bar, so have to find overlap
## among adjacent
where <- factor(
    sapply(split(as.character(Top10Table$Sample),
                 as.character(Top10Table$original.clones)),
           paste, collapse=':'),
    levels=c(
        "Adj-Brain:Blood:Tumor",
        "Adj-Brain:Tumor",
        "Blood:Tumor",
        "Tumor"
    ))
myclonesInOrder <- names(sort(where))

## clonalCompare orders by alphanumeric sort so trick it by adding
## numbers in order we want clone names
X <- do.call('rbind', combined)
cdr3 <- unique(X$cdr3_nt2)
bigclonesInOrder <- c(myclonesInOrder, setdiff(cdr3, myclonesInOrder))
newNames <- setNames(
    paste0(sprintf('%02d', 1:length(bigclonesInOrder)), bigclonesInOrder),
    bigclonesInOrder)
newmyclonesInOrder <- newNames[myclonesInOrder]

## Change order in combined object to force our ordering of factors
combined.factorOrder <- lapply(
    combined,
    function(x) {
        ## Clonal compare uses CTnt when clone.call == 'nt'
        CTnt = as.character(x$CTnt)
        CTnt.B = gsub('.+\\_([^\\_]+)', '\\1',CTnt)
        CTnt.A = gsub('(.+\\_)[^\\_]+', '\\1',CTnt)
        x$CTnt = paste0(CTnt.A, newNames[CTnt.B])
        # Error in strsplit(dat[, cloneCall], "_") : non-character argument
        return(x)
    })

## Colors chosen to emphasize clones not among all 3
alluvial.colors <- 
    c(darkblue   = "#1F77B4FF", lightblue  = "#AEC7E8FF",
      darkpink   = "#E377C2FF", lightpink  = "#F7B6D2FF",
      darkgreen  = "#2CA02CFF", lightgreen = "#98DF8AFF",
      darkpurple = "#9467BDFF", lightpurple= "#C5B0D5FF",
      darkbrown  = "#8C564BFF", lightbrown = "#C49C94FF",
      darkteal   = "#499894FF", lightteal  = "#86BCB6FF",
      darkgrey   = "#7F7F7FFF", lightgrey  = "#C7C7C7FF",
      darkyellow = "#BCBD22FF", lightyellow= "#DBDB8DFF",
      darkturq   = "#17BECFFF", lightturq  = "#9EDAE5FF",
      lightorange= "#FFBB78FF", darkorange = "#FF7F0EFF",
      darkred    = "#D62728FF", lightred   = "#FF9896FF"
      )

p.alluvial <- 
    clonalCompare(combined.factorOrder,
              top.clones=10,
              clone.call='nt',
              chain = 'TRB',
              relabel.clones=FALSE,
              palette= 'viridis', ## 'inferno' ## default
              samples=c('Blood', 'Tumor', 'Adj-Brain'),
              order.by=c('Blood', 'Tumor', 'Adj-Brain'),
              graph = 'alluvial',
              exportTable=FALSE
              ) + ggtitle('CD8+ CDR3 TRB nt\nTop 10 Clones') +
    scale_fill_manual(values = unname(alluvial.colors)) +
    theme_classic(base_size=16) + NoLegend()


ggsave(filename="PLOT_Alluvial_CD8.CDR3.TRB.nt.png", p.alluvial,
       width=6, height=8, units="in")
ggsave(filename="PLOT_Alluvial_CD8.CDR3.TRB.nt.eps", p.alluvial,
       width=6, height=8, units="in", dpi=300, device=cairo_ps)

##===================================================================
## Venn diagram showing all shared TCR clonotypes across Patient
## 018 tissues.
##===================================================================

X <- unique(obj@meta.data[,c('tissue_classification',
                            'original.clones',
                            'WhereString')])
cbind(table(X$WhereString))
##                   [,1]
## Blood              195
## Blood.Brain         12
## Blood.Tumor         36
## Blood.Tumor.Brain  108
## Brain              153
## Tumor              399
## Tumor.Brain         86

x <- split(as.character(X$original.clones),
           as.character(X$tissue_classification))
sapply(x, length)
## Adj-Brain     Blood     Tumor
##       238       255       496

p.venn <- ggVennDiagram(x) +
    scale_fill_gradient(low="#D1E5F0FF", high="#2166ACFF")

ggsave(filename="PLOT_Venn_CD8.CDR3.TRB.nt.png", p.venn,
       width=6, height=6, units="in")
ggsave(filename="PLOT_Venn_CD8.CDR3.TRB.nt.eps", p.venn,
       width=6, height=6, units="in", dpi=300, device= cairo_ps)

##===================================================================
## Violin plots comparing GZMA expression in tumor-infiltrating
## CD8+ T cells belonging to blood-tumor shared clonotypes
## vs. peritumoral brain shared clonotypes.
##===================================================================

vlnColors <- c(gold =rgb(228,169,70, maxColorValue=255),
              slate=rgb(97,141,183, maxColorValue=255))

## Get expression of gene as metadata for cells with single Tumor-shared clones
gene <- 'GZMA'
temp <- subset(obj, subset=WhereString %in% c('Blood.Tumor', 'Tumor.Brain'))
Idents(temp) <- 'WhereString'
temp$GENE <- temp[["SCT"]]$data[gene,]
ctvals <- temp@meta.data[,c('WhereString', 'GENE')]


library(ggpubr)
library(rstatix)

pwc <- ctvals %>%
    pairwise_wilcox_test(as.formula('GENE ~ WhereString'),
                         p.adjust.method="BH") %>%
    add_xy_position(x <- "WhereString")  ## already has significance
pwc$y.position <- 1.1 * pwc$y.position  ## was too squished

p.vln <- 
    ggplot(ctvals, aes(x=WhereString, y=GENE)) +
    theme_classic(base_size=20) +
    geom_violin(aes(fill=WhereString)) +
    stat_summary(fun="median", geom = "crossbar",
                 width = .5, color = "grey30") +
    stat_pvalue_manual(pwc, hide.ns = TRUE, label = "p", label.size=6) +
    ggtitle(gene) + xlab('') +
    ylab('Expression Level')  +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)
          ) +
    scale_fill_manual(values = alpha(unname(vlnColors), alpha=0.8))


ggsave(filename="PLOT_VlnGZMA_CD8.CDR3.TRB.nt.png", p.vln,
       width=6, height=8, units="in")
ggsave(filename="PLOT_VlnGZMA_CD8.CDR3.TRB.nt.eps", p.vln,
       width=6, height=8, units="in", dpi=300, device= cairo_ps)

##===================================================================
## Save processed data for GEO supplement
##===================================================================
obj
## An object of class Seurat
## 45007 features across 2571 samples within 3 assays
## Active assay: SCT (20368 features, 2921 variable features)
##  3 layers present: counts, data, scale.data
##  2 other assays present: RNA, tissue
##  2 dimensional reductions calculated: pca, umap

colnames(obj@meta.data)
##  [1] "orig.ident"                   "nCount_RNA"
##  [3] "nFeature_RNA"                 "percent.MT"
##  [5] "nCount_tissue"                "nFeature_tissue"
##  [7] "tissue_maxID"                 "tissue_secondID"
##  [9] "tissue_margin"                "tissue_classification"
## [11] "tissue_classification.global" "hash.ID"
## [13] "nCount_SCT"                   "nFeature_SCT"
## [15] "SCT_snn_res.0.3"              "seurat_clusters"
## [17] "CTgene"                       "CTnt"
## [19] "CTaa"                         "CTstrict"
## [21] "clonalProportion"             "clonalFrequency"
## [23] "cloneSize"                    "CDR3TCRB"
## [25] "Blood"                        "Tumor"
## [27] "Adj.Brain"                    "Blood_proportion"
## [29] "Tumor_proportion"             "Adj.Brain_proportion"
## [31] "original.clones"              "original.clones.aa"
## [33] "WhereString"

names(combined)
##[1] "Adj-Brain" "Blood"     "Tumor"
sapply(combined, dim)
##     Adj-Brain Blood Tumor
##[1,]      1062   454  1209
##[2,]        12    12    12

colnames(combined[[1]])
## [1] "barcode"  "sample"   "TCR1"     "cdr3_aa1" "cdr3_nt1" "TCR2"
## [7] "cdr3_aa2" "cdr3_nt2" "CTgene"   "CTnt"   "CTaa"  "CTstrict"

saveRDS(obj, 'seurat_obj.humanGBM.Kleist2026.RDS')
write.csv(obj@meta.data, 'seurat_obj_metadata.humanGBM.Kleist2026.csv')
saveRDS(combined, 'seurat_TCRinfo_combined.humanGBM.Kleist2026.RDS')


##===================================================================
sessionInfo()

## R version 4.5.3 (2026-03-11)
## Platform: x86_64-conda-linux-gnu
## Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

## Matrix products: default
## BLAS/LAPACK: /R453/lib/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0

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
##  [1] future_1.70.0       rstatix_0.7.3       ggpubr_0.6.3       
##  [4] ggVennDiagram_1.5.7 ragg_1.5.2          tidyr_1.3.2        
##  [7] scales_1.4.0        patchwork_1.3.2     Seurat_5.4.0       
## [10] SeuratObject_5.4.0  sp_2.2-1            dplyr_1.2.1        
## [13] scRepertoire_2.7.3  ggplot2_4.0.3      

## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.23            splines_4.5.3              
##   [3] later_1.4.8                 tibble_3.3.1               
##   [5] R.oo_1.27.1                 polyclip_1.10-7            
##   [7] fastDummies_1.7.5           lifecycle_1.0.5            
##   [9] globals_0.19.1              lattice_0.22-9             
##  [11] MASS_7.3-65                 backports_1.5.1            
##  [13] magrittr_2.0.5              plotly_4.12.0              
##  [15] httpuv_1.6.17               otel_0.2.0                 
##  [17] glmGamPoi_1.22.0            sctransform_0.4.3          
##  [19] spam_2.11-3                 spatstat.sparse_3.1-0      
##  [21] reticulate_1.45.0           cowplot_1.2.0              
##  [23] pbapply_1.7-4               RColorBrewer_1.1-3         
##  [25] abind_1.4-8                 Rtsne_0.17                 
##  [27] GenomicRanges_1.62.1        purrr_1.2.2                
##  [29] R.utils_2.13.0              ggraph_2.2.2               
##  [31] BiocGenerics_0.56.0         tweenr_2.0.3               
##  [33] evmix_2.12                  IRanges_2.44.0             
##  [35] S4Vectors_0.48.1            ggrepel_0.9.8              
##  [37] irlba_2.3.7                 listenv_0.10.1             
##  [39] spatstat.utils_3.2-2        iNEXT_3.0.2                
##  [41] MatrixModels_0.5-4          goftest_1.2-3              
##  [43] RSpectra_0.16-2             spatstat.random_3.4-5      
##  [45] fitdistrplus_1.2-6          parallelly_1.47.0          
##  [47] DelayedMatrixStats_1.32.0   codetools_0.2-20           
##  [49] DelayedArray_0.36.1         ggforce_0.5.0              
##  [51] tidyselect_1.2.1            farver_2.1.2               
##  [53] viridis_0.6.5               matrixStats_1.5.0          
##  [55] stats4_4.5.3                spatstat.explore_3.8-0     
##  [57] Seqinfo_1.0.0               jsonlite_2.0.0             
##  [59] tidygraph_1.3.1             progressr_0.19.0           
##  [61] Formula_1.2-5               ggridges_0.5.7             
##  [63] ggalluvial_0.12.6           survival_3.8-6             
##  [65] systemfonts_1.3.2           tools_4.5.3                
##  [67] ica_1.0-3                   Rcpp_1.1.1-1               
##  [69] glue_1.8.1                  gridExtra_2.3              
##  [71] SparseArray_1.10.10         MatrixGenerics_1.22.0      
##  [73] withr_3.0.2                 fastmap_1.2.0              
##  [75] SparseM_1.84-2              digest_0.6.39              
##  [77] R6_2.6.1                    mime_0.13                  
##  [79] textshaping_1.0.5           scattermore_1.2            
##  [81] tensor_1.5.1                spatstat.data_3.1-9        
##  [83] R.methodsS3_1.8.2           generics_0.1.4             
##  [85] data.table_1.17.8           graphlayouts_1.2.3         
##  [87] httr_1.4.8                  htmlwidgets_1.6.4          
##  [89] S4Arrays_1.10.1             uwot_0.2.4                 
##  [91] pkgconfig_2.0.3             gtable_0.3.6               
##  [93] lmtest_0.9-40               S7_0.2.2                   
##  [95] SingleCellExperiment_1.32.0 XVector_0.50.0             
##  [97] htmltools_0.5.9             carData_3.0-6              
##  [99] dotCall64_1.2               Biobase_2.70.0             
## [101] png_0.1-9                   spatstat.univar_3.1-7      
## [103] ggdendro_0.2.0              reshape2_1.4.5             
## [105] rjson_0.2.23                nlme_3.1-169               
## [107] cachem_1.1.0                zoo_1.8-15                 
## [109] stringr_1.6.0               KernSmooth_2.23-26         
## [111] parallel_4.5.3              miniUI_0.1.2               
## [113] pillar_1.11.1               grid_4.5.3                 
## [115] vctrs_0.7.3                 RANN_2.6.2                 
## [117] promises_1.5.0              car_3.1-5                  
## [119] beachmat_2.26.0             xtable_1.8-8               
## [121] cluster_2.1.8.2             cli_3.6.6                  
## [123] compiler_4.5.3              rlang_1.2.0                
## [125] future.apply_1.20.2         ggsignif_0.6.4             
## [127] labeling_0.4.3              immApex_1.5.4              
## [129] plyr_1.8.9                  stringi_1.8.7              
## [131] viridisLite_0.4.3           deldir_2.0-4               
## [133] gsl_2.1-9                   lazyeval_0.2.3             
## [135] spatstat.geom_3.7-3         quantreg_6.1               
## [137] Matrix_1.7-4                RcppHNSW_0.6.0             
## [139] sparseMatrixStats_1.22.0    shiny_1.13.0               
## [141] SummarizedExperiment_1.40.0 ROCR_1.0-12                
## [143] igraph_2.1.4                broom_1.0.12               
## [145] memoise_2.0.1

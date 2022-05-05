
####======Enrichment analysis=============================PCOS DEP=================================

#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html

#====INSTALL PACKAGES================================================================================

# Lista de paquetes de funciones a instalar
.packages = c("BiocManager","devtools","dplyr","msigdbr","data.table","plyr","ggplot2", "reshape2","Hmisc","corrgram", 
              "psych","corrplot","faraway","tidyverse","tidyr","tibble","org.Hs.eg.db","emapplot","parallel","doParallel",
              "mixOmics","phyloseq","GOsummaries","ggupset","enrichplot","BiocGenerics",
              "survcomp","qusage","ReactomePA","clusterProfiler","DOSE","pathview",
              "MeSH.Hsa.eg.db","AnnotationHub","MeSHDbi","meshes","ggnewscale","svglite")
#"RFunrichWebService",
# Instala los paquetes sin? los tienes instalados
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  devtools::install_github("vqv/ggbiplot")
  devtools::install_github("fxia22/RadaFDR")
  BiocManager::install(c("mixOmics","phyloseq","GOsummaries","ReactomePA","enrichplot",'MeSHDbi',
                         "clusterProfiler","DOSE","pathview","viewPathway","MeSH.Hsa.eg.db",
                         "AnnotationHub","MeSHDbi","meshes","emapplot","msigdbr"),force = TRUE)
}

# Carga los paquetes sin? los tienes cargados
list.pk <- lapply(.packages, require, character.only=TRUE)
names(list.pk) <- .packages
list.pk
#=====Loading dataset===================================================================== 

DEP_prot <- as.data.frame(readxl::read_excel("data/20200229_FF_Dep_small_R2_abundances_Statistic_PCOSipp1_normalized in R1.xlsx")) # Open the data
rownames(DEP_prot) <- paste(DEP_prot$Accession,".",DEP_prot$`Gene`)
class(DEP_prot)


DEP_mean <- as.data.frame(readxl::read_excel("data/mean groups from perseu2_a.xlsx")) # Open the data
rownames(DEP_mean) <- paste(DEP_mean$Accession,".",DEP_mean$`Gene`)
class(DEP_mean)

main_DEP <- DEP_mean %>% dplyr::select(contains("F")) # Select columns whose names contains "Control" and "PCOS" (Intenstity columns)


Annotations <- as.data.frame(readxl::read_excel("data/Copy of PCOSprojectDataforLund (3).xlsx")) # Open the data
row.names(Annotations) <- Annotations$ID
head(Annotations)

row.Annotations <- as.data.frame(readxl::read_excel("data/row annotation_mean groups from perseu3.xlsx")) # Open the data
row.names(row.Annotations) <- paste(row.Annotations$protein,".",row.Annotations$GeneSymbol)

heatmap.clusters <- as.data.frame(readxl::read_excel("data/DEP_clusters.xlsx")) # Open the data
row.names(heatmap.clusters) <- paste(heatmap.clusters$Accession,".",heatmap.clusters$SYMBOL)
#========================================================================================

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

genes.DEP <- DEP_mean$Gene

df.ENTREZID = clusterProfiler::bitr(genes.DEP, fromType="SYMBOL", 
                                               toType="ENTREZID", 
                                               OrgDb="org.Hs.eg.db") #  convert biological IDs to Entrez gene ID.
#write.csv(df.ENTREZID,"hSAF_PCOS_ENTREZID.csv")

#======== GO enrichment analysis ===================

genes.DEP <- df.ENTREZID$ENTREZID

genes.DEP.FC <- row.Annotations %>% dplyr::select("GeneSymbol","Log2.FC_PCOSvsControl")
names(genes.DEP.FC)[1] <- "SYMBOL"

genes.DEP.FC1 <- plyr::join_all(list(df.ENTREZID,genes.DEP.FC),by="SYMBOL")
genes.DEP.FC2 <- as.numeric(genes.DEP.FC1$Log2.FC_PCOSvsControl)
names(genes.DEP.FC2) <- genes.DEP.FC1$ENTREZID

genes.DEP.FC2.sort <- BiocGenerics::sort(genes.DEP.FC2, decreasing=T)
genes.DEP.FC2.sort2<- as.data.frame(genes.DEP.FC2.sort)

genes.DEP.FC2.sort2$FC2<- apply(genes.DEP.FC2.sort2,1,
                                function(x){ifelse(x<0,-1,1)})

genes.DEP.FC2.sort2v<-genes.DEP.FC2.sort2$FC2
names(genes.DEP.FC2.sort2v) <- names(genes.DEP.FC2.sort)


#----Cellular component----GO classification
ggo.CC <- clusterProfiler::groupGO(gene     = genes.DEP,
                                   OrgDb    = org.Hs.eg.db,
                                   ont      = "CC",
                                   level    = 4,
                                   readable = TRUE)

head(ggo.CC)
ggo.CC.df<-as.data.frame(ggo.CC)
View(ggo.CC.df)

barplot(ggo.CC, showCategory=500)

enrichplot::goplot(ggo.CC)

#----Cellular component----GO over-representation analysis
eGO_CC <- clusterProfiler::enrichGO(gene          = genes.DEP,    #universe      = genes.DEP,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
head(eGO_CC)
eGO_CC.df <- as.data.frame(eGO_CC)

enrichplot::goplot(eGO_CC)
barplot(eGO_CC, showCategory=10) 

heatplot(eGO_CC, foldChange=genes.DEP.FC2.sort, showCategory=10) 

enrichplot::pairwise_termsim(eGO_CC,method = "JC", semData = NULL, showCategory = 10)
emapplot(eGO_CC)

eGO_CCx <- setReadable(eGO_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(eGO_CCx, foldChange = genes.DEP.FC2.sort,
         colorEdge = T,circular = T)
cnetplot(eGO_CCx, foldChange = genes.DEP.FC2.sort2v,
         colorEdge = T,circular = T)
#----Cellular component----GO Gene Set Enrichment Analysis
# here we need a numeric variable as a vector (e.g. Fold changes)

genes.DEP.FC <- row.Annotations %>% dplyr::select("GeneSymbol","Log2.FC_PCOSvsControl")
names(genes.DEP.FC)[1] <- "SYMBOL"

genes.DEP.FC1 <- plyr::join_all(list(df.ENTREZID,genes.DEP.FC),by="SYMBOL")
genes.DEP.FC2 <- as.numeric(genes.DEP.FC1$Log2.FC_PCOSvsControl)
names(genes.DEP.FC2) <- genes.DEP.FC1$SYMBOL

genes.DEP.FC2.sort <- BiocGenerics::sort(genes.DEP.FC2, decreasing=T)


cl <- parallel::makePSOCKcluster(6)                                    #  nCores=6 ??? For parallelism
doParallel::registerDoParallel(cl)

eGO_CC2 <- clusterProfiler::gseGO(geneList     = genes.DEP.FC2.sort,
                                  OrgDb        = org.Hs.eg.db,
                                  ont          = "CC",
                                  keyType = "SYMBOL",
                                  minGSSize    = 10,
                                  maxGSSize    = 500,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 1,
                                  verbose      = T,
                                  by = 'fgsea') # or  'DOSE','fgsea'
head(eGO_CC2)
eGO_CC2.df <- as.data.frame(eGO_CC2)

 

#----Mollecular Function----GO classification
ggo.MF <- clusterProfiler::groupGO(gene     = genes.DEP,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "MF",
                  level    = 3,
                  readable = TRUE)

head(ggo.MF)

ggo.MF.df <- as.data.frame(ggo.MF)

barplot(ggo.MF, showCategory=10) 

#----Mollecular Function----GO over-representation analysis
eGO_MF <- clusterProfiler::enrichGO(gene          = genes.DEP,    #universe      = genes.DEP,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1,
                   readable      = TRUE)
head(eGO_MF)
eGO_MF.df <- as.data.frame(eGO_MF)

barplot(eGO_MF, showCategory=15) 
heatplot(eGO_MF, foldChange=genes.DEP.FC2.sort, showCategory=15)

enrichplot::goplot(eGO_MF,
                   showCategory = 15,
                   color = "p.adjust",
                   layout = "sugiyama", #"kk", "sugiyama"
                   geom = "text")

eGO_MFx <- setReadable(eGO_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(eGO_MFx, foldChange = genes.DEP.FC2.sort, colorEdge = T)
cnetplot(eGO_MFx, foldChange = genes.DEP.FC2.sort2v, colorEdge = T)


#----Mollecular Function----GO Gene Set Enrichment Analysis
# here we need a numeric variable as a vector (e.g. Fold changes)

genes.DEP.FC <- row.Annotations %>% dplyr::select("GeneSymbol","Log2.FC_PCOSvsControl")
names(genes.DEP.FC)[1] <- "SYMBOL"

genes.DEP.FC1 <- plyr::join_all(list(df.ENTREZID,genes.DEP.FC),by="SYMBOL")
genes.DEP.FC2 <- as.numeric(genes.DEP.FC1$Log2.FC_PCOSvsControl)
names(genes.DEP.FC2) <- genes.DEP.FC1$SYMBOL

genes.DEP.FC2.sort <- BiocGenerics::sort(genes.DEP.FC2, decreasing=T)


eGO_MF2 <- clusterProfiler::gseGO(geneList     = genes.DEP.FC2.sort,
                                  OrgDb        = org.Hs.eg.db,
                                  ont          = "MF",
                                  keyType = "SYMBOL",
                                  minGSSize    = 10,
                                  maxGSSize    = 500,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.3,
                                  verbose      = T,
                                  by = 'fgsea') # or  'DOSE','fgsea'
head(eGO_MF2)
eGO_MF2.df <- as.data.frame(eGO_MF2)


#----Biological Process---GO classification
ggo.BP <- clusterProfiler::groupGO(gene     = genes.DEP,
                                   OrgDb    = org.Hs.eg.db,
                                   ont      = "BP",
                                   level    = 3,
                                   readable = TRUE)

head(ggo.BP)


#----Biological Process----GO over-representation analysis
eGO_BP <- clusterProfiler::enrichGO(gene          = genes.DEP,    #universe      = genes.DEP,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.1,
                                    readable      = TRUE)
head(eGO_BP)
eGO_BP.df <- as.data.frame(eGO_BP)

heatplot(eGO_BP, foldChange=genes.DEP.FC2.sort, showCategory=10)

enrichplot::goplot(eGO_BP,
                   showCategory = 10,
                   color = "p.adjust",
                   layout = "sugiyama", #"kk", "sugiyama"
                   geom = "text")

eGO_BPx <- setReadable(eGO_BP, 'org.Hs.eg.db', 'ENTREZID')

barplot(eGO_BPx, showCategory=10)
cnetplot(eGO_BPx, foldChange = genes.DEP.FC2.sort,  colorEdge = T,circular = T)
cnetplot(eGO_BPx, foldChange = genes.DEP.FC2.sort2v,  
         colorEdge = T,circular = T)

#----Biological Process----GO Gene Set Enrichment Analysis
# here we need a numeric variable as a vector (e.g. Fold changes)

genes.DEP.FC <- row.Annotations %>% dplyr::select("GeneSymbol","Log2.FC_PCOSvsControl")
names(genes.DEP.FC)[1] <- "SYMBOL"

genes.DEP.FC1 <- plyr::join_all(list(df.ENTREZID,genes.DEP.FC),by="SYMBOL")
genes.DEP.FC2 <- as.numeric(genes.DEP.FC1$Log2.FC_PCOSvsControl)
names(genes.DEP.FC2) <- genes.DEP.FC1$SYMBOL

genes.DEP.FC2.sort <- BiocGenerics::sort(genes.DEP.FC2, decreasing=T)


eGO_BP2 <- clusterProfiler::gseGO(geneList     = genes.DEP.FC2.sort,
                                  OrgDb        = org.Hs.eg.db,
                                  ont          = "BP",
                                  keyType = "SYMBOL",
                                  minGSSize    = 10,
                                  maxGSSize    = 500,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.3,
                                  verbose      = T,
                                  by = 'fgsea') # or  'DOSE','fgsea'
head(eGO_BP2)
eGO_BP2.df <- as.data.frame(eGO_BP2)



#================ KEGG enrichment analysis =============
#---KEGG pathway over-representation analysis

clusterProfiler::search_kegg_organism('hsa', by='kegg_code') # https://www.genome.jp/kegg/catalog/org_list.html 

kegg <- clusterProfiler::enrichKEGG(gene         = genes.DEP,
                                    organism     = 'hsa',
                                    pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH",  
                                    minGSSize = 10,
                                    maxGSSize = 500,
                                    qvalueCutoff = 0.1)
head(kegg)
kegg.df <- as.data.frame(kegg)

barplot(kegg, showCategory=10) 
cnetplot(kegg, foldChange = genes.DEP.FC2.sort2v,  
         colorEdge = T,circular = F)

#clusterProfiler::browseKEGG(ekegg, 'hsa01100')


#---KEGG pathway gene set enrichment analysis---

# here we need a numeric variable as a vector (e.g. Fold changes)

genes.kk.FC <- row.Annotations %>% dplyr::select("GeneSymbol","Log2.FC_PCOSvsControl")
names(genes.kk.FC)[1] <- "SYMBOL"

genes.kk.FC1 <- plyr::join_all(list(df.ENTREZID,genes.kk.FC),by="SYMBOL")
genes.kk.FC2 <- as.numeric(genes.kk.FC1$Log2.FC_PCOSvsControl)
names(genes.kk.FC2) <- genes.kk.FC1$ENTREZID

genes.kk.FC2.sort <- BiocGenerics::sort(genes.kk.FC2, decreasing=T)

ekegg <- clusterProfiler::gseKEGG(geneList     = genes.kk.FC2.sort,
                                  organism     = 'hsa',
                                  minGSSize    = 10,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.8,
                                  verbose      = FALSE)
head(ekegg)
df.ekegg <- as.data.frame(ekegg)

# hsa01100 <- pathview(gene.data  = genes.kk.FC2.sort,
#                      pathway.id = "hsa01100",
#                      species    = "hsa",
#                      limit      = list(gene=max(abs(genes.kk.FC2.sort)), cpd=1))

#---KEGG module over-representation analysis---------

mkegg <- clusterProfiler::enrichMKEGG(gene = genes.DEP,
                                      organism = 'hsa',
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)
head(mkegg)      
df.mkegg <- as.data.frame(mkegg)

dotplot(mkegg)
cnetplot(mkegg, foldChange = genes.DEP.FC2.sort2v,  
         colorEdge = T,circular = T)
#===========================================================
#              Reactome Pathways enrichment analysis 
#===========================================================
library(ReactomePA)

#--Reactome pathway over-representation analysis--
React.PA <- ReactomePA::enrichPathway(gene=df.ENTREZID$ENTREZID,
                                pAdjustMethod = "BH",   # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
                                pvalueCutoff=0.05, 
                                qvalueCutoff = 0.1,
                                readable=T)
head(as.data.frame(React.PA))

React.PA.df <- as.data.frame(React.PA)


emapplot(React.PA,line_scale = 0.5,
         layout = "kk",
         showCategory = 30,
         color = "p.adjust")

cnetplot(React.PA)
#--Reactome pathway gene set enrichment analysis--

genes.DEP.FC <- row.Annotations %>% dplyr::select("GeneSymbol","Log2.FC_PCOSvsControl")
names(genes.DEP.FC)[1] <- "SYMBOL"

genes.DEP.FC1 <- plyr::join_all(list(df.ENTREZID,genes.DEP.FC),by="SYMBOL")
genes.DEP.FC2 <- as.numeric(genes.DEP.FC1$Log2.FC_PCOSvsControl)
names(genes.DEP.FC2) <- genes.DEP.FC1$ENTREZID

genes.DEP.FC2.sort <- BiocGenerics::sort(genes.DEP.FC2, decreasing=T)


gse.Reactome <- ReactomePA::gsePathway(genes.DEP.FC2.sort,
                                       organism = "human",
                                       pvalueCutoff = 1,
                                       minGSSize = 10, 
                                       maxGSSize = 500,
                                       pAdjustMethod = "BH", 
                                       verbose = FALSE,
                                       by = "DOSE")
head(gse.Reactome)
gse.Reactome.df <- as.data.frame(gse.Reactome)

ReactomePA::viewPathway("Signal Transduction", 
                        organism = "human",
                        readable = TRUE,
                        keyType = "ENTREZID",
                        layout = "kk",
                        foldChange = genes.DEP.FC2.sort)

#======Disease enrichment analysis====================
library(DOSE)

#--Disease over-representation analysis----
DO <- DOSE::enrichDO(gene          = genes.DEP,
              ont           = "DO",
              pvalueCutoff  = 0.1,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.1,
              readable      = T)
head(DO)
DO.df <- as.data.frame(DO)

#--- Over-representation analysis for the network of cancer gene ----

ncg <- DOSE::enrichNCG(gene = genes.DEP,
                       pvalueCutoff  = 1,
                       pAdjustMethod = "BH",
                       minGSSize     = 5,
                       maxGSSize     = 500,
                       qvalueCutoff  = 1,
                       readable      = T) 
head(ncg)
ncg.df <- as.data.frame(ncg)

#------ Over-representation analysis for the disease gene network ----
 
dgn <- DOSE::enrichDGN(gene = genes.DEP,
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       minGSSize     = 5,
                       maxGSSize     = 500,
                       qvalueCutoff  = 0.1,
                       readable      = T) 
head(dgn)
dgn.df <- as.data.frame(dgn)

#--Disease gene set enrichment analysis------

gseDO <- gseDO(genes.DEP.FC2.sort,
               minGSSize     = 10,
               pvalueCutoff  = 1,
               pAdjustMethod = "BH",
               verbose       = FALSE,
               by = "fgsea")
head(gseDO, 3)
gseDO.df <- as.data.frame(gseDO)


#======= MeSH over-representation analysis====MeSH:semantic similarity analysis==============
# First, we need to load/fetch species-specific MeSH annotation database:
  
  #############################
## BioC 2.14 to BioC 3.13  ##
#############################
##
library(MeSH.Hsa.eg.db)
db <- "MeSH.Hsa.eg.db"
##---------------------------

# From BioC 3.14 (Nov. 2021, with R-4.2.0)
library(AnnotationHub)
library(MeSHDbi)
ah <- AnnotationHub(localHub=TRUE)
hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
file_hsa <- hsa[[1]]
db <- MeSHDbi::MeSHDb(file_hsa)

# In this example, we use data source from gendoo and C (Diseases) category.

library(meshes)
genes.DEP.num <- as.numeric(genes.DEP)

MeSH <- meshes::enrichMeSH(genes.DEP, 
                        MeSHDb =  "MeSH.Hsa.eg.db", 
                        database='gene2pubmed',         #gene2pubmed; gendoo: Text-mining
                       category = 'C')
MeSH.df.pubmed <- as.data.frame(MeSH) 

MeSH.df.gendoo <- as.data.frame(MeSH) 


#====== Cell Marker  ==========
cell_marker_data <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt')

## instead of `cellName`, users can use other features (e.g. `cancerType`)
cells <- cell_marker_data %>%
  dplyr::select(cellName, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', ')) %>%
  tidyr::unnest(cols = c(geneID))

#---Cell Marker over-presentaton analysis---
ecell <- clusterProfiler::enricher(genes.DEP, 
                                   pvalueCutoff = 1,
                                   qvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   TERM2GENE = cells)

ecell <- setReadable(ecell, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

head(ecell)

ecell.df <- as.data.frame(ecell)

dotplot(ecell,showCategory = 15,color = "p.adjust",
        title = "Cell Marker over-presentation analysis")
barplot(ecell)
cnetplot(ecell,foldChange = genes.DEP.FC2.sort2v, showCategory = 15, 
         colorEdge = T,circular = T)

#--Cell Marker gene set enrichment analysis
genes.DEP.FC1 <- plyr::join_all(list(df.ENTREZID,genes.DEP.FC),by="SYMBOL")
genes.DEP.FC2 <- as.numeric(genes.DEP.FC1$Log2.FC_PCOSvsControl)
names(genes.DEP.FC2) <- genes.DEP.FC1$ENTREZID

genes.DEP.FC2.sort <- BiocGenerics::sort(genes.DEP.FC2, decreasing=T)

gse.cell <- GSEA(genes.DEP.FC2.sort, 
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 TERM2GENE = cells,
                 by = "fgsea")
head(gse.cell)
gse.cell.df <- as.data.frame(gse.cell)

#=============MSigDb analysis========================
library(msigdbr)
msigdbr_species()

# select specific collection. 
# Here we use C8,	cell type signature gene sets

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame


m_t2g <- msigdbr(species = "Homo sapiens", category = "C8") %>% 
  dplyr::select(gs_name, entrez_gene)                                 # "gs_cat","gs_subcat","gs_name","gene_symbol","entrez_gene","ensembl_gene","human_gene_symbol","human_entrez_gene","human_ensembl_gene","gs_id","gs_pmid","gs_geoid","gs_exact_source","gs_url","gs_description" 
head(m_t2g)


#--MSigDb over-presentaton analysis--
em <- clusterProfiler::enricher(genes.DEP, TERM2GENE = m_t2g)

em <- setReadable(em, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

em.df <- as.data.frame(em)
#write.csv(em.df,"MSigDb_cell type signature gene sets_over-presentaton analysis.df.csv")

showCategory.v <- c("FAN_OVARY_CL3_MATURE_CUMULUS_GRANULOSA_CELL_1",
                    "FAN_OVARY_CL11_MURAL_GRANULOSA_CELL",
                    "FAN_OVARY_CL3_MATURE_CUMULUS_GRANULOSA_CELL_1",
                    "FAN_OVARY_CL10_PUTATIVE_EARLY_ATRESIA_GRANULOSA_CELL",
                    "FAN_OVARY_CL6_PUTATIVE_EARLY_ATRETIC_FOLLICLE_THECAL_CELL_2",
                    "FAN_OVARY_CL5_HEALTHY_SELECTABLE_FOLLICLE_THECAL_CELL",
                    "FAN_OVARY_CL2_PUTATIVE_EARLY_ATRETIC_FOLLICLE_THECAL_CELL_1",
                    "DESCARTES_FETAL_THYMUS_STROMAL_CELLS",
                    "FAN_OVARY_CL9_PUTATIVE_APOPTOTIC_ENDOTHELIAL_CELL",
                    "FAN_OVARY_CL7_ANGEIOGENIC_ENDOTHELIAL_CELL",
                    "FAN_OVARY_CL1_GPRC5A_TNFRS12A_HIGH_SELECTABLE_FOLLICLE_STROMAL_CELL")
dotplot(em,showCategory = showCategory.v)


#==
heatmap.clusters1<-heatmap.clusters
heatmap.clusters1$Clust.num <- apply(as.matrix(heatmap.clusters1[,"Cluster"]),1,function(x){ifelse(x=="C1",4,ifelse(x=="C2",-4,ifelse(x=="C3",-1,1)))}) 

genes.DEP.cl <- heatmap.clusters1 %>% dplyr::select("SYMBOL","Clust.num")
names(genes.DEP.cl)[1] <- "SYMBOL"

genes.DEP.cl1 <- plyr::join_all(list(df.ENTREZID,genes.DEP.cl),by="SYMBOL")
genes.DEP.cl2 <- as.numeric(genes.DEP.cl1$Clust.num)
names(genes.DEP.cl2) <- genes.DEP.cl1$ENTREZID

genes.DEP.cl2.sort <- BiocGenerics::sort(genes.DEP.cl2, decreasing=T)



cneC8<-cnetplot(em,foldChange = genes.DEP.cl2.sort,  
         colorEdge = T,circular = T,layout="kk",
         showCategory = showCategory.v,
         cex_label_category = 1,shadowtext = "none",
         cex_label_gene = 0.7)
cneC8
ggplot2::ggsave(file="MSigDb_Cells2.svg",plot=cneC8, width=10, height=8)


emapplot(em,line_scale = 0.5,layout = "kk")
#--MSigDb gene set enrichment analysis--



em.gsea <- GSEA(genes.DEP.FC2.sort, 
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            TERM2GENE = m_t2g,
            by = "fgsea")

head(em.gsea)
em.gsea.df <- as.data.frame(em.gsea)

gseadist(em.gsea, IDs,type = "density")  #type one of ’density’ or ’boxplot’

gseaplot(em.gsea, geneSetID, by = "all", title = "")

#write.csv(em1,"C8_cell type hSAF_PCOS.csv")
emapplot(em,line_scale = 0.5,layout = "kk")

#==============================================================================
#                      COMPARING DEP Heatmap clusters
#==============================================================================

cluster.m <- plyr::join_all(list(DEP_mean, heatmap.clusters), by="Accession")
cluster.m <- plyr::join_all(list(cluster.m,
                                 as.data.frame(row.Annotations[,c("Accession","Log2.FC_PCOSvsControl")])), by = "Accession")
cluster.m1 <- cluster.m %>% dplyr::select(c("SYMBOL", "Cluster", "Log2.FC_PCOSvsControl"))
rownames(cluster.m1) <- cluster.m1$SYMBOL

cl.ENTREZID = clusterProfiler::bitr(cluster.m1$SYMBOL, fromType="SYMBOL", 
                                    toType="ENTREZID", 
                                    OrgDb="org.Hs.eg.db") #  convert biological IDs to Entrez gene ID.
cluster.m2 <- plyr::join_all(list(cluster.m1,cl.ENTREZID),by="SYMBOL")

#kegg
kegg_DEP.cluster <- compareCluster(ENTREZID ~ Cluster, data=cluster.m2, 
                                   organism = "hsa",fun="enrichKEGG")
kegg_DEP.cluster <- setReadable(kegg_DEP.cluster, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

kegg_DEP.cluster.df <- as.data.frame(kegg_DEP.cluster)
#write.csv(kegg_DEP.cluster.df,"kegg_DEP.4cluster.df.csv")

p1<-dotplot(kegg_DEP.cluster,showCategory=5)+ ggtitle("KEGG")
p1
#ggplot2::ggsave(file="KEGG_4cluster.svg",plot=p1, width=10, height=10)

cnetplot(kegg_DEP.cluster)

#enrichPathway-ReactomePA
Reac_DEP.cluster <- compareCluster(ENTREZID ~ Cluster, data=cluster.m2, 
                                   organism = "human",fun="enrichPathway")
Reac_DEP.cluster <- DOSE::setReadable(Reac_DEP.cluster, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

Reac_DEP.cluster.df <- as.data.frame(Reac_DEP.cluster)
#write.csv(Reac_DEP.cluster.df,"Reac_DEP.4cluster.df.csv")



p2<-dotplot(Reac_DEP.cluster,showCategory=5)+ ggtitle("Reactome Pathways")
ggplot2::ggsave(file="Reactome_4cluster.svg",plot=p2, width=13, height=10)
cnetplot(Reac_DEP.cluster)


#GO
GO_DEP.cluster <- compareCluster(ENTREZID ~ Cluster, data=cluster.m2, fun = enrichGO, 
                                 OrgDb = org.Hs.eg.db)
GO_DEP.cluster <- setReadable(GO_DEP.cluster, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
GO_DEP.cluster.df <- as.data.frame(GO_DEP.cluster)
#write.csv(GO_DEP.cluster.df,"GO_DEP.4cluster.df.csv")

p3<- dotplot(GO_DEP.cluster,showCategory=5)+ ggtitle("GO")
#ggplot2::ggsave(file="GO_4cluster.svg",plot=p3, width=11, height=10)

cnetplot(GO_DEP.cluster)

cowplot::plot_grid(p2, p1,ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, .6))

#DO
# DO_DEP.cluster <- compareCluster(ENTREZID ~ Cluster, data=cluster.m2, fun = DOSE::enrichDO)
# DO_DEP.cluster <- setReadable(DO_DEP.cluster, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
# DO_DEP.cluster.df <- as.data.frame(DO_DEP.cluster)
# 
# dotplot(DO_DEP.cluster)
# cnetplot(DO_DEP.cluster)

#==============================================================================
#    COMPARING DEP FROM AMBEKAR 2015 (BIG FOLLICLES) AND OUR DEP FROM hSAF
#==============================================================================

genelist <- as.list(readxl::read_excel("data/DEP in Ambekar2015 and hSAF_PCOS.xlsx"))
BigFF <- clusterProfiler::bitr(genelist$`Ambekar.A 2015`, fromType="SYMBOL", 
                                   toType="ENTREZID", 
                                   OrgDb="org.Hs.eg.db")
hSAF <- clusterProfiler::bitr(genelist$hSAF_PCOS, fromType="SYMBOL", 
                              toType="ENTREZID", 
                              OrgDb="org.Hs.eg.db")
genes.gc <- list(BigFF$ENTREZID,hSAF$ENTREZID)
names(genes.gc) <- c("BigFF", "hSAF")

#KEgg
gc.kegg <- compareCluster(geneCluster = genes.gc, fun = enrichKEGG)
gc.kegg <- setReadable(gc.kegg, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

gc.kegg.df <- as.data.frame(gc.kegg)

remov <- c("Coronavirus disease - COVID-19") 

showCategory.k <- as.data.frame(gc.kegg.df[!gc.kegg.df$Description %in% remov,])
showCategory.k <- showCategory.k$Description

p5.dot <- dotplot(gc.kegg, showCategory = showCategory.k)
p5.dot
#ggplot2::ggsave(file="KEGG_Big_vs_smallFF_dot.svg",plot=p5.dot, width=10, height=10)

p5.cne <- cnetplot(gc.kegg,circular = F,showCategory=15,colorEdge = T,
                   cex_category = 1.6,
                   cex_gene = 1.6,max.overlaps =50,
                   cex_label_category = 2,
                   cex_label_gene = 1.6,
                   shadowtext = "none",
                   color_category="#E6C494")
p5.cne
ggplot2::ggsave(file="KEGG_Big_vs_smallFF_cne1.svg",plot=p5.cne, width=15, height=15)

#Reactome pathways
gc.ReactPA <- compareCluster(geneCluster = genes.gc,fun="enrichPathway",
                          organism = "human")
gc.ReactPA <- setReadable(gc.ReactPA, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

gc.ReactPA.df <- as.data.frame(gc.ReactPA)

p6.dot <- dotplot(gc.ReactPA,showCategory =10)
ggplot2::ggsave(file="ReactomePA_Big_vs_smallFF_dot.svg",plot=p6.dot, width=10, height=10)

p6.cne <- cnetplot(gc.ReactPA,circular = F,
                   cex_category = 1.6,
                   cex_gene = 1.6,
                   cex_label_category = 2,
                   cex_label_gene = 1.6,
                   shadowtext = "none",
                   color_category="#E6C494")
p6.cne
ggplot2::ggsave(file="ReactomePA_Big_vs_smallFF_cne.svg",plot=p5.cne, width=10, height=10)



# GO
gc.GO <- compareCluster(geneCluster = genes.gc, fun = enrichGO, OrgDb = org.Hs.eg.db)
gc.GO.df <- as.data.frame(gc.GO)

dotplot(gc.GO,showCategory =10)
cnetplot(gc.GO, shadowtext = "none",cex_label_category = 2)

# DO
gc.DO <- compareCluster(geneCluster = genes.gc, fun = enrichDO)
gc.DO.df <- as.data.frame(gc.DO)

dotplot(gc.DO)
cnetplot(gc.GO)

#==============================================================================
#    COMPARING DEP FROM BIG FOLLICLES papers AND OUR DEP FROM hSAF
#==============================================================================

genelist <- readxl::read_excel("data/DEP in BigF and hSAF_PCOS.xlsx")
LargeFF <- clusterProfiler::bitr(genelist$`Large Foll.`, fromType="UNIPROT", 
                               toType="ENTREZID", 
                               OrgDb="org.Hs.eg.db")
LargeFF$cond <- "LargeFF"
hSAF <- clusterProfiler::bitr(genelist$hSAF_PCOS, fromType="SYMBOL", 
                              toType="ENTREZID", 
                              OrgDb="org.Hs.eg.db")
hSAF$cond <- "hSAF"


genes.df<- rbind(LargeFF[,c("ENTREZID","cond")],hSAF[,c("ENTREZID","cond")])
genes.gc <- list(LargeFF$ENTREZID,hSAF$ENTREZID)
names(genes.gc) <- c("LargeFF", "hSAF")


#KEgg
gc.kegg <- compareCluster(geneCluster = genes.gc, fun = enrichKEGG)

gc.kegg <- compareCluster(ENTREZID ~ cond, data=genes.df,fun = enrichKEGG)

gc.kegg <- setReadable(gc.kegg, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

gc.kegg.df <- as.data.frame(gc.kegg)

remov <- c("Coronavirus disease - COVID-19") 

showCategory.k <- as.data.frame(gc.kegg.df[!gc.kegg.df$Description %in% remov,])
showCategory.k <- showCategory.k$Description

p5.dot <- dotplot(gc.kegg, showCategory = showCategory.k)
p5.dot
#ggplot2::ggsave(file="KEGG_Big_vs_smallFF_dot.svg",plot=p5.dot, width=10, height=10)

p5.cne <- cnetplot(gc.kegg,circular = F,showCategory=15,colorEdge = T,
                   cex_category = 1.6,
                   cex_gene = 1.6,max.overlaps =50,
                   cex_label_category = 2,
                   cex_label_gene = 1.6,
                   shadowtext = "none",
                   color_category="#E6C494")
p5.cne
#ggplot2::ggsave(file="KEGG_Big_vs_smallFF_cne1.svg",plot=p5.cne, width=15, height=15)

#Reactome pathways
gc.ReactPA <- compareCluster(geneCluster = genes.gc,fun="enrichPathway",pvalueCutoff=0.05,
                             organism = "human")
gc.ReactPA <- setReadable(gc.ReactPA, OrgDb = org.Hs.eg.db,keyType="ENTREZID")

gc.ReactPA.df <- as.data.frame(gc.ReactPA)

p6.dot <- dotplot(gc.ReactPA,showCategory =20)
ggplot2::ggsave(file="ReactomePA_Big_vs_smallFF_dot.svg",plot=p6.dot, width=10, height=10)

p6.cne <- cnetplot(gc.ReactPA,circular = F,
                   cex_category = 1.6,
                   cex_gene = 1.6,
                   cex_label_category = 2,
                   cex_label_gene = 1.6,
                   shadowtext = "none",
                   color_category="#E6C494")
p6.cne
ggplot2::ggsave(file="ReactomePA_Big_vs_smallFF_cne.svg",plot=p5.cne, width=10, height=10)



# GO
gc.GO <- compareCluster(geneCluster = genes.gc, fun = enrichGO, OrgDb = org.Hs.eg.db)
gc.GO.df <- as.data.frame(gc.GO)

dotplot(gc.GO,showCategory =10)
cnetplot(gc.GO, shadowtext = "none",cex_label_category = 2)

# DO
gc.DO <- compareCluster(geneCluster = genes.gc, fun = enrichDO)
gc.DO.df <- as.data.frame(gc.DO)

dotplot(gc.DO)
cnetplot(gc.GO)


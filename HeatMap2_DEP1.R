
#===Heatmap of the differentially expressed proteins=================================================

#====INSTALL PACKAGES================================================================================

# Lista de paquetes de funciones a instalar
.packages = c("BiocManager","devtools","ggplot2","ggbiplot", "tidyverse","circlize","grid",
              "GOsummaries","randomcoloR","parallel","doParallel","pheatmap","ComplexHeatmap")
#"RFunrichWebService",
# Instala los paquetes sin? los tienes instalados
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  install_github("vqv/ggbiplot")
  install_github("fxia22/RadaFDR")
  BiocManager::install(c("mixOmics","phyloseq","GOsummaries"))
}

# Carga los paquetes sin? los tienes cargados
lapply(.packages, require, character.only=TRUE)

#=====Loading dataset===================================================================== 

DEP_proteins <- as.data.frame(readxl::read_excel("data/20200229_FF_Dep_small_R2_abundances_Statistic_PCOSipp1_normalized in R1_zscore.xlsx")) # Open the data
rownames(DEP_proteins) <- paste(DEP_proteins$Accession,".",DEP_proteins$`Gene`)
class(DEP_proteins)


DEP_mean <- as.data.frame(readxl::read_excel("data/mean groups from perseu2_a.xlsx")) # Open the data
rownames(DEP_mean) <- paste(DEP_mean$Accession,".",DEP_mean$`Gene`)
class(DEP_mean)

#main_table <- Big.table %>% dplyr::select(contains(c("_1","_2"))) # Select columns whose names contains "Control" and "PCOS" (Intenstity columns)
main_DEP <- DEP_proteins[, 4:ncol(DEP_proteins)]
main_DEP[1:5,1:4]

Annotations <- as.data.frame(readxl::read_excel("data/Copy of PCOSprojectDataforLund (3).xlsx")) # Open the data
row.names(Annotations) <- Annotations$ID
head(Annotations)

row.Annotations <- DEP_proteins %>% dplyr::select(contains(c("Analysis")))
row.Annotations$Protein.ID <- rownames(row.Annotations)

row.Annotations <- as.data.frame(readxl::read_excel("data/row annotation_mean groups from perseu2a.xlsx")) # Open the data
row.names(row.Annotations) <- paste(row.Annotations$protein,".",row.Annotations$GeneSymbol)
head(row.Annotations)
#========================================================================================

##==========HEATMAP================================
library(ComplexHeatmap)
library(pheatmap)

cim.matrix1.m <- as.matrix(main_DEP)
rownames(cim.matrix1.m) <- rownames(main_DEP)

heatmap(cim.matrix1.m, scale = "none")

# col<- colorRampPalette(c("blue", "white", "red"))(256)
# library("RColorBrewer")
# col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
# 
# # Use RColorBrewer color palette names
# library("RColorBrewer")
# col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
# heatmap(cim.matrix1.m, scale = "none", col =  col, 
#         ColSideColors = c(rep("purple", 6), rep("orange", 7)))
# 
# library("gplots")
# heatmap.2(cim.matrix1.m, scale = "none", col = bluered(100), 
#           trace = "none", density.info = "none")
# 
# library("d3heatmap")
# d3heatmap(cim.matrix1.m, colors = "RdYlBu",
#           k_row = 4, # Number of groups in rows
#           k_col = 2 # Number of groups in columns
# )
# 
# library(ComplexHeatmap)
# Heatmap(cim.matrix1.m, 
#         name = "intensity", #title of legend
#         column_title = "patients", row_title = "proteins",
#         row_names_gp = gpar(fontsize = 7) # Text size for row names
# )
# 
# library(circlize)
# mycols <- colorRamp2(breaks = c(-2, 0, 2), 
#                      colors = c("green", "white", "red"))
# Heatmap(cim.matrix1.m, name = "sPLS-DA", col = mycols)
# 
# library("circlize")
# library("RColorBrewer")
# Heatmap(cim.matrix1.m, name = "mtcars",
#         col = colorRamp2(c(-2, 0, 2), brewer.pal(n=3, name="RdBu")))
# 
# library(dendextend)
# row_dend = hclust(dist(cim.matrix1.m)) # row clustering
# col_dend = hclust(dist(t(cim.matrix1.m))) # column clustering
# Heatmap(cim.matrix1.m, name = "mtcars", 
#         row_names_gp = gpar(fontsize = 6.5),
#         cluster_rows = color_branches(row_dend, k = 4),
#         cluster_columns = color_branches(col_dend, k = 2))
library(tidyverse)

cim.matrix1.m1 <- as.data.frame(cim.matrix1.m)
cim.matrix1.m1$Prot.ID <- rownames(cim.matrix1.m)

row.Annotations1 <- row.Annotations %>% dplyr::select(c("Analysis","Secreted","LOG10p.value",
                                                        "Granul.Cell","more concentrated in hSAF","Occyte", 
                                                        "occyte competence(Pla.I)","Ambekar A. 2015","Zhang X 2019"))
row.Annotations1$Prot.ID <- rownames(row.Annotations1)
  
  
cim.matrix1.m1<- plyr::join_all(list(cim.matrix1.m1, row.Annotations1),by="Prot.ID")
row.names(cim.matrix1.m1)<- cim.matrix1.m1$Prot.ID
cim.matrix1.m1 <- cim.matrix1.m1[,!(names(cim.matrix1.m1) %in% "Prot.ID")]

ff <- as.data.frame(t(cim.matrix1.m))
ff$sample <- rownames(ff)
# ff <- plyr::join_all(list(ff,col.Annotation),by = "sample")
# col.Annotation1 <- as.data.frame(ff[,c("sample","condition","Follicular size","patient")])

matrix1 <- as.matrix(cim.matrix1.m1[,1:20])
matrix2 <- as.data.frame(cim.matrix1.m1[,21:ncol(cim.matrix1.m1)])




# TOP ANNOTATIONS. 
# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
library(circlize)


col.top = list(condition = c("non-PCO" = "sienna1","PCO" = "steelblue3"))

# Create the heatmap annotation
ha.top <- HeatmapAnnotation(
  condition = Annotations$Cond,
  col = col.top, show_legend = T,simple_anno_size = unit(0.3, "cm"))

# row.annot <- rowAnnotation(`-LOG10p-value` = anno_lines(matrix2$LOG10p.value,
#                                                         gp = gpar(col = 2:3), 
#                                                         add_points = TRUE, 
#                                                         pt_gp = gpar(col = 5:6), 
#                                                         pch = c(1, 16)),
#                            prot = anno_text(matrix2$Gene, gp = gpar(fontsize = 5)))

row.annot <- rowAnnotation(`-LOG10p-value` = anno_lines(cbind(matrix2$LOG10p.value,rep(1.3,nrow(matrix2))),
                                                        gp = gpar(col = c("blue","red")), 
                                                        pt_gp = gpar(col = "blue"), add_points = TRUE,
                                                        size = unit(2, "mm"),
                                                        pch = c('o','|'),width = unit(1.8, "cm")))


col.ann_rows = list(Secreted = c("+" = "green"),
                   `Granul.Cell`= c("+" = "cyan3"),
                    Occyte = c("+" = "magenta"),
                   `occyte competence(Pla.I)`= c("+" = "brown1","-"= "deepskyblue2"),
                   `Ambekar A. 2015`= c("+" = "brown1","-"= "deepskyblue2"),
                   `Zhang X 2019`= c("+" = "brown1","-"= "deepskyblue2"))



row.annot1 <- ComplexHeatmap::rowAnnotation(Secreted = matrix2$Secreted,
                                            `Granul.Cell`= matrix2$Granul.Cell,
                                             Occyte = matrix2$Occyte,
                                            `occyte competence(Pla.I)`= matrix2$`occyte competence(Pla.I)`,
                                            `Ambekar A. 2015` = matrix2$`Ambekar A. 2015`,
                                            `Zhang X 2019` = matrix2$`Zhang X 2019`,
                                            show_legend = T,simple_anno_size = ggplot2::unit(0.3, "cm"),
                                            col=col.ann_rows,na_col = "white",
                                            border = T, 
                                            gap = unit(2, "points"),
                                            show_annotation_name = F,
                                            gp = gpar(col = "white"))


col3<- colorRampPalette(c("mediumblue","dodgerblue3", "white","coral2","red3"))(100)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# Combine the heatmap and the annotation

heatm <- ComplexHeatmap::Heatmap(matrix1,top_annotation = ha.top,
                        name = "s",
                        row_names_gp = grid::gpar(fontsize = 5),
                        column_names_gp = grid::gpar(fontsize = 8),
                        show_row_names = T,show_column_names = T,
                        col = col3,
                        left_annotation = row.annot)
draw(heatm)

pa = cluster::pam(matrix1, k = 2)

heatm <- ComplexHeatmap::Heatmap(matrix1,
                                 name = "s",
                                 column_split = c(rep(".",10),rep("",10)),
                                 row_split = paste0("C", pa$clustering),
                                 row_names_gp = grid::gpar(fontsize = 5),
                                 column_names_gp = grid::gpar(fontsize = 8),
                                 show_row_names = F,show_column_names = T,
                                 col = col_fun,na_col = "grey",
                                 right_annotation = row.annot1,
                                 top_annotation = ha.top,
                                 clustering_distance_rows = "euclidean",
                                 clustering_method_rows = "complete",
                                 clustering_distance_columns = "euclidean",
                                 clustering_method_columns = "complete",
                                 row_dend_side = c("right"),
                                 row_dend_width = unit(8, "mm"),
                                 column_dend_side = c("top"),
                                 column_dend_height = unit(8, "mm"),
                                 show_column_dend = TRUE,
                                 show_row_dend = F,
                                 border = TRUE,
                                 border_gp = gpar(col = "grey", lty = 1),
                                 row_title = "Differentially expressed proteins",
                                 column_title_side = "bottom",
                                 row_title_side = "left",
                                 column_title = "non-PCO   PCO")
draw(heatm)

# SVG graphics device
svg("DEP_heatmap.svg")

# Code of the plot
plot(heatm)

# Close the graphics device
dev.off() 

pheatmap(matrix1,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "average")


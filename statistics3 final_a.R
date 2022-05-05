##=====Follicular fluid=====PCOS project
##http://mixomics.org/graphics/sample-plots/plotindiv/

#setwd("G:/PhD/work/02 Follicle study/PCOS samples/sPLS-DA/PCOS_Project/PCOS_data_analysis/data")

#====INSTALL PACKAGES================================================================================

# Lista de paquetes de funciones a instalar
.packages = c("BiocManager","devtools","ggplot2","ggbiplot", "pca3d","pcadapt","outliers","igraph",
              "rgl","graphics","reshape2","dplyr","ggpubr","remotes",
              "FactoMineR", "factoextra","corrplot","ggpubr","fpc", "NbClust","mixOmics", 
              "phyloseq","lme4","nlme","car","plotly","RadaFDR","GOplot","tidyverse",
              "GOsummaries","randomcoloR","parallel","doParallel","afex","ggbeeswarm",
              "emmeans","psych", "pheatmap","lmerTest","nlme","Seurat","proBatch")
#"RFunrichWebService",
# Instala los paquetes sin? los tienes instalados
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  install_github("vqv/ggbiplot")
  install_github("fxia22/RadaFDR")
  BiocManager::install(c("mixOmics","phyloseq","GOsummaries","Seurat"))
  devtools::install_github('symbioticMe/proBatch', build_vignettes = TRUE)
}

# Carga los paquetes sin? los tienes cargados
lapply(.packages, require, character.only=TRUE)

#=========================================================================================

# Big_Tabla <- read.delim("data/20200229_FF_Dep_small_R2_abundances_Statistic_PCOSipp1.txt", header = TRUE,row.names = NULL)
# 
# head(Big_Tabla)
# row.names(Big_Tabla)<- paste(Big_Tabla$Accession,".",Big_Tabla$Gene.name)
# 
# Big_Tabla1 <- Big_Tabla[,14:ncol(Big_Tabla)] 
# head(Big_Tabla1)
# 
# annotation_colum <- read.delim("data/annotation.txt", header = TRUE,row.names = 1)
# head(annotation_colum)
# annotation_colum$sample <- row.names(annotation_colum)


#=====Loading dataset===================================================================== 

Big.table <- as.data.frame(readxl::read_excel("data/20210628_FF_Quantified proteins.xlsx")) # Open the data
rownames(Big.table) <- paste(Big.table$Accession,".",Big.table$`Gene name`)
class(Big.table)

#main_table <- Big.table %>% dplyr::select(contains(c("_1","_2"))) # Select columns whose names contains "Control" and "PCOS" (Intenstity columns)
main_table <- Big.table[, 6:ncol(Big.table)]
main_table[1:5,1:5]

Annotations <- as.data.frame(readxl::read_excel("data/Copy of PCOSprojectDataforLund (3).xlsx")) # Open the data
row.names(Annotations) <- Annotations$ID
head(Annotations)
#=========================================================================================

#=========================FUNCTIONS======================================================

# Loading from the current directory (run getwd())the Rfile that contains all the functions created for the analyses

source('PCO_RFunctions.R')

#====================Data PRE-PROCESSING=================================

# checking the data
library(proBatch)
library(ggpubr)

color.list <- sample_annotation_to_colors(Annotations,
                                          factor_columns = c('Cond',"group_Size2", 'PCA batch'),
                                          numeric_columns = c('Sample no'))

pca1 = plot_PCA(main_table, Annotations, color_by = 'Cond',sample_id_col="ID",
                plot_title = 'Condition', color_scheme = color.list[['Cond']])
pca2 = plot_PCA(main_table, Annotations, color_by = 'group_Size2',sample_id_col="ID",
                plot_title = 'group_Size2', color_scheme = color.list[['group_Size2']])
pca3 = plot_PCA(main_table, Annotations, color_by = 'PCA batch',sample_id_col="ID",
                plot_title = 'PCA batch', color_scheme = color.list[['PCA batch']])

ggpubr::ggarrange(pca1, pca2, pca3, ncol = 3, nrow = 1)

boxplot(main_table,las = 2) 




## (1) Log2 transformation

main_table.Log2 <- log2(main_table)
boxplot(main_table.Log2,las = 2)      # Checking normalization of the data


pca1 = plot_PCA(main_table.Log2, Annotations, color_by = 'Cond',sample_id_col="ID",
                plot_title = 'Condition', color_scheme = color.list[['Cond']])
pca2 = plot_PCA(main_table.Log2, Annotations, color_by = 'group_Size2',sample_id_col="ID",
                plot_title = 'group_Size2', color_scheme = color.list[['group_Size2']])
pca3 = plot_PCA(main_table.Log2, Annotations, color_by = 'PCA batch',sample_id_col="ID",
                plot_title = 'PCA batch', color_scheme = color.list[['PCA batch']])


ggpubr::ggarrange(pca1, pca2, pca3, ncol = 3, nrow = 1)


main_table.Log2.topca <- as.data.frame(t(na.omit(main_table.Log2)))

results.pca <- mixOmics::pca(main_table.Log2.topca,ncomp = 3, scale = T)

plotIndiv(results.pca,ind.names = TRUE, cex=1.5)





#whitout missing values
pca1 = plot_PCA(na.omit(main_table.Log2), Annotations, color_by = 'Cond',sample_id_col="ID",
                plot_title = 'Condition', color_scheme = color.list[['Cond']])
pca2 = plot_PCA(na.omit(main_table.Log2), Annotations, color_by = 'group_Size2',sample_id_col="ID",
                plot_title = 'group_Size2', color_scheme = color.list[['group_Size2']])plot_PCA(data_matrix, sample_annotation,
  feature_id_col = "peptide_group_label",
  sample_id_col = "FullRunName", color_by = "MS_batch",
  PC_to_plot = c(1, 2), fill_the_missing = -1,
  color_scheme = "brewer", filename = NULL, width = NA,
  height = NA, units = c("cm", "in", "mm"), plot_title = NULL,
  theme = "classic", base_size = 20)
pca3 = plot_PCA(na.omit(main_table.Log2), Annotations, color_by = 'PCA batch',sample_id_col="ID",
                plot_title = 'PCA batch', color_scheme = color.list[['PCA batch']])

ggpubr::ggarrange(pca1, pca2, pca3, ncol = 3, nrow = 1)


## (2) Apply filtering (At least 1 Valid value in Total)
main_table.Log2.F = filter_valids(main_table.Log2,
                                  conditions = c("C", "P"),
                                  min_count = c(3, 3),
                                  at_least_one = T)

main_table.Log2.F1 = filter(main_table.Log2.F, KEEP)   ## Remove rows where KEEP is FALSE

main_table.Log2.F1 <- dplyr::select(main_table.Log2.F1, -c(KEEP, C,P))

## (3) standardization (Substract median)

main_table.Log2.F.N = median_centering(main_table.Log2.F1)
boxplot(main_table.Log2.F.N,las = 2)      # Checking normalization of the data

pca1 = plot_PCA(na.omit(main_table.Log2.F.N), Annotations, color_by = 'Cond',sample_id_col="ID",
                plot_title = 'Condition', color_scheme = color.list[['Cond']])
pca2 = plot_PCA(na.omit(main_table.Log2.F.N), Annotations, color_by = 'group_Size2',sample_id_col="ID",
                plot_title = 'group_Size2', color_scheme = color.list[['group_Size2']])
pca3 = plot_PCA(na.omit(main_table.Log2.F.N), Annotations, color_by = 'PCA batch',sample_id_col="ID",
                plot_title = 'PCA batch', color_scheme = color.list[['PCA batch']])

ggpubr::ggarrange(pca1, pca2, pca3, ncol = 3, nrow = 1)


results.pca <- mixOmics::pca(t(main_table.Log2.F.N),ncomp = 3, scale = T)

plotIndiv(results.pca,ind.names = TRUE, cex=1.5)


# Coeficient of variation among samples

cv.prot <- apply(main_table.Log2.F.N,1,function(x){sd(x)/mean(x)*100})
cv.prot

main_table.Log2.F.N_cv <- main_table.Log2.F.N
main_table.Log2.F.N_cv$CV <- cv.prot
main_table.Log2.F.N_cv50 <- subset(main_table.Log2.F.N_cv,main_table.Log2.F.N_cv$CV>50)
main_table.Log2.F.N_cv50 <- main_table.Log2.F.N_cv50 %>% dplyr::select(-CV)
main_table.Log2.F.N_cv50.pca <- as.data.frame(t(main_table.Log2.F.N_cv50))
main_table.Log2.F.N_cv50.pca$ID <- rownames(main_table.Log2.F.N_cv50.pca)

#factorMine object
short.Annotations <- Annotations[,c("ID","Diameter","Cond","group_Size")]

facMine <- plyr::join_all(list(main_table.Log2.F.N_cv50.pca,
                               short.Annotations),by="ID")
rownames(facMine)<- facMine$ID
facMine1<- facMine %>% dplyr::select(-ID)

#facMine.active <- facMine[,1:nrow(main_table.Log2.F.N_cv50)]


res.pca <- FactoMineR::PCA(facMine1, 
                           quanti.sup = nrow(main_table.Log2.F.N_cv50)+1,
                           quali.sup = (nrow(main_table.Log2.F.N_cv50)+2):(nrow(main_table.Log2.F.N_cv50)+3),
                           graph = F)

fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) # Avoid text overlapping (slow if many points))


fviz_pca_biplot(res.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = facMine$Diameter,
                fill.ind = facMine$Cond,
                col.ind = "black",
                addEllipses = F,
                invisible = "var",
                col.quanti.sup = "red",
                col.quali.sup = "green",
                # Color variable by groups
                #col.var = factor(c("sepal", "sepal", "petal", "petal")),
                
                legend.title = list(fill = "Cond", pointsize = "Diameter"),
                repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors



#write.csv(main_table.Log2.F.N_cv50,"proteins with more than 50 perc of CV.csv")

## (4) Filtering missing values (at least 7 valid values in at least one condition)

df.log2.N = filter_valids(main_table.Log2.F.N,
                          conditions = c("C", "P"),
                          min_count = c(7, 7),
                          at_least_one = T)

df.log2.N1 = filter(df.log2.N, KEEP)   ## Remove rows where KEEP is FALSE

df.log2.N1 <- filter(df.log2.N1, C > 0) # removing proteins with 0 values in C
df.log2.N1 <- filter(df.log2.N1, P > 0) # removing proteins with 0 values in P


df.log2.N1 <- dplyr::select(df.log2.N1, -c(KEEP, C,P))
boxplot(df.log2.N1,las = 2)      # Checking normalization of the data

pca1 = plot_PCA(df.log2.N1, Annotations, color_by = 'Cond',sample_id_col="ID",
                plot_title = 'Condition', color_scheme = color.list[['Cond']])
pca2 = plot_PCA(df.log2.N1, Annotations, color_by = 'group_Size2',sample_id_col="ID",
                plot_title = 'group_Size2', color_scheme = color.list[['group_Size2']])
pca3 = plot_PCA(df.log2.N1, Annotations, color_by = 'PCA batch',sample_id_col="ID",
                plot_title = 'PCA batch', color_scheme = color.list[['PCA batch']])

ggpubr::ggarrange(pca1, pca2, pca3, ncol = 3, nrow = 1)
#write.csv(df.log2.N1, "data/20200229_FF_Dep_small_R2_abundances_Statistic_PCOSipp1_normalized in R.csv")

# (4.1) ploting plots to check normality
gg<- melt(df.log2.N1)

p1 <- ggplot(gg, aes(x=variable, y=value)) +
  geom_boxplot()+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title="Normalization: Log2 + sample median substraction",
       x = "samples", y= "Intensity (Log2 scaled)")

p2 <- ggplot(gg, aes(x=value, color=variable)) +
  geom_density()+theme_bw()+theme(legend.position="none")+
  scale_x_continuous(limits = c(-20, 20))+
  labs(title="Normalization: Log2 + sample median substraction",
       x = "Intensity (Log2 scaled)")

multiplot(p1, p2, cols=2)


## (5) selecting on-off proteins with more than 3 values in at least one condition

on_offs = filter_valids(main_table.Log2,
                        conditions = c("C", "P"),
                        min_count = c(3, 3),
                        at_least_one = T)
on_offs1 <- filter(on_offs, KEEP)
on_offs1_C <-  filter(on_offs1, C==0)
on_offs1_P <-  filter(on_offs1, P==0)

on_offs2 <- rbind(on_offs1_C,on_offs1_P)
on_offs2 <- dplyr::select(on_offs2, -KEEP)
on_offs3 <- dplyr::select(on_offs2, c(C,P))
dim(on_offs3)
head(on_offs3)
##===================================================================================

#========================DATA ANALAYSIS==============================================

# (1) UNIVARIATE ANALYSIS ------------------non-PCOS vs PCOS

# (1.1) Checking differences in Folliculle Diameter btw both non-PCOS vs PCOS 

ggplot(Annotations, aes(x = Diameter, fill = Cond)) + geom_density(alpha = 0.5) + theme_bw()

ggplot(Annotations, aes(y = Diameter,x=Cond, fill = Cond)) + 
  geom_boxplot(alpha = 0.5)

# coefficient of variation (Follicular Diameter)
CV.foll <- (sd(Annotations$Diameter)/mean(Annotations$Diameter))*100 
CV.foll


# (2)----------Linear Regression-------------non-PCOS vs PCOS

df.regres <- df.log2.N1

groups.ttest <- Annotations

GR <- as.factor(groups.ttest$Cond)
pat.GR<-as.factor(groups.ttest$Patient) 

####    Function DEP_prot      ####
# m = matrix with samples as columns and proteins as rows
# annotation.table = Annotation matrix
# groups = name of the column that contain the groups that should be compared (factor)  
# cofounders =  vector of characteres with the names of the columns that contain the cofounders variables (ex: c("Diameter", "sex","age").

# m=df.regres
# annotation.table = groups.ttest
# groups = "Cond"
# id_ind = "Patient"
# cofounders <- c("Diameter")
# i=1

DEP_prot <- function(m, annotation.table, groups, graphs=FALSE,cofounders=NULL, id_ind=NULL) {
  library(dplyr)
  library(tidyverse)
  
  groups1 <- as.factor(annotation.table[,groups])
  
  #--- making tables where the results will be stored---
  
  residuals.LM <- matrix(nrow=nrow(m),ncol = ncol(m), dimnames = list(row.names(m),paste("R_",colnames(m), sep = "")))  ######
  
  ttest <- matrix(nrow=nrow(m),ncol = 6, dimnames = list(rownames(m),c("p.val_ttest","Low.L", "Upp.L",
                                                                       paste("mean.",levels(groups1)[1],sep=""),
                                                                       paste("mean.",levels(groups1)[2],sep=""),
                                                                       "Log2.FC")))
  LR.results1 <- matrix(nrow=nrow(m),ncol = 3, dimnames = list(rownames(m),c("Beta","R2", paste("p.val_LR.",groups,sep = ""))))
  LR.results2 <- matrix(nrow=nrow(m),ncol = 5, dimnames = list(rownames(m),c(paste("B.",groups,sep = ""),
                                                                             paste("B.",cofounders,sep = ""),
                                                                             "R2",
                                                                             paste("p-val.",groups,sep = ""),
                                                                             paste("p-val.",cofounders,sep = ""))))
  LR.Follsize <- matrix(nrow = nrow(m), ncol = 4, dimnames = list(rownames(m), c("p.val_small.Cond", "Beta_small.Cond","p.val_med.Cond", "Beta_med.Cond")))
  correlation.Size <- matrix(nrow = nrow(m), ncol = 4, dimnames = list(rownames(m), c("p.val_PCOS", "r_PCOS","p.val_Ctrl", "r_Ctrl")))
  
  #-----
  
  m1 = t(m)
  groups2 <- annotation.table %>% dplyr::select(all_of(c(groups,id_ind,cofounders)))
  groups2$sample <- rownames(groups2)
  
  for (i in 1:ncol(m1)) {
    
    prot.data <- as.data.frame(m1[,i])
    
    prot.data$sample <- rownames(prot.data)
    
    data <- plyr::join_all(list(prot.data,groups2),by="sample")
    row.names(data) <- data$sample
    colnames(data)[1] <- "protein"
    data <- data[,-2]
    data <- na.omit(data)
    
    
    #====Student t-test=========
    
    ttest1 <- t.test(protein ~ Cond, data = data, 
                     paired = FALSE, var.equal = FALSE,
                     alternative = "two.sided")
    
    ttest[i,1] <- ttest1$p.value                         # p-value
    ttest[i,2] <- ttest1$conf.int[1]                     # DL
    ttest[i,3] <- ttest1$conf.int[2]                     # UppL
    ttest[i,4] <- ttest1$estimate[1]                     # mean cond 1
    ttest[i,5] <- ttest1$estimate[2]                     # mean cond 2
    ttest[i,6] <- ttest1$estimate[2]-ttest1$estimate[1]  # Fold change
    
    #===Linear regression=======
    
    # (1) working with residuals
    
    #LR01 <- lmerTest::lmer(protein ~  Diameter + (1|Patient), data = data)
    
    LR01 <- lm(protein ~ Diameter, data = data)
    summary(LR01)
    resid <- residuals(LR01)
    
    data$residuals <- resid
    

    LR <- lm(residuals ~ Cond, data = data)
    summary(LR)
    
    LR.results1[i,1] <- LR$coefficients[2]
    LR.results1[i,2] <- summary(LR)$r.squared
    LR.results1[i,3] <- summary(LR)$coefficients[2,4]
    
    # (2) working cofounders
    
    LR1 <- lm(protein ~  Cond + Diameter, data = data)
    
    #LR1 <- lmerTest::lmer(protein ~  Cond + Diameter + Diameter*Cond + (1|Patient), data = data)
    
    summary(LR1)
    
    LR.results2[i,1] <- LR1$coefficients[2]               # beta Cond
    LR.results2[i,2] <- LR1$coefficients[3]               # beta Diameter

    LR.results2[i,3] <- summary(LR1)$r.squared            # R2
    LR.results2[i,4] <- summary(LR1)$coefficients[2,4]    # p.val Cond
    LR.results2[i,5] <- summary(LR1)$coefficients[3,4]    # p.val Diameter
    
    
   
    # Stratifing by Cond
    
    # data_PCOS <- subset(data,data$Cond=="PCOS")
    # data_C <- subset(data,data$Cond=="non-PCOS")
    # 
    # i
    # correlation.PCOS <- cor.test(data_PCOS$protein,data_PCOS$Diameter,data=data_PCOS)
    # 
    # correlation.Size[i,1] <- correlation.PCOS$p.value
    # correlation.Size[i,2] <- correlation.PCOS$estimate
    # 
    # correlation.Ctrl <- cor.test(data_C$protein,data_C$Diameter,data=data_C)
    # correlation.Size[i,3] <- correlation.Ctrl$p.value
    # correlation.Size[i,4] <- correlation.Ctrl$estimate

    #tuk.sig <- cut(pv[i,], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("p<.001", "p<.01", "p<.05","N/S"))
  }
  
  
  #__________plots
  #library(ggpmisc)
  
  if(graphs == TRUE){
    dir.create("Plots.Reg")
    
    for (i in 1:ncol(m1)){
      png(filename = paste("Plots.Reg/","correlation_","plot_",rownames(m)[i],".png",sep=""),width = 500, height = 480, units = "px", pointsize = 16,
          bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
      
      p <- ggplot(data=data, mapping=aes(x=Diameter, y=protein, colour=Cond)) +
        geom_point() + theme_bw()+
        geom_smooth(method='lm', se=TRUE) +
        scale_colour_discrete(guide=guide_legend(title.position='left', title.hjust=1))+
        #geom_hline(yintercept=0, colour='gray35', linetype='dashed')+
        labs(subtitle=rownames(m)[i],
             caption = paste("PCOS_p.val =  ", round(correlation.PCOS$p.value,3), "   ",
                             "PCOS_r = ",round(correlation.PCOS$estimate,2),
                             " ;   Ctrl_p.val = ",round(correlation.Ctrl$p.value,3), "   ",
                             "Ctrl_r = ",round(correlation.Ctrl$estimate,2),sep=""),
             x = "Follicular Diameter (mm)", y = "Intensity (Log2 scaled)")
      print(p)
      dev.off()
      }
  }
  ttest <- as.data.frame(ttest)
  ttest$q.value_ttest <- p.adjust(ttest$p.val_ttest, method = "fdr")
  #ttest$q.v.sig <- cut(ttest$q.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("p<.001", "p<.01", "p<.05","N/S"))
  
  LR.results1 <- as.data.frame(LR.results1)
  LR.results1$q.value_LR1 <- p.adjust(LR.results1$p.val_LR.Cond, method = "fdr")
  
  LR.results2 <- as.data.frame(LR.results2)
  LR.results2$q.value_LR2.Cond <- p.adjust(LR.results2$`p-val.Cond`, method = "fdr")
  LR.results2$q.value_LR2.Diam <- p.adjust(LR.results2$`p-val.Diameter`, method = "fdr")
    
  results <- as.data.frame(cbind(ttest,LR.results1,LR.results2))#,correlation.Size))
 
  list.results <- list(results,residuals.LM)
  return(results)
}


DEP <- as.data.frame(DEP_prot(m=df.regres,
                              cofounders=c("Diameter"),
                              annotation.table=groups.ttest,
                              groups="Cond",
                              id_ind="Patient",
                              graphs = FALSE))
DEP$protein <- rownames(DEP)
######################################################################################

# (2) MULTIVARIATE ANALYSIS ---(mixOmics package)---non-PCOS vs PCOS------
 
# This analysis is not multilevel

library(mixOmics)

X <- as.data.frame(t(df.log2.N1))
Y <- as.factor(Annotations$Cond)

# Patient indicates the repeated measurements
# setup the design matrix by indicating the repeated measurements
design <- data.frame(Patient = Annotations$Patient)


# (1.1)----PCA-------------------------------------
## Apply filtering to keep proteins with 100% of valid values
    # X.pca0 <- as.data.frame(t(X))
    # X.pca = filter_valids(X.pca0,
    #                       conditions = c("C", "P"),
    #                       min_count = c(10, 10),
    #                       at_least_one = F)
    # X.pca <- filter(X.pca, KEEP)
    # X.pca <- dplyr::select(X.pca, -c(KEEP,C,P))
    # X.pca <- as.data.frame(t(X.pca))

X.pca <- X
# pca
pca <- mixOmics::pca(X.pca, ncomp = 3, scale = F, center = T)
pca

#2D plot per condition (non-PCOS vs PCOS)
pca1 <- mixOmics::plotIndiv(pca, ind.names = F, legend = T,comp=c(1,2),
                    ellipse = T,cex = 1, style="ggplot2",
                    group = Y,
                    title = 'PCA')

#3D plot
# mixOmics::plotIndiv(pca.multilevel, ind.names = F, legend = T,comp=c(1,2,3),
#                     ellipse = F,cex = 1, style="3d",
#                     group = Annotations$group_Size,
#                     title = 'PCA multilevel')

# (1.2) ----PLS-DA multilevel model----------------
## Apply filtering to keep proteins with 100% of valid values
    # X.PLS0 <- as.data.frame(t(X))
    # X.PLS = filter_valids(X.PLS0,
    #                       conditions = c("C", "P"),
    #                       min_count = c(10, 10),
    #                       at_least_one = F)
    # X.PLS <- filter(X.PLS, KEEP)
    # X.PLS <- dplyr::select(X.PLS, -c(KEEP,C,P))
    # X.PLS <- as.data.frame(t(X.PLS))

X.PLS <- X

# pls
plsda <- plsda(X.PLS, Y = Y, ncomp = 3)

#2D plot per condition (non-PCOS vs PCOS)
pls1 <- plotIndiv(plsda, ind.names = F, 
                  ellipse =T, ellipse.level=.95,
                  star = F, legend=T,group = Y,
                  style = 'ggplot2',cex=1,title = "PLS-DA")

# 3D plot
# plotIndiv(TP.plsda.multilevel, ind.names = F, ellipse =T,ellipse.level=.95,star = F, legend=T,
#           style = '3d',cex=0.5,title = "PLS-DA")

# Heatmap

col.ID <- randomcoloR::distinctColorPalette(k=10)[as.factor(Annotations$Patient)]
col.groups <- color.mixo(1:2)[as.factor(Annotations$Cond)]

cim(plsda,comp=c(1:2),
    margins = c(4,10),zoom=F,keysize = c(1,1),transpose = T,
    row.sideColors = cbind(col.groups, col.ID),
    clust.method=c("complete","complete"),
    row.cex = 0.8,col.cex = 0.8,
    row.names = Annotations$ID,
    col.names = FALSE, legend=list(legend = c(levels(as.factor(Annotations$Cond))), 
                                   col = color.mixo(1:3),
                                   title = "Condition", cex = 0.8))

# The plot Loading function displays the loading weights, 
plotLoadings(plsda, comp = 1, title = 'Loadings on comp 1', size.name = 0.6,
             contrib = 'max', method = 'mean',ndisplay=50)#,ndisplay=50)



# (1.3) ------sPLS-DA model -------before Tuning----------------------------
## Apply filtering to keep proteins with 100% of valid values

    # X.sPLS0 <- as.data.frame(t(X))
    # X.sPLS = filter_valids(X.sPLS0,
    #                       conditions = c("C", "P"),
    #                       min_count = c(10, 10),
    #                       at_least_one = F)
    # X.sPLS <- filter(X.sPLS, KEEP)
    # X.sPLS <- dplyr::select(X.sPLS, -c(KEEP,C,P))
    # boxplot(X.sPLS)
    # X.sPLS <- as.data.frame(t(X.sPLS))

X.sPLS <- X

# sPLS-DA
splsda <- splsda(X.sPLS, Y = Y,ncomp = 5, scale = T)

#----Tuning parameters and numerical outputs___sPLS-DA

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat

## 1 - The number of components to retain

cl <- makePSOCKcluster(6)                                    #  nCores=6 ??? For parallelism
registerDoParallel(cl)
MyPerf.plsda1 <- perf(splsda, validation = "Mfold", folds = 5, dist="all",         # folds :: https://machinelearningmastery.com/k-fold-cross-validation/#:~:text=Cross%2Dvalidation%20is%20a%20resampling,k%2Dfold%20cross%2Dvalidation. 
                      progressBar = T, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda1, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda1$error.rate
MyPerf.plsda1$error.rate.class
MyPerf.plsda1$choice.ncomp
g<-MyPerf.plsda1$predict         # Prediction values for each component

## 2 - The number of variables keepX

list.keepX <- seq(1, ncol(X),1)
head(list.keepX)


cl <- makePSOCKcluster(6)                                    #  nCores=6 ??? For parallelism
registerDoParallel(cl)

tune.splsda <- tune.splsda(X=X.sPLS, Y=Y, ncomp = 4, # ncomp = 3 as suggested in the previous step
                           validation = 'Mfold',
                           folds = 5, dist = 'max.dist',
                           measure = "BER", test.keepX = list.keepX,
                           nrepeat = 50,
                           progressBar = T)   # we suggest nrepeat = 50

error <- tune.splsda$error.rate
error
ncomp <- tune.splsda$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda$choice.keepX  # optimal number of variables to select
select.keepX  #     

plot(tune.splsda)

##------ sPLS-DA after tunning-----------------

splsda.tunned <- splsda(X = X.sPLS, Y = Y,
                        mode = "classic", 
                        ncomp = ncomp, max.iter=100,
                        keepX = select.keepX[1:ncomp])    # select.keepX but instead of 31 proteins in comp1 we selected 50

spls1 <- plotIndiv(splsda.tunned, comp=c(1,2),ind.names = T, 
                  ellipse =T, ellipse.level=.95,
                  star = F, legend=T, group = Y,
                  style = 'ggplot2',cex=4,
                  title = "sPLS-DA")

#---------Stability

perf.obj <- perf(splsda.tunned, 
                 dist = "max.dist",
                 validation = "Mfold",
                 folds = 5, nrepeat =50, 
                 auc = F, progressBar = TRUE)
perf.obj$error.rate
plot(perf.obj)


# We can also examine the stability of the variables selected across 
# the different cross-validation folds. Each variable's stability that is 
# selected across the CV runs is represented with a vertical bar. We often 
# observe a decrease in stability when more components are added in the model.

par(mfrow=c(1,3))
plot(perf.obj$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.obj$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)

# plot(perf.obj$features$stable[[3]], type = 'h', ylab = 'Stability', 
#      xlab = 'Features', main = 'Comp 3', las =2)

par(mfrow=c(1,1))

# The selectVar function outputs the selected variables along with their loading weight value
# here we match the selected variables to the stable features

#comp1____selectVar and stability
ind.match = match(mixOmics::selectVar(splsda.tunned, comp = 1)$name, 
                  names(perf.obj$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq1 = as.numeric(perf.obj$features$stable[[1]][ind.match])

df.selectVar_C1 <- data.frame(mixOmics::selectVar(splsda.tunned, comp = 1)$value, Freq1)
df.selectVar_C1$protein <- rownames(df.selectVar_C1)
df.selectVar_C1$comp <- "Comp1"

df.selectVar_C1.F <- subset(df.selectVar_C1, df.selectVar_C1$Freq1>=0.7)
dim(df.selectVar_C1.F)


#comp2____selectVar and stability

ind.match = match(mixOmics::selectVar(splsda.tunned, comp = 2)$name, 
                  names(perf.obj$features$stable[[2]]))
#extract the frequency of selection of those selected variables
Freq2 = as.numeric(perf.obj$features$stable[[2]][ind.match])

df.selectVar_C2 <- data.frame(mixOmics::selectVar(splsda.tunned, comp = 2)$value, Freq2)
df.selectVar_C2$protein <- rownames(df.selectVar_C2)
df.selectVar_C2$comp <- "Comp2"

df.selectVar_C2.F <- subset(df.selectVar_C2, df.selectVar_C2$Freq2>=0.7)


dim(df.selectVar_C2.F)

#comp3____selectVar and stability

  # ind.match = match(mixOmics::selectVar(splsda.tunned, comp = 3)$name, 
  #                   names(perf.obj$features$stable[[3]]))
  # #extract the frequency of selection of those selected variables
  # Freq3 = as.numeric(perf.obj$features$stable[[3]][ind.match])
  # 
  # df.selectVar_C3 <- data.frame(mixOmics::selectVar(splsda.tunned, comp = 3)$value, Freq3)
  # df.selectVar_C3$protein <- rownames(df.selectVar_C3)
  # df.selectVar_C3$comp <- "Comp3"
  # 
  # df.selectVar_C3.F <- subset(df.selectVar_C3, df.selectVar_C3$Freq3>=0.7)
  # dim(df.selectVar_C3.F)


df.selectVar <- plyr::join_all(list(df.selectVar_C1.F,
                                    df.selectVar_C2.F),
                               by="protein",
                               match = "all",
                               type = "full")

selectVar <- df.selectVar$protein
length(selectVar)

#df.selectVar <- rbind(df.selectVar_C1.F,df.selectVar_C2.F)

# saving matrixs

loading.matrix.X <- as.data.frame(splsda.tunned$loadings$X)
loading.matrix.X$protein <- rownames(loading.matrix.X)

vip.feature <- as.data.frame(vip(splsda.tunned))
colnames(vip.feature) <- paste("VIP_", names(vip.feature),sep = "")
vip.feature$protein <- row.names(vip.feature)

sPLS.DA_result <- plyr::join_all(list(loading.matrix.X, 
                                      vip.feature,
                                      df.selectVar_C1,
                                      df.selectVar_C2),
                                 by="protein")

rownames(sPLS.DA_result) <- sPLS.DA_result$protein

Big.table$protein <- rownames(Big.table)
Big.table.VV <- Big.table[,c("ValidV_C","ValidV_PCOS","protein")]
  
Total_results <- plyr::join_all(list(DEP,
                                     Big.table.VV,
                                     sPLS.DA_result),
                                 by="protein")
#write.csv(Total_results,"Total_results_non-PCOS vs PCOS2_210823.csv")

# 2D plot
plotIndiv(splsda.tunned, ind.names = F, ellipse =T,ellipse.level=.95,star = F, legend=T,
          style = 'ggplot2',cex=1,title = "sPLS-DA")

#The plot Loading function displays the loading weights, 
mixOmics::plotLoadings(splsda.tunned, comp = 1, 
             title = 'Loadings on comp 1', 
             size.name = 0.6,
             contrib = 'max', 
             method = 'mean',
             ndisplay=select.keepX[3])


# Heatmap sPLS-DA

col.ID <- randomcoloR::distinctColorPalette(k=10)[as.factor(Annotations$Patient)]
col.groups <- color.mixo(1:2)[as.factor(Annotations$Cond)]

cim.results <- cim(splsda.tunned,comp=c(1:2),#threshold = 0.9,
                   margins = c(4,10),zoom=F,keysize = c(1,1),transpose = T,
                   row.sideColors = cbind(col.groups, col.ID),
                   dist.method = c("euclidean","euclidean"),
                   clust.method=c("ward","ward"),
                   row.cex = 0.7,col.cex = 0.8,
                   row.names = Annotations$ID,
                   col.names = T, legend=list(legend = c(levels(as.factor(Annotations$Cond))), 
                                                   col = color.mixo(1:3),
                                                   title = "Condition", cex = 0.8))
cim.matrix <- cim.results$mat
cim.mat.short <- cim.matrix[ ,colnames(cim.matrix) %in% selectVar]
dim(cim.mat.short)

#write.csv(cim.mat.short,"cim.mat.short_sPLS-DA_CvsPCOS.csv")

# pheatmap package plot
library(pheatmap)

annotation_col <- as.data.frame(Annotations[,c("Patient","Cond")])
rownames(annotation_col)<-rownames(Annotations)
colnames(annotation_col) <- c("Patient","Cond")

ann_colors = list(Cond = c(`non-PCOS`="deepskyblue3", PCOS="darkorange2"))

#____________

col3<- colorRampPalette(c("dodgerblue4","dodgerblue3", "white","coral2","brown3"))(100)

hm<-pheatmap(t(cim.mat.short),  color = col3,
             cutree_cols = 2,
             cutree_rows = 2, show_rownames = T, 
             cluster_rows = T,
             cluster_cols = T, show_colnames = T,
             clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
             legend = T, 
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             annotation_legend = T,
             border_color="azure2",
             fontsize_row = 6.0,
             fontsize_col = 7.0,
             treeheight_row = 0, breaks = round(seq(-3,3,length=101 ),2)) #, filename = "test.pdf", cellwidth = 10, cellheight = 12
dev.off()
#====
dendog <- as.dendrogram(cim.results$ddc)
plot(dendog)
dendog1 <- as.hclust(cim.results$ddc)
clusters <-as.data.frame(cutree(dendog1, h=3))
colnames(clusters) <- "cluster"
clusters$protein <- row.names(clusters)


#####################################################################
                # small follicle vs Big follicles #
#####################################################################

# Patient 827 should be removed because she has 2 follicles with the same sizes.

df.log2.N1a <- df.log2.N1 %>% dplyr::select(!contains("C827"))

Annotations1 <- Annotations[!Annotations$Patient=="C827",]

#== (1) UNIVARIATE ANALYSIS --------(t-test paired)----------Foll 1 vs Foll 2---

df.ttest.Foll <- df.log2.N1a

groups.ttest <- Annotations1

GR <- as.factor(groups.ttest$group_Size2)
pat.GR<-as.factor(groups.ttest$Patient) 

####    Function DEP_prot      ####
# m = matrix with samples as columns and proteins as rows
# annotation.table = Annotation matrix
# groups = name of the column that contain the groups that should be compared (factor)  
# cofounders =  vector of characteres with the names of the columns that contain the cofounders variables (ex: c("Diameter", "sex","age").
# n_pat = number of patients

# m=df.ttest.Foll
# annotation.table = subset(Annotations1,Annotations1$Cond=="non-PCOS")
# groups = "group_Size2"
# n_pat=4
# id_ind = "Patient"
# cofounders <- c("Diameter")
# i=385
# i=1

DEP_ttest.Foll <- function(m, annotation.table, n_pat, groups) {

  groups1 <- as.factor(annotation.table[,groups])
  
  #--- making tables where the results will be stored---
  
  ttest <- matrix(nrow=nrow(m),ncol = 5, dimnames = list(rownames(m),c("p.val_ttest","Low.L", "Upp.L",
                                                                       "mean of the differences","DF")))
  #-----
  
  m1 = t(m)
  groups2 <- annotation.table %>% dplyr::select(all_of(c(groups)))
  groups2$sample <- rownames(groups2)
  
  for (i in 1:ncol(m1)) {
    
    prot.data <- as.data.frame(m1[,i])
    
    prot.data$sample <- rownames(prot.data)
    
    data <- plyr::join_all(list(groups2,prot.data),by="sample")
    row.names(data) <- data$sample
    colnames(data)[3] <- "protein"
    data <- data[,-2]
    data <- na.omit(data)
    data[,groups] <- as.factor(data[,groups])
    
    #====
    g1 <- as.data.frame(subset(data,data[,groups]==levels(data[,groups])[1]))
    g2 <- as.data.frame(subset(data,data[,groups]==levels(data[,groups])[2]))
    
    
    if(nrow(g1)==n_pat && nrow(g1)==nrow(g2)){
    
    
    ttest1 <- t.test(protein ~ group_Size2, data = data, 
                     paired = T, var.equal = FALSE,
                     alternative = "two.sided")
    
    ttest[i,1] <- ttest1$p.value                         # p-value
    ttest[i,2] <- ttest1$conf.int[1]                     # DL
    ttest[i,3] <- ttest1$conf.int[2]                     # UppL
    ttest[i,4] <- ttest1$estimate[1]   # mean of the differences
    ttest[i,5] <- ttest1$parameter[1]  # degree of freedom (Df)
    }
    else {ttest[i,] <- "N/S"}
  }
  ttest.Foll <- as.data.frame(ttest)
  
  ttest.Foll$p.val_ttest <- as.character(ttest.Foll$p.val_ttest)
  ttest.Foll$p.val_ttest <- as.numeric(ttest.Foll$p.val_ttest)
  ttest.Foll$q.value_ttest <- p.adjust(ttest.Foll$p.val_ttest, method = "fdr")
  
  
  
  #ttest$q.v.sig <- cut(ttest$q.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("p<.001", "p<.01", "p<.05","N/S"))
  
  return(ttest.Foll)
}

# t-Test Foll1 vs Foll2 --non-PCOS
DEP_ttest.Foll.non.PCOS <- DEP_ttest.Foll(m=df.ttest.Foll,n_pat=4,
                                               annotation.table=subset(Annotations1,Annotations1$Cond=="non-PCOS"),
                                               groups="group_Size2")
colnames(DEP_ttest.Foll.non.PCOS) <- paste(colnames(DEP_ttest.Foll.non.PCOS),".C",sep = "")
DEP_ttest.Foll.non.PCOS$protein <- rownames(DEP_ttest.Foll.non.PCOS)


# t-Test Foll1 vs Foll2 --PCOS
DEP_ttest.Foll.PCOS <- DEP_ttest.Foll(m=df.ttest.Foll,n_pat=5,
                                               annotation.table=subset(Annotations1,Annotations1$Cond=="PCOS"),
                                               groups="group_Size2")
colnames(DEP_ttest.Foll.PCOS) <- paste(colnames(DEP_ttest.Foll.PCOS),".pcos",sep = "")
DEP_ttest.Foll.PCOS$protein <- rownames(DEP_ttest.Foll.PCOS)


Total_results1 <- plyr::join_all(list(Total_results,
                                      DEP_ttest.Foll.non.PCOS,
                                      DEP_ttest.Foll.PCOS),
                                by="protein")

#write.csv(Total_results1, "Total_results1_Foll1vsFoll2_210823.csv")


#== (2) MULTIVARIATE ANALYSIS ---(mixOmics package)---Foll1 vs Foll2------non-PCOS

Annotations1a <- Annotations1[Annotations1$Cond=="non-PCOS",]
prot.matrix <- df.log2.N1a %>% dplyr::select(contains("C"))
names(prot.matrix)

# This analysis is multilevel

library(mixOmics)

X <- as.data.frame(t(prot.matrix))
Y <- as.factor(Annotations1a$group_Size2)

# Patient indicates the repeated measurements
# setup the design matrix by indicating the repeated measurements
design <- data.frame(Patient = Annotations1a$Patient)


# (1.1)----PCA-------------------------------------
## Apply filtering to keep proteins with 100% of valid values
# X.pca0 <- as.data.frame(t(X))
# X.pca = filter_valids(X.pca0,
#                       conditions = c("C", "P"),
#                       min_count = c(10, 10),
#                       at_least_one = F)
# X.pca <- filter(X.pca, KEEP)
# X.pca <- dplyr::select(X.pca, -c(KEEP,C,P))
# X.pca <- as.data.frame(t(X.pca))

X.pca <- X
# pca
pca.multilevel <- mixOmics::pca(X.pca, ncomp = 3, 
                                scale = F, center = F, 
                                multilevel = design)
pca

#2D plot per condition (non-PCOS vs PCOS)
pca1 <- mixOmics::plotIndiv(pca.multilevel, ind.names = T, legend = T,comp=c(1,2),
                            ellipse = T,cex = 3.5, style="ggplot2",
                            group = Y, 
                            title = 'PCA multilevel-PCOS')
pca1
#3D plot
# mixOmics::plotIndiv(pca.multilevel, ind.names = F, legend = T,comp=c(1,2,3),
#                     ellipse = F,cex = 1, style="3d",
#                     group = Annotations$group_Size,
#                     title = 'PCA multilevel')

# (1.3) ------sPLS-DA model -----------------------------

## Apply filtering to keep proteins with 100% of valid values

X.sPLS0 <- as.data.frame(t(X))
X.sPLS = filter_valids(X.sPLS0,
                      conditions = c("_1", "_2"),
                      min_count = c(4,4),
                      at_least_one = F)
X.sPLS <- filter(X.sPLS, KEEP)
X.sPLS <- dplyr::select(X.sPLS, -c(KEEP,"_1","_2"))
boxplot(X.sPLS)
X.sPLS.F <- as.data.frame(t(X.sPLS))
dim(X.sPLS.F)

#X.sPLS <- X

# sPLS-DA ------before Tuning----------
splsda.F.multilevel <- splsda(X.sPLS.F, Y = Y,ncomp = 2,
                              multilevel = design)

plotIndiv(splsda.F.multilevel, ind.names = F, 
          ellipse =T,ellipse.level=.95,
          star = F, legend=T,
          style = 'ggplot2',
          cex=1,title = "sPLS-DA multilevel")
#----Tuning parameters and numerical outputs___sPLS-DA

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat

## 1 - The number of components to retain

cl <- makePSOCKcluster(6)                                    #  nCores=6 ??? For parallelism
registerDoParallel(cl)
MyPerf.plsda.F <- perf(splsda.F.multilevel, 
                       validation = "Mfold", folds = 4, dist="all",         # folds :: https://machinelearningmastery.com/k-fold-cross-validation/#:~:text=Cross%2Dvalidation%20is%20a%20resampling,k%2Dfold%20cross%2Dvalidation. 
                       progressBar = T, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda.F, col = color.mixo(5:7), 
     sd = TRUE, legend.position = "horizontal")

MyPerf.plsda.F$error.rate
MyPerf.plsda.F$error.rate.class
MyPerf.plsda.F$choice.ncomp
g<-MyPerf.plsda.F$predict         # Prediction values for each component

## 2 - The number of variables keepX

list.keepX <- seq(1, ncol(X.sPLS.F),5)
head(list.keepX)


cl <- makePSOCKcluster(8)                                    #  nCores=6 ??? For parallelism
registerDoParallel(cl)

tune.splsda.F <- tune.splsda(X=X.sPLS.F, Y=Y, ncomp = 2, # ncomp = 3 as suggested in the previous step
                             multilevel = design,
                             validation = 'Mfold',
                             folds = 4, dist = 'max.dist',
                             measure = "BER", test.keepX = list.keepX,
                             nrepeat = 50,
                             progressBar = T)   # we suggest nrepeat = 50

error <- tune.splsda.F$error.rate
error
ncomp <- tune.splsda.F$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX.F <- tune.splsda.F$choice.keepX  # optimal number of variables to select
select.keepX.F  #    56  16  1  PCOS    #    56  16  1  non-PCOS 

plot(tune.splsda.F)

##------ sPLS-DA after tunning-----------------

splsda.tunned.F.multilevel <- splsda(X = X.sPLS.F, Y = Y,
                                     multilevel = design,
                                     ncomp = 2,
                                     keepX = c(50,50))   #select.keepX[1:2]

#---------Stability

perf.obj.F <- perf(splsda.tunned.F.multilevel, 
                   dist = "max.dist",
                   validation = "Mfold",
                   folds = 4, nrepeat =50, 
                   auc = F, progressBar = TRUE)
perf.obj.F$error.rate
plot(perf.obj.F)


# We can also examine the stability of the variables selected across 
# the different cross-validation folds. Each variable's stability that is 
# selected across the CV runs is represented with a vertical bar. We often 
# observe a decrease in stability when more components are added in the model.

par(mfrow=c(1,2))
plot(perf.obj.F$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.obj.F$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)

par(mfrow=c(1,1))

# The selectVar function outputs the selected variables along with their loading weight value
# here we match the selected variables to the stable features

#comp1____selectVar and stability
ind.match.F = match(mixOmics::selectVar(splsda.tunned.F.multilevel, comp = 1)$name, 
                  names(perf.obj.F$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq1.F = as.numeric(perf.obj.F$features$stable[[1]][ind.match.F])

df.selectVar_C1.F.multilevel <- data.frame(mixOmics::selectVar(splsda.tunned.F.multilevel, comp = 1)$value, Freq1.F)
df.selectVar_C1.F.multilevel$protein <- rownames(df.selectVar_C1.F.multilevel)
df.selectVar_C1.F.multilevel$comp <- "Comp1"

df.selectVar_C1.F.multilevel1 <- subset(df.selectVar_C1.F.multilevel, 
                                        df.selectVar_C1.F.multilevel$Freq1.F>=0.5)
dim(df.selectVar_C1.F.multilevel1)


#comp2____selectVar and stability

ind.match.F = match(mixOmics::selectVar(splsda.tunned.F.multilevel, comp = 2)$name, 
                  names(perf.obj.F$features$stable[[2]]))
#extract the frequency of selection of those selected variables
Freq2.F = as.numeric(perf.obj.F$features$stable[[1]][ind.match.F])

df.selectVar_C2.F.multilevel <- data.frame(mixOmics::selectVar(splsda.tunned.F.multilevel, comp = 2)$value, Freq2.F)
df.selectVar_C2.F.multilevel$protein <- rownames(df.selectVar_C2.F.multilevel)
df.selectVar_C2.F.multilevel$comp <- "Comp2"


df.selectVar_C2.F.multilevel1 <- subset(df.selectVar_C2.F.multilevel, 
                                        df.selectVar_C2.F.multilevel$Freq2.F>=0.5)
dim(df.selectVar_C2.F.multilevel1)


df.selectVar.F <- plyr::join_all(list(df.selectVar_C1.F.multilevel1,
                                      df.selectVar_C2.F.multilevel1),
                                 by="protein",
                                 match = "all",
                                 type = "full")

selectVar.F <- rownames(df.selectVar.F)
length(selectVar.F)

#df.selectVar <- rbind(df.selectVar_C1.F,df.selectVar_C2.F)

# saving matrixs

loading.matrix.X.F <- as.data.frame(splsda.tunned.F.multilevel$loadings$X)
loading.matrix.X.F$protein <- rownames(loading.matrix.X.F)

vip.feature.F <- as.data.frame(vip(splsda.tunned.F.multilevel))
colnames(vip.feature.F) <- paste("VIP_", names(vip.feature.F),sep = "")
vip.feature.F$protein <- row.names(vip.feature.F)

sPLS.DA_result.F <- plyr::join_all(list(loading.matrix.X.F, 
                                      vip.feature.F,
                                      df.selectVar_C1.F.multilevel,
                                      df.selectVar_C2.F.multilevel),
                                 by="protein")

rownames(sPLS.DA_result.F) <- sPLS.DA_result.F$protein

Total_results2 <- plyr::join_all(list(Total_results1,
                                      sPLS.DA_result.F),
                                 by="protein")

#write.csv(Total_results2,"Total_results2_Foll1vsFoll2a_210823a.csv")

# 2D plot
plotIndiv(splsda.tunned.F.multilevel, ind.names = F, ellipse =T,ellipse.level=.95,star = F, legend=T,
          style = 'ggplot2',cex=1,title = "sPLS-DA")

#The plot Loading function displays the loading weights, 
mixOmics::plotLoadings(splsda.tunned.F.multilevel, comp = 2, 
                       title = 'Loadings on comp 1 (non-PCOS vs PCOS)', 
                       size.name = 0.6,
                       contrib = 'max', 
                       method = 'mean',
                       ndisplay=45)


# Heatmap sPLS-DA

col.ID <- randomcoloR::distinctColorPalette(k=10)[as.factor(Annotations$Patient)]
col.groups <- color.mixo(1:2)[as.factor(Annotations$Cond)]

cim.results <- cim(splsda.tunned.F.multilevel,#comp=c(1:2),
                   margins = c(4,10),zoom=F,keysize = c(1,1),transpose = T,
                   row.sideColors = cbind(col.groups, col.ID),
                   dist.method = c("euclidean","euclidean"),
                   clust.method=c("ward","ward"),
                   row.cex = 0.8,col.cex = 0.8,
                   row.names = Annotations$ID,
                   col.names = T, legend=list(legend = c(levels(as.factor(Annotations$Cond))), 
                                              col = color.mixo(1:3),
                                              title = "Condition", cex = 0.8))
cim.matrix <- cim.results$mat
cim.mat.short <- cim.matrix[ ,colnames(cim.matrix) %in% selectVar]
dim(cim.mat.short)

#write.csv(cim.mat.short,"cim.mat.short_sPLS-DA_CvsPCOS.csv")

# pheatmap package plot
library(pheatmap)

annotation_col <- as.data.frame(Annotations[,c("Patient","Cond")])
rownames(annotation_col)<-rownames(Annotations)
colnames(annotation_col) <- c("Patient","Cond")

ann_colors = list(Cond = c(`non-PCOS`="deepskyblue3", PCOS="darkorange2"))

#____________

col3<- colorRampPalette(c("dodgerblue4","dodgerblue3", "white","coral2","brown3"))(100)

hm<-pheatmap(t(cim.mat.short),  color = col3,
             cutree_cols = 2,
             cutree_rows = 2, show_rownames = T, 
             cluster_rows = T,
             cluster_cols = T, show_colnames = T,
             clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
             legend = T, 
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             annotation_legend = T,
             border_color="azure2",
             fontsize_row = 6.0,
             fontsize_col =7.0,
             treeheight_row = 0, breaks = round(seq(-3,3,length=101 ),2)) #, filename = "test.pdf", cellwidth = 10, cellheight = 12
dev.off()
#====
dendog <- as.dendrogram(cim.results$ddc)
plot(dendog)
dendog1 <- as.hclust(cim.results$ddc)
clusters <-as.data.frame(cutree(dendog1, h=3))
colnames(clusters) <- "cluster"
clusters$protein <- row.names(clusters)


#########################################################################################

#== (2) MULTIVARIATE ANALYSIS ---(mixOmics package)---Foll1 vs Foll2------PCOS

Annotations1b <- Annotations1[Annotations1$Cond=="PCOS",]
prot.matrix2 <- df.log2.N1a %>% dplyr::select(contains("P"))
names(prot.matrix2)

# This analysis is multilevel

library(mixOmics)

X <- as.data.frame(t(prot.matrix2))
Y <- as.factor(Annotations1b$group_Size2)

# Patient indicates the repeated measurements
# setup the design matrix by indicating the repeated measurements
design <- data.frame(Patient = Annotations1b$Patient)


# (1.1)----PCA-------------------------------------
## Apply filtering to keep proteins with 100% of valid values
# X.pca0 <- as.data.frame(t(X))
# X.pca = filter_valids(X.pca0,
#                       conditions = c("C", "P"),
#                       min_count = c(10, 10),
#                       at_least_one = F)
# X.pca <- filter(X.pca, KEEP)
# X.pca <- dplyr::select(X.pca, -c(KEEP,C,P))
# X.pca <- as.data.frame(t(X.pca))

X.pca <- X
# pca
pca.multilevel <- mixOmics::pca(X.pca, ncomp = 3, 
                                scale = F, center = F, 
                                multilevel = design)
pca

#2D plot per condition (non-PCOS vs PCOS)
pca1 <- mixOmics::plotIndiv(pca.multilevel, ind.names = T, legend = T,comp=c(1,2),
                            ellipse = T,cex = 3.5, style="ggplot2",
                            group = Y, 
                            title = 'PCA multilevel-PCOS')
pca1
#3D plot
# mixOmics::plotIndiv(pca.multilevel, ind.names = F, legend = T,comp=c(1,2,3),
#                     ellipse = F,cex = 1, style="3d",
#                     group = Annotations$group_Size,
#                     title = 'PCA multilevel')

# (1.3) ------sPLS-DA model -----------------------------

## Apply filtering to keep proteins with 100% of valid values

X.sPLS0 <- as.data.frame(t(X))
X.sPLS = filter_valids(X.sPLS0,
                       conditions = c("_1", "_2"),
                       min_count = c(5,5),
                       at_least_one = F)
X.sPLS <- filter(X.sPLS, KEEP)
X.sPLS <- dplyr::select(X.sPLS, -c(KEEP,"_1","_2"))
boxplot(X.sPLS)
X.sPLS.F <- as.data.frame(t(X.sPLS))
dim(X.sPLS.F)

#X.sPLS <- X

# sPLS-DA ------before Tuning----------
splsda.F.multilevel <- splsda(X.sPLS.F, Y = Y,ncomp = 2,
                              multilevel = design)

plotIndiv(splsda.F.multilevel, ind.names = F, 
          ellipse =T,ellipse.level=.95,
          star = F, legend=T,
          style = 'ggplot2',
          cex=1,title = "sPLS-DA multilevel")
#----Tuning parameters and numerical outputs___sPLS-DA

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat

## 1 - The number of components to retain

cl <- makePSOCKcluster(8)                                    #  nCores=6 ??? For parallelism
registerDoParallel(cl)
MyPerf.plsda.F <- perf(splsda.F.multilevel, 
                       validation = "Mfold", folds = 4, dist="all",         # folds :: https://machinelearningmastery.com/k-fold-cross-validation/#:~:text=Cross%2Dvalidation%20is%20a%20resampling,k%2Dfold%20cross%2Dvalidation. 
                       progressBar = T, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda.F, col = color.mixo(5:7), 
     sd = TRUE, legend.position = "horizontal")

MyPerf.plsda.F$error.rate
MyPerf.plsda.F$error.rate.class
MyPerf.plsda.F$choice.ncomp
g<-MyPerf.plsda.F$predict         # Prediction values for each component

## 2 - The number of variables keepX

list.keepX <- seq(1, ncol(X.sPLS.F),5)
head(list.keepX)


cl <- makePSOCKcluster(8)                                    #  nCores=6 ??? For parallelism
registerDoParallel(cl)

tune.splsda.F <- tune.splsda(X=X.sPLS.F, Y=Y, ncomp = 2, # ncomp = 3 as suggested in the previous step
                             multilevel = design,
                             validation = 'Mfold',
                             folds = 4, dist = 'max.dist',
                             measure = "BER", test.keepX = list.keepX,
                             nrepeat = 50,
                             progressBar = T)   # we suggest nrepeat = 50

error <- tune.splsda.F$error.rate
error
ncomp <- tune.splsda.F$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX.F <- tune.splsda.F$choice.keepX  # optimal number of variables to select
select.keepX.F  #    56  16  1  PCOS    #    56  16  1  non-PCOS 

plot(tune.splsda.F)

##------ sPLS-DA after tunning-----------------

splsda.tunned.F.multilevel <- splsda(X = X.sPLS.F, Y = Y,
                                     multilevel = design,
                                     ncomp = 2,
                                     keepX = c(50,50))   #select.keepX[1:2]

#---------Stability

perf.obj.F <- perf(splsda.tunned.F.multilevel, 
                   dist = "max.dist",
                   validation = "Mfold",
                   folds = 4, nrepeat =50, 
                   auc = F, progressBar = TRUE)
perf.obj.F$error.rate
plot(perf.obj.F)


# We can also examine the stability of the variables selected across 
# the different cross-validation folds. Each variable's stability that is 
# selected across the CV runs is represented with a vertical bar. We often 
# observe a decrease in stability when more components are added in the model.

par(mfrow=c(1,2))
plot(perf.obj.F$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.obj.F$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)

par(mfrow=c(1,1))

# The selectVar function outputs the selected variables along with their loading weight value
# here we match the selected variables to the stable features

#comp1____selectVar and stability
ind.match.F = match(mixOmics::selectVar(splsda.tunned.F.multilevel, comp = 1)$name, 
                    names(perf.obj.F$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq1.F = as.numeric(perf.obj.F$features$stable[[1]][ind.match.F])

df.selectVar_C1.F.multilevel <- data.frame(mixOmics::selectVar(splsda.tunned.F.multilevel, comp = 1)$value, Freq1.F)
df.selectVar_C1.F.multilevel$protein <- rownames(df.selectVar_C1.F.multilevel)
df.selectVar_C1.F.multilevel$comp <- "Comp1"

df.selectVar_C1.F.multilevel1 <- subset(df.selectVar_C1.F.multilevel, 
                                        df.selectVar_C1.F.multilevel$Freq1.F>=0.5)
dim(df.selectVar_C1.F.multilevel1)


#comp2____selectVar and stability

ind.match.F = match(mixOmics::selectVar(splsda.tunned.F.multilevel, comp = 2)$name, 
                    names(perf.obj.F$features$stable[[2]]))
#extract the frequency of selection of those selected variables
Freq2.F = as.numeric(perf.obj.F$features$stable[[1]][ind.match.F])

df.selectVar_C2.F.multilevel <- data.frame(mixOmics::selectVar(splsda.tunned.F.multilevel, comp = 2)$value, Freq2.F)
df.selectVar_C2.F.multilevel$protein <- rownames(df.selectVar_C2.F.multilevel)
df.selectVar_C2.F.multilevel$comp <- "Comp2"


df.selectVar_C2.F.multilevel1 <- subset(df.selectVar_C2.F.multilevel, 
                                        df.selectVar_C2.F.multilevel$Freq2.F>=0.5)
dim(df.selectVar_C2.F.multilevel1)


df.selectVar.F <- plyr::join_all(list(df.selectVar_C1.F.multilevel1,
                                      df.selectVar_C2.F.multilevel1),
                                 by="protein",
                                 match = "all",
                                 type = "full")

selectVar.F <- rownames(df.selectVar.F)
length(selectVar.F)

#df.selectVar <- rbind(df.selectVar_C1.F,df.selectVar_C2.F)

# saving matrixs

loading.matrix.X.F <- as.data.frame(splsda.tunned.F.multilevel$loadings$X)
loading.matrix.X.F$protein <- rownames(loading.matrix.X.F)

vip.feature.F <- as.data.frame(vip(splsda.tunned.F.multilevel))
colnames(vip.feature.F) <- paste("VIP_", names(vip.feature.F),sep = "")
vip.feature.F$protein <- row.names(vip.feature.F)

sPLS.DA_result.F <- plyr::join_all(list(loading.matrix.X.F, 
                                        vip.feature.F,
                                        df.selectVar_C1.F.multilevel,
                                        df.selectVar_C2.F.multilevel),
                                   by="protein")

rownames(sPLS.DA_result.F) <- sPLS.DA_result.F$protein

Total_results3 <- plyr::join_all(list(Total_results2,
                                      sPLS.DA_result.F),
                                 by="protein")

#write.csv(Total_results3,"Total_results3_Foll1vsFoll2a_210823.csv")

# 2D plot
plotIndiv(splsda.tunned.F.multilevel, ind.names = F, ellipse =T,ellipse.level=.95,star = F, legend=T,
          style = 'ggplot2',cex=1,title = "sPLS-DA")

#The plot Loading function displays the loading weights, 
mixOmics::plotLoadings(splsda.tunned.F.multilevel, comp = 2, 
                       title = 'Loadings on comp 1 (non-PCOS vs PCOS)', 
                       size.name = 0.6,
                       contrib = 'max', 
                       method = 'mean',
                       ndisplay=45)


# Heatmap sPLS-DA

col.ID <- randomcoloR::distinctColorPalette(k=10)[as.factor(Annotations$Patient)]
col.groups <- color.mixo(1:2)[as.factor(Annotations$Cond)]

cim.results <- cim(splsda.tunned.F.multilevel,#comp=c(1:2),
                   margins = c(4,10),zoom=F,keysize = c(1,1),transpose = T,
                   row.sideColors = cbind(col.groups, col.ID),
                   dist.method = c("euclidean","euclidean"),
                   clust.method=c("ward","ward"),
                   row.cex = 0.8,col.cex = 0.8,
                   row.names = Annotations$ID,
                   col.names = T, legend=list(legend = c(levels(as.factor(Annotations$Cond))), 
                                              col = color.mixo(1:3),
                                              title = "Condition", cex = 0.8))
cim.matrix <- cim.results$mat
cim.mat.short <- cim.matrix[ ,colnames(cim.matrix) %in% selectVar]
dim(cim.mat.short)

#write.csv(cim.mat.short,"cim.mat.short_sPLS-DA_CvsPCOS.csv")

# pheatmap package plot
library(pheatmap)

annotation_col <- as.data.frame(Annotations[,c("Patient","Cond")])
rownames(annotation_col)<-rownames(Annotations)
colnames(annotation_col) <- c("Patient","Cond")

ann_colors = list(Cond = c(`non-PCOS`="deepskyblue3", PCOS="darkorange2"))

#____________

col3<- colorRampPalette(c("dodgerblue4","dodgerblue3", "white","coral2","brown3"))(100)

hm<-pheatmap(t(cim.mat.short),  color = col3,
             cutree_cols = 2,
             cutree_rows = 2, show_rownames = T, 
             cluster_rows = T,
             cluster_cols = T, show_colnames = T,
             clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
             legend = T, 
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             annotation_legend = T,
             border_color="azure2",
             fontsize_row = 6.0,
             fontsize_col =7.0,
             treeheight_row = 0, breaks = round(seq(-3,3,length=101 ),2)) #, filename = "test.pdf", cellwidth = 10, cellheight = 12
dev.off()
#====
dendog <- as.dendrogram(cim.results$ddc)
plot(dendog)
dendog1 <- as.hclust(cim.results$ddc)
clusters <-as.data.frame(cutree(dendog1, h=3))
colnames(clusters) <- "cluster"
clusters$protein <- row.names(clusters)

#======================
#== (1) UNIVARIATE ANALYSIS --------(Linear Model Using Generalized Least Squares)------
library(nlme)

m = df.log2.N1
annotation.table = Annotations
sample.ID="ID"
var1 = "Cond"
var2 = "group_Size"
patient.col = "Patient"


#--- making matrix where the results will be stored---

test <- matrix(nrow=nrow(m),ncol = 4, 
               dimnames = list(rownames(m),
                               c("Intercept",var1,var2, "Interaccion")))

#-----
#df.log2.N1
m1 = as.data.frame(t(as.data.frame(m)))


for (i in 1:ncol(m1)) {
  
  prot.data <- as.data.frame(m1[,i])
  
  prot.data$ID <- rownames(m1)
  
  
  data <- plyr::join_all(list(annotation.table[,c(sample.ID,var1,var2,patient.col)],prot.data),by=sample.ID)
  row.names(data) <- data$ID
  colnames(data)[5] <- "protein"
  data <- data[,-1]
  data <- na.omit(data)
  for (j in 1:ncol(data)) {if(is.character(data[,j])){data[,j] <- as.factor(data[,j])}}
  
  
  fit.compsym <- gls(protein ~ Cond * group_Size, data=data, 
                     corr=corCompSymm(form= ~ 1 | Patient))
  
  avo <- as.data.frame(anova(fit.compsym)) 
  
  test[i,] <- t(avo[,3])
  
  
  #plots
  # 
  # attach(data)
  # 
  # #--
  # png(filename = paste("plots/interaction plots/",colnames(m1)[i],"_interaction plot_",".png",sep=""),
  #     width = 500, height = 480, units = "px", pointsize = 16,
  #     bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
  # 
  # plot<-interaction.plot(group_Size, Cond, protein, fun = function(x) mean(x, na.rm = TRUE),
  #                        lty=c(1:2),lwd=2,
  #                        type = "c",
  #                        main=colnames(m1)[i],
  #                        ylab="Log2 Intensity (scaled)",
  #                        xlab="Follicular Size",
  #                        trace.label="")
  # print(plot)
  # dev.off()
  # 
  # #--
  # svg(filename = paste("plots/interaction plots/",colnames(m1)[i],"_interaction plot_",".svg",sep=""),
  #     width = 4, height = 4, pointsize = 12, onefile = F, family = "sans", bg = "white")
  #     #antialias = c("default", "none", "gray", "subpixel"),
  #     
  # 
  # plot<-interaction.plot(group_Size, Cond, protein, fun = function(x) mean(x, na.rm = TRUE),
  #                        lty=c(1:2),lwd=2,
  #                        main=colnames(m1)[i],
  #                        ylab="Log2 Intensity (scaled)",
  #                        xlab="Follicular Size",
  #                        trace.label="")
  # print(plot)
  # dev.off()
  # detach(data)
  # 
  #--boxplot
  
  attach(data)
  #tuk.sig <- cut(pv[i,], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("p<.001", "p<.01", "p<.05","N/S"))

  png(filename = paste("plots/boxplot/",colnames(m1)[i],"_boxplot_",".png",sep=""),width = 500, height = 480, units = "px", pointsize = 16,
      bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))

  #https://datavizpyr.com/how-to-make-boxplot-with-a-line-connecting-mean-values-in-r/#google_vignette
  
  df1_median <- data %>% 
    group_by(Cond) %>% 
    summarize(average = median(protein)) %>%
    ungroup()
  
  plot1 <- data %>% 
                ggplot(mapping = aes(x = Cond, y = protein, fill = Cond)) + 
                geom_boxplot(alpha=0.5,show.legend = FALSE) +
                geom_beeswarm(priority='density',cex=1.5, alpha=0.5,show.legend = FALSE)+
                theme_bw(base_size = 16)+
                geom_point(data = df1_median, show.legend = FALSE,
                           mapping = aes(x = Cond, y = average),
                           color="red", size=3) +
                geom_line(data = df1_median, colour="red",
                          mapping = aes(x = Cond, y = average, group=1))+
                labs(subtitle=colnames(m1)[i],
                     caption = "",
                     x = "", y= "Intensity (Log2 scaled)")+
                scale_fill_manual(breaks = Cond,na.value = "red",
                                  values = c("#1b98e0", "#1b99e0"))
              
  print(plot1)
  dev.off()
  detach(data)
  
}
test1 <- as.data.frame(test)
test1$Accession <- row.names(test)

result.m <- plyr::join_all(list(prot.exp,test1),by="Accession")

#write.csv(result.m,"results.csv")

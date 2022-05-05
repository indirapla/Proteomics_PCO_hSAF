
##=====Follicular fluid=====PCOS project=============================================================
# Indira Pla Parada 

#====INSTALL PACKAGES================================================================================

# List of packages to install
.packages = c("BiocManager","devtools","ggplot2","ggbiplot", "pca3d","pcadapt","outliers","igraph",
              "rgl","graphics","reshape2","dplyr","ggpubr","remotes",
              "FactoMineR", "factoextra","corrplot","ggpubr","fpc", "NbClust","mixOmics", 
              "phyloseq","lme4","nlme","car","plotly","RadaFDR","GOplot","tidyverse",
              "GOsummaries","randomcoloR","parallel","doParallel","afex","ggbeeswarm",
              "emmeans","psych", "pheatmap","lmerTest","nlme","Seurat","proBatch")
#"RFunrichWebService",
# Install packages if not installed
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  install_github("vqv/ggbiplot")
  install_github("fxia22/RadaFDR")
  BiocManager::install(c("mixOmics","phyloseq","GOsummaries","Seurat"))
  devtools::install_github('symbioticMe/proBatch', build_vignettes = TRUE)
}

# Loading packages
lapply(.packages, require, character.only=TRUE)


#=========================FUNCTIONS======================================================

# -------Data filtering -----------------
# This function filters proteins based on the amoung of missing values per sample/condition

filter_valids = function(df, conditions, min_count, at_least_one = FALSE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  
  log2.names = names(df)
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  values.count = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
  })
  colnames(values.count) <- conditions
  
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  
  df <- cbind(df,values.count)
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}


## ------Data Normalization (standardization) --------------

# Subtracting the sample median from each value

median_centering = function(df) {
  # df = data frame containing LOG2 columns for normalization
  LOG2.names = names(df)
  
  df[, LOG2.names] = lapply(LOG2.names, 
                            function(x) {
                              LOG2 = df[[x]]
                              LOG2[!is.finite(LOG2)] = NA   # Exclude missing values from median calculation
                              gMedian = median(LOG2, na.rm = TRUE)
                              LOG2 - gMedian
                            }
  )
  
  return(df)
}
# -------Multiple plot function--------------------------

# Multiple plot function http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/#:~:text=multiplot%20function%20This%20is%20the%20definition%20of%20multiplot.,a%20list%20of%20plot%20objects%20passed%20to%20plotlist.
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# -------Coeficient of variation (CV)-----------------------------------
# calculate CV

CV <- function(x){sd(x)/mean(x)*100}

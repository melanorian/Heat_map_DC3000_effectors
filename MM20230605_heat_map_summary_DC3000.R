# Name: Melanie Mendel
# Date created: 13-01-2023
# Date continues: 18-01-2023
# Updated with new data: 30-01-2023

# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# clean up
rm(list=ls())

# laod libraries
library(dplyr)
library(pheatmap)
library(ggplot2)
library(openxlsx)
library(plotly)
library(gridExtra)
library(svglite)
library(ggpubr)
library("factoextra")

# set working directory
setwd("C:/Users/Mende012/Documents/Bioinformatics/Heat_map_effectors ")

# load data

cfu <- read.csv("./Median_CFU_counts.csv") # change names in csv file wiht excel before loading, + causes trouble in changing names
ros <- read.xlsx ("./summary_negative_median_ROS_percentage.xlsx", startRow = 2) # ROS script give strange output that requires editing before loading
death <- read.xlsx("./MM20230601_disease_score.xlsx")

# match column names and names of effectors 
# effector lines: D36E, DC3000, HopA1,...
# heading: effector lines = Pst_strain


# modify CFU data
cfu <- cbind.data.frame("Pst_strain" = cfu$Pst_strain, "CFU" = cfu$median)
cfu[cfu == "D36"] <- "AAAA_D36E"
cfu[cfu == "DC3000"] <- "ZZZZ_DC3000"
cfu[cfu == "HopAA1_1"] <- "HopAA1-1"
cfu[cfu == "HopAA1_2"] <- "HopAA1-2"
cfu[cfu == "HopO1"] <- "HopO1-1"
cfu[cfu == "HopQ1"] <- "HopQ1-1"
cfu <- cfu[!(cfu$Pst_strain == "dCEL" ),]
cfu <- cfu[!(cfu$Pst_strain == "hrcC-" ),]

length(sort(cfu$Pst_strain))

# modify ROS data
ros <- cbind.data.frame("Pst_strain" = ros$Effector, "ROS" = ros$Median)
length(ros$Pst_strain)

# modify cell-death data
death <- data.frame("Pst_strain" = death$Effector, "Cell-death_positive" = death$Median)
length(death$Pst_strain)

# merge the three different data sets by row 
m1 <- merge(cfu, ros, all.x = TRUE) #by = "Pst_strain")
m2 <- merge(m1, death, all.x = TRUE)

# drop hrcC-/dCEL for final figure
m2 <- m2[!(m2$Pst_strain == "dCEL"),]
m2 <- m2[!(m2$Pst_strain == "hrcC-"),]

rownames(m2) <- m2$Pst_strain
m2 <- m2[,-1]

# convert data frame into matrix
matrix <- data.matrix(m2)
class(matrix)

# basic heatmap  interactive heatmap d3heatmap ??
p1 <- pheatmap(m2, # input matrix
         scale = "column", # scale by column since the input values are very different
         cluster_rows = T,
         cluster_cols = F,
         treeheight_row = F, # hide dendogramm of row
         treeheight_col = F,
         cellwidth = 20,
         cellheight = 20,
         fontsize = 12,
         border_color = "NA",
         filename = paste(pre, sep = "","_summary_heatmap_1.pdf"))

p1
dev.off()
         
# pdf to svg https://convertio.co/download/ee6526671d9b1efd74aa2da5f9bf6e5723ded5/

# experiment with clustered heatmap (acutally cut the heat map)
p2 <- pheatmap(m2, # input matrix
               scale = "column", # scale by column since the input values are very different
               cluster_rows = T,
               cluster_cols = F,
               cutree_rows=5,
               cellwidth = 20,
               cellheight = 20,
               fontsize = 12,
               border_color = "NA",
               filename = paste(pre, sep = "","_summary_heatmap_scaled_by_column_5_cluseters.pdf"))

p2
dev.off()

################################################################################
#######################  Comparisons betweeen experiments    ###################
################################################################################

# lm ROS - CFU
lm_r_c <- lm(CFU ~ ROS, data = m2)
lm_r_c_summary <- summary(lm_r_c)

sink(paste(pre, sep = "","statsummary_CFU_ROS_lm.txt"))
print(lm_r_c_summary)
sink()

# lm CFU - cell death
lm_c_d <- lm(m2$CFU ~ m2$Cell.death_positive, data = m2)
lm_c_d_summary <- summary(lm_c_d)

sink(paste(pre, sep = "","statsummary_CUF_celldeath_lm.txt"))
print(lm_c_d_summary)
sink()

# lm CFU - cell death
lm_r_d <- lm(m2$ROS ~ m2$Cell.death_positive, data = m2)
lm_r_d_summary <- summary(lm_r_d)

sink(paste(pre, sep = "","statsummary_ROS_celldeath_lm.txt"))
print(lm_r_d_summary)
sink()

# compare with scatter plot the different experimental outcomes
# use m2 df 

#"#E69F00", "#faebcc"

g1 <- ggplot(m2, aes(x=CFU, y=ROS)) + 
      geom_point() +
 
  geom_text(aes(label=rownames(m2))) +
  
  geom_smooth(method=lm) +
  stat_regline_equation(label.x=-0.5, label.y=0) +
  stat_cor(aes(label=..rr.label..), label.x=-0.5, label.y=-5) +
  theme_classic()  
  
ggplotly(g3) 


g2 <- ggplot(m2, aes(x=m2$CFU, y= m2$Cell.death_positive)) + 
  geom_point() +
  
  geom_text(aes(label=rownames(m2))) +
  
  stat_regline_equation(label.x=-0.5, label.y=0) +
  stat_cor(aes(label=..rr.label..), label.x=-0.5, label.y=-5) +
  geom_smooth(method=lm) +
  theme_classic()  

g2

g3 <- ggplot(m2, aes(x=m2$ROS, y= m2$Cell.death_positive)) + 
  geom_point() +
  
  geom_text(aes(label=rownames(m2))) +
  
  geom_smooth(method=lm) +
  stat_regline_equation(label.x=50, label.y=0) +
  stat_cor(aes(label=..rr.label..), label.x=3, label.y=-1) +
  theme_classic()  

g3

# multipanel figure
g4 <- grid.arrange(g1, g2, g3, nrow = 3)

# safe multipanel
ggsave(filename =  paste(pre, sep = "","_lm_figure_labeled.svg"), 
       plot = g4,
       device = "svg")


######## plot without label ###########

g5 <- ggplot(m2, aes(x=CFU, y=ROS)) + 
  geom_point() +
  geom_smooth(method=lm) +
  theme_classic()  

ggplotly(g5) 


g6 <- ggplot(m2, aes(x=m2$CFU, y= m2$Cell.death_positive)) + 
  geom_point() +
  geom_smooth(method=lm) +
  theme_classic()  

g6

g7 <- ggplot(m2, aes(x=m2$ROS, y= m2$Cell.death_positive)) + 
  geom_point() +
  geom_smooth(method=lm) +
  theme_classic()  

g7

# multipanel figure
g8 <- grid.arrange(g5, g6, g7, nrow = 3)

ggsave(filename =  paste(pre, sep = "","_lm_figure_nolabel.svg"), 
       plot = g8,
       device = "svg")

dev.off()

#################################################################################
####################   new heatmap with comparison columns     ##################
#################################################################################

### add additional columns to m2 that calculating the differendes between results
# calculate differences between the results
m2["Celldeath-ROS"] <- m2$Cell.death_positive - m2$ROS
m2 ["Celldeath-CFU"] <- m2$Cell.death_positive - m2$CFU
m2 ["CFU-ROS"] <- m2$CFU - m2$ROS

# rename the dataframe to avoid confustions
m3 <- m2

# convert data frame into matrix
matrix <- data.matrix(m3)
class(matrix)

# basic heatmap  interactive heatmap d3heatmap ??
p2 <- pheatmap(m3, # input matrix
               scale = "column", # scale by column since the input values are very different
               cluster_cols = F,
               cluster_rows = T,
               treeheight_row = 0, # hide dendogramm of row
               treeheight_col = 0,
               cutree_rows=5,
               cellwidth = 20,
               cellheight = 20,
               fontsize = 12,
               border_color = "NA",
               filename = paste(pre, sep = "","_summary_heatmap_comparison_columns_1.pdf")) #scaled_by_column

p2
dev.off()

# try methods to find optimal number of clusters for heat map

wss <- fviz_nbclust(m2, FUNcluster = kmeans, method = "wss")
sil <- fviz_nbclust(m2, FUNcluster = kmeans, method = "silhouette")
hcut <- fviz_nbclust(m2, FUNcluster = kmeans, method = "gap_stat")


gr_cluster <- grid.arrange(wss, sil, hcut, nrow = 3)
ggsave(filename = paste(pre, sep = "","optimal_nr_clusters_heatmap_DC3000_effectors_1.pdf"), 
       plot = gr_cluster, 
       device = "svg")


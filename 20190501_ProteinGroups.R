setwd("K:/Collaborations/Ethiraj_Ravindran/20190410_Ethiraj_NucFrac_Rep1-6/output")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(limma)
library(ggrepel)
library(VennDiagram)
library(gplots)
library(data.table)
library(GGally)







#### proteinGroups table load and preparation ####

PG <- fread("../txt/proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)

#filter out contaminants, reverse and only identified by site
PG <- subset(PG, Reverse != "+")
PG <- subset(PG, Potential.contaminant != "+")
PG <- subset(PG, Only.identified.by.site != "+")


IBAQ <- c("iBAQ.Control_Cytosol_Rep1",
          "iBAQ.Control_Cytosol_Rep2",
          "iBAQ.Control_Cytosol_Rep3",
          "iBAQ.Control_Cytosol_Rep4",
          "iBAQ.Control_Cytosol_Rep5",
          "iBAQ.Control_Cytosol_Rep6",
          
          "iBAQ.Patient_Cytosol_Rep1",
          "iBAQ.Patient_Cytosol_Rep2",
          "iBAQ.Patient_Cytosol_Rep3",
          "iBAQ.Patient_Cytosol_Rep4",
          "iBAQ.Patient_Cytosol_Rep5",
          "iBAQ.Patient_Cytosol_Rep6",
          
          "iBAQ.Control_Input_Rep1",
          "iBAQ.Control_Input_Rep2",
          "iBAQ.Control_Input_Rep3",
          "iBAQ.Control_Input_Rep4",
          "iBAQ.Control_Input_Rep5",
          "iBAQ.Control_Input_Rep6",
          
          "iBAQ.Patient_Input_Rep1",
          "iBAQ.Patient_Input_Rep2",
          "iBAQ.Patient_Input_Rep3",
          "iBAQ.Patient_Input_Rep4",
          "iBAQ.Patient_Input_Rep5",
          "iBAQ.Patient_Input_Rep6",
          
          "iBAQ.Control_Nucleus_Rep1",
          "iBAQ.Control_Nucleus_Rep2",
          "iBAQ.Control_Nucleus_Rep3",
          "iBAQ.Control_Nucleus_Rep4",
          "iBAQ.Control_Nucleus_Rep5",
          "iBAQ.Control_Nucleus_Rep6",
          
          "iBAQ.Patient_Nucleus_Rep1",
          "iBAQ.Patient_Nucleus_Rep2",
          "iBAQ.Patient_Nucleus_Rep3",
          "iBAQ.Patient_Nucleus_Rep4",
          "iBAQ.Patient_Nucleus_Rep5",
          "iBAQ.Patient_Nucleus_Rep6")



LFQ <- c("LFQ.intensity.Control_Cytosol_Rep1",
         "LFQ.intensity.Control_Cytosol_Rep2",
         "LFQ.intensity.Control_Cytosol_Rep3",
         "LFQ.intensity.Control_Cytosol_Rep4",
         "LFQ.intensity.Control_Cytosol_Rep5",
         "LFQ.intensity.Control_Cytosol_Rep6",
         
         "LFQ.intensity.Patient_Cytosol_Rep1",
         "LFQ.intensity.Patient_Cytosol_Rep2",
         "LFQ.intensity.Patient_Cytosol_Rep3",
         "LFQ.intensity.Patient_Cytosol_Rep4",
         "LFQ.intensity.Patient_Cytosol_Rep5",
         "LFQ.intensity.Patient_Cytosol_Rep6",
         
         "LFQ.intensity.Control_Input_Rep1",
         "LFQ.intensity.Control_Input_Rep2",
         "LFQ.intensity.Control_Input_Rep3",
         "LFQ.intensity.Control_Input_Rep4",
         "LFQ.intensity.Control_Input_Rep5",
         "LFQ.intensity.Control_Input_Rep6",
         
         "LFQ.intensity.Patient_Input_Rep1",
         "LFQ.intensity.Patient_Input_Rep2",
         "LFQ.intensity.Patient_Input_Rep3",
         "LFQ.intensity.Patient_Input_Rep4",
         "LFQ.intensity.Patient_Input_Rep5",
         "LFQ.intensity.Patient_Input_Rep6",
         
         "LFQ.intensity.Control_Nucleus_Rep1",
         "LFQ.intensity.Control_Nucleus_Rep2",
         "LFQ.intensity.Control_Nucleus_Rep3",
         "LFQ.intensity.Control_Nucleus_Rep4",
         "LFQ.intensity.Control_Nucleus_Rep5",
         "LFQ.intensity.Control_Nucleus_Rep6",
         
         "LFQ.intensity.Patient_Nucleus_Rep1",
         "LFQ.intensity.Patient_Nucleus_Rep2",
         "LFQ.intensity.Patient_Nucleus_Rep3",
         "LFQ.intensity.Patient_Nucleus_Rep4",
         "LFQ.intensity.Patient_Nucleus_Rep5",
         "LFQ.intensity.Patient_Nucleus_Rep6")




#transform LFQ and ratios to Log2 (or 10 if you prefer)
PG[c(IBAQ, LFQ)] = log2(PG[c(IBAQ, LFQ)])
# change Inf values for na
is.na(PG[c(IBAQ, LFQ)]) <- sapply(PG[c(IBAQ, LFQ)], is.infinite)
is.na(PG[c(IBAQ, LFQ)]) <- sapply(PG[c(IBAQ, LFQ)], is.nan)






PG$Gene.names <- sapply(strsplit(PG$Gene.names, ";"), "[", 1)
PG$Majority.protein.IDs <- sapply(strsplit(PG$Majority.protein.IDs, ";"), "[", 1)






# how many proteins were identified (have LFQ values not NA) in each L-H group pair
for (i in 1:length(LFQ)) {
  cat(LFQ[i])
  cat("\t")
  cat("\t")
  cat(nrow(PG[!is.na(PG[LFQ[i]]),]))
  cat("\t")
  cat(mean(PG[,LFQ[i]], na.rm = T))
  cat("\n")
  rm(i)
}

for (i in 1:length(IBAQ)) {
  cat(IBAQ[i])
  cat("\t")
  cat("\t")
  cat(nrow(PG[!is.na(PG[IBAQ[i]]),]))
  cat("\t")
  cat(mean(PG[,IBAQ[i]], na.rm = T))
  cat("\n")
  rm(i)
}







#### Define PG_summ as working table ####

PG_summ <- subset(PG, select = c("Majority.protein.IDs","Gene.names",
                                 LFQ, IBAQ,
                                 "id"))



PG_summ <- subset(PG_summ,
                  apply(PG_summ[LFQ], 1,
                        function(x) (sum(!is.na(x[1:6])) >= 2 & sum(!is.na(x[7:12])) >= 2) |
                          (sum(!is.na(x[13:18])) >= 2 & sum(!is.na(x[19:24])) >= 2) |
                          (sum(!is.na(x[25:30])) >= 2 & sum(!is.na(x[31:36])) >= 2)))

#### Fold Change calculation ####
#Cytosol
PG_summ$FC_LFQ_Cytosol_Rep1 <- PG_summ[,LFQ[7]]-PG_summ[,LFQ[1]]
PG_summ$FC_LFQ_Cytosol_Rep2 <- PG_summ[,LFQ[8]]-PG_summ[,LFQ[2]]
PG_summ$FC_LFQ_Cytosol_Rep3 <- PG_summ[,LFQ[9]]-PG_summ[,LFQ[3]]
PG_summ$FC_LFQ_Cytosol_Rep4 <- PG_summ[,LFQ[10]]-PG_summ[,LFQ[4]]
PG_summ$FC_LFQ_Cytosol_Rep5 <- PG_summ[,LFQ[11]]-PG_summ[,LFQ[5]]
PG_summ$FC_LFQ_Cytosol_Rep6 <- PG_summ[,LFQ[12]]-PG_summ[,LFQ[6]]


#Input
PG_summ$FC_LFQ_Input_Rep1 <- PG_summ[,LFQ[19]]-PG_summ[,LFQ[13]]
PG_summ$FC_LFQ_Input_Rep2 <- PG_summ[,LFQ[20]]-PG_summ[,LFQ[14]]
PG_summ$FC_LFQ_Input_Rep3 <- PG_summ[,LFQ[21]]-PG_summ[,LFQ[15]]
PG_summ$FC_LFQ_Input_Rep4 <- PG_summ[,LFQ[22]]-PG_summ[,LFQ[16]]
PG_summ$FC_LFQ_Input_Rep5 <- PG_summ[,LFQ[23]]-PG_summ[,LFQ[17]]
PG_summ$FC_LFQ_Input_Rep6 <- PG_summ[,LFQ[24]]-PG_summ[,LFQ[18]]


#Nucleus
PG_summ$FC_LFQ_Nucleus_Rep1 <- PG_summ[,LFQ[31]]-PG_summ[,LFQ[25]]
PG_summ$FC_LFQ_Nucleus_Rep2 <- PG_summ[,LFQ[32]]-PG_summ[,LFQ[26]]
PG_summ$FC_LFQ_Nucleus_Rep3 <- PG_summ[,LFQ[33]]-PG_summ[,LFQ[27]]
PG_summ$FC_LFQ_Nucleus_Rep4 <- PG_summ[,LFQ[34]]-PG_summ[,LFQ[28]]
PG_summ$FC_LFQ_Nucleus_Rep5 <- PG_summ[,LFQ[35]]-PG_summ[,LFQ[29]]
PG_summ$FC_LFQ_Nucleus_Rep6 <- PG_summ[,LFQ[36]]-PG_summ[,LFQ[30]]


# Filter out cases of not pared maching in a repicate
FC_LFQ <- c("FC_LFQ_Cytosol_Rep1",
            "FC_LFQ_Cytosol_Rep2",
            "FC_LFQ_Cytosol_Rep3",
            "FC_LFQ_Cytosol_Rep4",
            "FC_LFQ_Cytosol_Rep5",
            "FC_LFQ_Cytosol_Rep6",
            
            "FC_LFQ_Input_Rep1",
            "FC_LFQ_Input_Rep2",
            "FC_LFQ_Input_Rep3",
            "FC_LFQ_Input_Rep4",
            "FC_LFQ_Input_Rep5",
            "FC_LFQ_Input_Rep6",
            
            "FC_LFQ_Nucleus_Rep1",
            "FC_LFQ_Nucleus_Rep2",
            "FC_LFQ_Nucleus_Rep3",
            "FC_LFQ_Nucleus_Rep4",
            "FC_LFQ_Nucleus_Rep5",
            "FC_LFQ_Nucleus_Rep6")




#### Filter ####
# Filter for proteins with at least 2 FC values in the experimental groups
PG_summ_Cytosol <- PG_summ[apply(select(PG_summ, FC_LFQ[1:6]), 1, function(x) sum(!is.na(x)) > 1), ]
PG_summ_Cytosol <- subset(PG_summ_Cytosol, select = c(LFQ[1:12], FC_LFQ[1:6],"id"))

PG_summ_Input <- PG_summ[apply(select(PG_summ, FC_LFQ[7:12]), 1, function(x) sum(!is.na(x)) > 1), ]
PG_summ_Input <- subset(PG_summ_Input, select = c(LFQ[13:24], FC_LFQ[7:12],"id"))

PG_summ_Nucleus <- PG_summ[apply(select(PG_summ, FC_LFQ[13:18]), 1, function(x) sum(!is.na(x)) > 1), ]
PG_summ_Nucleus <- subset(PG_summ_Nucleus, select = c(LFQ[25:36], FC_LFQ[13:18],"id"))








#Merge back all filtered values into PG_summ_filtered
PG_summ_filtered <- merge(PG_summ, PG_summ_Cytosol, by = "id", suffixes = c("", "_Cytosol"), all = T)
PG_summ_filtered <- merge(PG_summ_filtered, PG_summ_Input, by = "id", suffixes = c("", "_Input"), all = T)
PG_summ_filtered <- merge(PG_summ_filtered, PG_summ_Nucleus, by = "id", suffixes = c("", "_Nucleus"), all = T)



# Pass values from filtered column to respective colum
PG_summ_filtered[,LFQ] <- PG_summ_filtered[,c(paste0(LFQ[1:12], "_Cytosol"),
                                              paste0(LFQ[13:24], "_Input"),
                                              paste0(LFQ[25:36], "_Nucleus"))]

PG_summ_filtered[,FC_LFQ] <- PG_summ_filtered[,c(paste0(FC_LFQ[1:6], "_Cytosol"),
                                                 paste0(FC_LFQ[7:12], "_Input"),
                                                 paste0(FC_LFQ[13:18], "_Nucleus"))]




#Leave only relevant columns
PG_summ_filtered <- subset(PG_summ_filtered, select = c("Majority.protein.IDs","Gene.names", LFQ, FC_LFQ))


#Remove columns of proteins that didn't pass my selection threshold
PG_summ_filtered <- PG_summ_filtered[apply(select(PG_summ_filtered, LFQ), 1, function(x) sum(!is.na(x)) > 1), ]







##### p.value and FDR ####
#Cytosol
PG_summ_filtered$LFQ.Control_Cytosol_Mean <- rowMeans(PG_summ_filtered[,LFQ[1:6]], na.rm = T)
PG_summ_filtered$LFQ.Patient_Cytosol_Mean <- rowMeans(PG_summ_filtered[,LFQ[7:12]], na.rm = T)

PG_summ_filtered$FC_LFQ_Cytosol_Mean <- rowMeans(PG_summ_filtered[,c("FC_LFQ_Cytosol_Rep1",
                                                                     "FC_LFQ_Cytosol_Rep2",
                                                                     "FC_LFQ_Cytosol_Rep3",
                                                                     "FC_LFQ_Cytosol_Rep4",
                                                                     "FC_LFQ_Cytosol_Rep5",
                                                                     "FC_LFQ_Cytosol_Rep6")], na.rm = T)



PG_summ_filtered$pval_Patient_CytosolXControl <- apply(PG_summ_filtered[,LFQ[1:12]], 1,
                                                       function(x) ifelse(sum(!is.na(x[1:6])) > 1 |
                                                                            sum(!is.na(x[7:12])) > 1,
                                                                          t.test(x[1:6], x[7:12],
                                                                                 paired = T)$p.value, NA))

PG_summ_filtered$minus_LOg10_pval_Patient_CytosolXControl <- -log10(PG_summ_filtered$pval_Patient_CytosolXControl)
PG_summ_filtered$padj_Patient_CytosolXControl <-  p.adjust(PG_summ_filtered$pval_Patient_CytosolXControl, method="BH")
PG_summ_filtered$minus_LOg10_padj_Patient_CytosolXControl <- -log10(PG_summ_filtered$padj_Patient_CytosolXControl)




#Input
PG_summ_filtered$LFQ.Control_Input_Mean <- rowMeans(PG_summ_filtered[,LFQ[13:18]], na.rm = T)
PG_summ_filtered$LFQ.Patient_Input_Mean <- rowMeans(PG_summ_filtered[,LFQ[19:24]], na.rm = T)

PG_summ_filtered$FC_LFQ_Input_Mean <- rowMeans(PG_summ_filtered[,c("FC_LFQ_Input_Rep1",
                                                                   "FC_LFQ_Input_Rep2",
                                                                   "FC_LFQ_Input_Rep3",
                                                                   "FC_LFQ_Input_Rep4",
                                                                   "FC_LFQ_Input_Rep5",
                                                                   "FC_LFQ_Input_Rep6")], na.rm = T)



PG_summ_filtered$pval_Patient_InputXControl <- apply(PG_summ_filtered[,LFQ[13:24]], 1,
                                                     function(x) ifelse(sum(!is.na(x[1:6])) > 1 |
                                                                          sum(!is.na(x[7:12])) > 1,
                                                                        t.test(x[1:6], x[7:12],
                                                                               paired = T)$p.value, NA))

PG_summ_filtered$minus_LOg10_pval_Patient_InputXControl <- -log10(PG_summ_filtered$pval_Patient_InputXControl)
PG_summ_filtered$padj_Patient_InputXControl <-  p.adjust(PG_summ_filtered$pval_Patient_InputXControl, method="BH")
PG_summ_filtered$minus_LOg10_padj_Patient_InputXControl <- -log10(PG_summ_filtered$padj_Patient_InputXControl)



#Nucleus
PG_summ_filtered$LFQ.Control_Nucleus_Mean <- rowMeans(PG_summ_filtered[,LFQ[25:30]], na.rm = T)
PG_summ_filtered$LFQ.Patient_Nucleus_Mean <- rowMeans(PG_summ_filtered[,LFQ[31:36]], na.rm = T)

PG_summ_filtered$FC_LFQ_Nucleus_Mean <- rowMeans(PG_summ_filtered[,c("FC_LFQ_Nucleus_Rep1",
                                                                     "FC_LFQ_Nucleus_Rep2",
                                                                     "FC_LFQ_Nucleus_Rep3",
                                                                     "FC_LFQ_Nucleus_Rep4",
                                                                     "FC_LFQ_Nucleus_Rep5",
                                                                     "FC_LFQ_Nucleus_Rep6")], na.rm = T)



PG_summ_filtered$pval_Patient_NucleusXControl <- apply(PG_summ_filtered[,LFQ[25:36]], 1,
                                                       function(x) ifelse(sum(!is.na(x[1:6])) > 1 |
                                                                            sum(!is.na(x[7:12])) > 1,
                                                                          t.test(x[1:6], x[7:12],
                                                                                 paired = T)$p.value, NA))

PG_summ_filtered$minus_LOg10_pval_Patient_NucleusXControl <- -log10(PG_summ_filtered$pval_Patient_NucleusXControl)
PG_summ_filtered$padj_Patient_NucleusXControl <-  p.adjust(PG_summ_filtered$pval_Patient_NucleusXControl, method="BH")
PG_summ_filtered$minus_LOg10_padj_Patient_NucleusXControl <- -log10(PG_summ_filtered$padj_Patient_NucleusXControl)










## Increased and decreased in each group of samples
PG_summ_filtered$Logic_Cytosol_increase_LFQ = (PG_summ_filtered$padj_Patient_CytosolXControl <= 0.05) & (PG_summ_filtered$FC_LFQ_Cytosol_Mean >= log2(2))

PG_summ_filtered$Logic_Cytosol_decrease_LFQ = (PG_summ_filtered$padj_Patient_CytosolXControl <= 0.05) & (PG_summ_filtered$FC_LFQ_Cytosol_Mean <= -log2(2))



PG_summ_filtered$Logic_Input_increase_LFQ = (PG_summ_filtered$padj_Patient_InputXControl <= 0.05) & (PG_summ_filtered$FC_LFQ_Input_Mean >= log2(2))

PG_summ_filtered$Logic_Input_decrease_LFQ = (PG_summ_filtered$padj_Patient_InputXControl <= 0.05) & (PG_summ_filtered$FC_LFQ_Input_Mean <= -log2(2))



PG_summ_filtered$Logic_Nucleus_increase_LFQ = (PG_summ_filtered$padj_Patient_NucleusXControl <= 0.05) & (PG_summ_filtered$FC_LFQ_Nucleus_Mean >= log2(2))

PG_summ_filtered$Logic_Nucleus_decrease_LFQ = (PG_summ_filtered$padj_Patient_NucleusXControl <= 0.05) & (PG_summ_filtered$FC_LFQ_Nucleus_Mean <= -log2(2))














#### GO annotation ####
#Use Uniprot to add annotations to genes
Uniprot <- fread("K:/Datasets/20160818_Uniprot_HS.tab", sep = "\t", stringsAsFactors = F, data.table = F)
Uniprot_Reviewed <- subset(Uniprot, Status == "reviewed")

PG_summ_filtered <- merge(PG_summ_filtered, Uniprot_Reviewed[,c("Gene names  (primary )", "Gene ontology (molecular function)", "Gene ontology (biological process)", "Gene ontology (cellular component)")], by.x = "Gene.names", by.y = "Gene names  (primary )", all.x = T, all.y = F)

rm(Uniprot, Uniprot_Reviewed)


# Proteins exclusivelly Nucleus
PG_summ_filtered$Nucleus.logical <- ifelse(grepl("GO:0005634", PG_summ_filtered$`Gene ontology (cellular component)`) & 
                                             !(grepl("GO:0005737", PG_summ_filtered$`Gene ontology (cellular component)`)), T, F)


# Proteins exclusivelly Cytosol
PG_summ_filtered$Cytosol.logical <- ifelse(!grepl("GO:0005634", PG_summ_filtered$`Gene ontology (cellular component)`) & 
                                             (grepl("GO:0005737", PG_summ_filtered$`Gene ontology (cellular component)`)), T, F)







fwrite(PG_summ_filtered, file = "PG_summ_filtered.txt", sep = "\t", na = "", quote = F, row.names = F)




#### Heatmap ####

# See: https://stackoverflow.com/questions/21983162/how-to-expand-the-dendogram-in-heatmap-2
foo <- PG_summ_filtered[FC_LFQ]
names(foo) <- c("Cytosol_Rep1",
                "Cytosol_Rep2",
                "Cytosol_Rep3",
                "Cytosol_Rep4",
                "Cytosol_Rep5",
                "Cytosol_Rep6",
                
                "Input_Rep1",
                "Input_Rep2",
                "Input_Rep3",
                "Input_Rep4",
                "Input_Rep5",
                "Input_Rep6",
                
                "Nucleus_Rep1",
                "Nucleus_Rep2",
                "Nucleus_Rep3",
                "Nucleus_Rep4",
                "Nucleus_Rep5",
                "Nucleus_Rep6")


mt <- data.matrix(foo)
rm(foo)
row.names(mt) <- sapply(strsplit(PG_summ_filtered$Gene.names, ";"), "[", 1)

#if scaling values:
mt <- t(scale(t(mt)))
breaks <- seq(min(mt, na.rm = T), max(mt, na.rm = T), length.out = 13)

#without scaling
breaks <- seq(mean(mt, na.rm = T) - 3*sd(mt, na.rm = T), mean(mt, na.rm = T) + 3*sd(mt, na.rm = T), length.out = 13)


#Change all NAs for a super high value then creates a breaking point for my heatmap
mt[is.na(mt)] <- -1000



png(filename = paste0("Heatmap_FC_LFQ.png"), width = 500, height = 1000, pointsize = 15)
heatmap.2(mt,
          dendrogram = c("none"),           #which dendograms to show
          Rowv = T,         #hclustering by row
          Colv = F,         #hclustering by column
          col = c(rgb(0,0,0,.5), rev(brewer.pal(11,"BrBG"))),           #heatmap plotting colors
          breaks=breaks,
          symm=F,symkey=F,symbreaks=T,
          na.color = rgb(0,0,0,0),
          trace = "none",  
          labRow = NA,
          key.xlab = "Fold Change (Log2)\nZ-Score",
          scale = "none",
          cexCol = 1,
          margins = c(9,1),
          keysize = 3,
          lhei=c(1,8), 
          lwid=c(2,5),
          key.title = NA,
          density.info = "none",           #remove super-posed trace in key
          srtCol = 45)             #angle column labels
dev.off()

rm(mt, breaks)









#### Heatmap means ####
foo1 <- subset(PG_summ_filtered, ifelse(apply(PG_summ_filtered[,c("Logic_Cytosol_increase_LFQ",
                                                                  "Logic_Cytosol_decrease_LFQ",
                                                                  "Logic_Nucleus_increase_LFQ",
                                                                  "Logic_Nucleus_decrease_LFQ",
                                                                  "Logic_Input_increase_LFQ",
                                                                  "Logic_Input_decrease_LFQ")], 1,
                                              function(x) sum(x, na.rm = T)) > 0, T, F))


foo <- data.frame("Cytosol" = rowMeans(foo1[FC_LFQ[1:6]], na.rm = T),
                  "Input" = rowMeans(foo1[FC_LFQ[7:12]], na.rm = T),
                  "Nucleus" = rowMeans(foo1[FC_LFQ[13:18]], na.rm = T))



mt <- data.matrix(foo)
row.names(mt) <- sapply(strsplit(foo1$Gene.names, ";"), "[", 1)
rm(foo, foo1)

#if scaling values:
#mt <- t(scale(t(mt)))
#breaks <- seq(min(mt, na.rm = T), max(mt, na.rm = T), length.out = 13)

#without scaling
breaks <- seq(mean(mt, na.rm = T) - 2*sd(mt, na.rm = T), mean(mt, na.rm = T) + 2*sd(mt, na.rm = T), length.out = 13)
breaks <- seq(min(mt, na.rm = T), max(mt, na.rm = T), length.out = 13)

#Change all NAs for a super high value then creates a breaking point for my heatmap
mt[is.na(mt)] <- -10000



png(filename = paste0("Heatmap_FC_LFQ_mean_significants.png"), width = 500, height = 1000, pointsize = 15)
heatmap.2(mt,
          dendrogram = c("none"),           #which dendograms to show
          Rowv = T,         #hclustering by row
          Colv = F,         #hclustering by column
          col = c(rgb(0,0,0,.5), rev(brewer.pal(11,"BrBG"))),           #heatmap plotting colors
          breaks=breaks,
          symm=F, symkey=F, symbreaks=T,
          na.color = rgb(0,0,0,0),
          trace = "none",  
          #   labRow = NA,
          key.xlab = "Fold Change (Log2)\nZ-Score",
          scale = "none",
          cexCol = 2,
          cexRow = .5,
          margins = c(9,5),
          keysize = 3,
          lhei=c(1,8), 
          lwid=c(2,5),
          key.title = NA,
          density.info = "none",           #remove super-posed trace in key
          srtCol = 45)             #angle column labels

dev.off()



rm(mt, breaks)















#### Heatmap mean LFQs ####
PG_summ_filtered$LFQ.intensity.Control_Cytosol_Mean <- rowMeans(PG_summ_filtered[,LFQ[1:6]], na.rm = T)
PG_summ_filtered$LFQ.intensity.Patient_Cytosol_Mean <- rowMeans(PG_summ_filtered[,LFQ[7:12]], na.rm = T)

PG_summ_filtered$LFQ.intensity.Control_Input_Mean <- rowMeans(PG_summ_filtered[,LFQ[13:18]], na.rm = T)
PG_summ_filtered$LFQ.intensity.Patient_Input_Mean <- rowMeans(PG_summ_filtered[,LFQ[19:24]], na.rm = T)

PG_summ_filtered$LFQ.intensity.Control_Nucleus_Mean <- rowMeans(PG_summ_filtered[,LFQ[25:30]], na.rm = T)
PG_summ_filtered$LFQ.intensity.Patient_Nucleus_Mean <- rowMeans(PG_summ_filtered[,LFQ[31:36]], na.rm = T)

foo <- subset(PG_summ_filtered, select = c("LFQ.intensity.Control_Cytosol_Mean",
                                           "LFQ.intensity.Patient_Cytosol_Mean",
                                           "LFQ.intensity.Control_Input_Mean",
                                           "LFQ.intensity.Patient_Input_Mean",
                                           "LFQ.intensity.Control_Nucleus_Mean",
                                           "LFQ.intensity.Patient_Nucleus_Mean"))

names(foo) <- c("Control_Cytosol",
                "Patient_Cytosol",
                "Control_Input",
                "Patient_Input",
                "Control_Nucleus",
                "Patient_Nucleus")


mt <- data.matrix(foo)
row.names(mt) <- sapply(strsplit(PG_summ_filtered$Gene.names, ";"), "[", 1)
rm(foo)

#if scaling values:
mt <- t(scale(t(mt)))
breaks <- seq(min(mt, na.rm = T), max(mt, na.rm = T), length.out = 11)

#Change all NAs for a super high value then creates a breaking point for my heatmap
mt[is.na(mt)] <- -1000



pdf(paste0("Heatmap_LFQ_mean.pdf"), width =5, height = 10)
heatmap.2(mt,
          dendrogram = c("none"),           #which dendograms to show
          Rowv = T,         #hclustering by row
          Colv = T,         #hclustering by column
          col = c(rgb(1,0,0,.5), brewer.pal(9,"Greys")),           #heatmap plotting colors
          breaks=breaks,
          symm=F, symkey=F, symbreaks=T,
          na.color = rgb(0,0,0,0),
          trace = "none",  
          #   labRow = NA,
          key.xlab = "LFQ intensity\nZ-Score",
          scale = "none",
          cexCol = 2,
          cexRow = .5,
          margins = c(12,.1),
          keysize = 3,
          lhei=c(1,8), 
          lwid=c(2,5),
          key.title = NA,
          density.info = "none",           #remove super-posed trace in key
          srtCol = 45)             #angle column labels

dev.off()



rm(mt, breaks)








#### boxplot LFQ ####
#Define plotting dataframe

MyBoxplot <- function(df, variables,
                      AxisName_y = ("Log2"), LabelNames = variables,
                      limits_y_min = round(min(df[variables], na.rm = T))-1,
                      limits_y_max = round(max(df[variables], na.rm = T))+1,
                      limits_breaks = round((limits_y_max+limits_y_min)/15)) {
  
  all <- melt(df,
              id.vars = c("Majority.protein.IDs"), measure.vars=variables)
  
  
  
  bplot <- ggplot(all, aes(variable,value)) +
    
    geom_violin(fill = rgb(0,0,0,.25), na.rm = T) +
    
    geom_boxplot(fill = rgb(0,0,0,0), width = 0.75, na.rm = T, notch = T) +
    
    ylab(AxisName_y) +
    coord_flip(ylim = c(limits_y_min, limits_y_max)) +
    scale_y_continuous(breaks = seq(limits_y_min, limits_y_max, limits_breaks)) +
    scale_x_discrete(labels = LabelNames) +
    
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=30),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=20),
          axis.title.x = element_text(face="bold",
                                      size=25,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=0.4,
                                      size=20))
  
  return(bplot)
  
}




plot_LFQ <- MyBoxplot(PG_summ_filtered, LFQ, AxisName_y = "LFQ (Log2)",
                      LabelNames = c("Control_Cytosol_1",
                                     "Control_Cytosol_2",
                                     "Control_Cytosol_3",
                                     "Control_Cytosol_4",
                                     "Control_Cytosol_5",
                                     "Control_Cytosol_6",
                                     
                                     "Patient_Cytosol_1",
                                     "Patient_Cytosol_2",
                                     "Patient_Cytosol_3",
                                     "Patient_Cytosol_4",
                                     "Patient_Cytosol_5",
                                     "Patient_Cytosol_6",
                                     
                                     "Control_Input_1",
                                     "Control_Input_2",
                                     "Control_Input_3",
                                     "Control_Input_4",
                                     "Control_Input_5",
                                     "Control_Input_6",
                                     
                                     "Patient_Input_1",
                                     "Patient_Input_2",
                                     "Patient_Input_3",
                                     "Patient_Input_4",
                                     "Patient_Input_5",
                                     "Patient_Input_6",
                                     
                                     "Control_Nucleus_1",
                                     "Control_Nucleus_2",
                                     "Control_Nucleus_3",
                                     "Control_Nucleus_4",
                                     "Control_Nucleus_5",
                                     "Control_Nucleus_6",
                                     
                                     "Patient_Nucleus_1",
                                     "Patient_Nucleus_2",
                                     "Patient_Nucleus_3",
                                     "Patient_Nucleus_4",
                                     "Patient_Nucleus_5",
                                     "Patient_Nucleus_6"))


png(paste("Boxplot_LFQ.png", sep = ""), width = 500, height = 500, pointsize = 25)
grid.arrange(plot_LFQ,
             nrow = 1)
dev.off()
rm(plot_LFQ)













#### My plot ####
MyScatterPlot <- function(df, X, Y,
                          xmin = round(min(df[X]))-1, xmax = round(max(df[X]))+1,
                          x_breaks = round((xmax+xmin)/20) + 1,
                          ymin = round(min(df[Y]))-1, ymax = round(max(df[Y]))+1,
                          y_breaks = round((ymax+ymin)/20) + 1,
                          mainTitle = NULL, xTitle = NULL, yTitle = NULL,
                          ShowDensityNuc = F, ShowDensityCyt = F,
                          GradientColorMin = rgb(.5,.5,.5), GradientColorMax = rgb(0,0,0),
                          FilterDF.by = NULL, NegFilterDF.by = NULL,
                          ShowSelected = F, Selected = NULL) {
  
  if (!is.null(FilterDF.by)) {
    if (length(FilterDF.by) == 1) {
      df <- subset(df, df[,FilterDF.by])
    } else {
      df <- subset(df, ifelse(apply(df[,FilterDF.by], 1, function(x) sum(x, na.rm = T)) > 0, T, F))
    }
  }
  
  
  if (!is.null(NegFilterDF.by)) {
    if (length(NegFilterDF.by) == 1) {
      df <- subset(df, !df[,NegFilterDF.by])
    } else {
      df <- subset(df, !ifelse(apply(df[,NegFilterDF.by], 1, function(x) sum(x, na.rm = T)) > 0, T, F))
    }
  }
  
  
  
  plot <- ggplot(data = df, aes(x = eval(parse(text = X)), y = eval(parse(text = Y)))) +
    
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    
    geom_point(na.rm = T,
               size = 3,
               colour = densCols(x = df[,X],
                                 y = df[,Y],
                                 colramp = colorRampPalette(c(GradientColorMin,
                                                              GradientColorMax)))) +
    
    ggtitle(mainTitle) +
    xlab(xTitle) +
    ylab(yTitle) +
    coord_cartesian(xlim = c(xmin, xmax),ylim = c(ymin, ymax)) +
    scale_x_continuous(breaks = seq(xmin, xmax, x_breaks)) +
    scale_y_continuous(breaks = seq(ymin, ymax, y_breaks)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(lineheight=.8,
                                    face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=30),
          axis.title.x = element_text(face="bold",
                                      size=30,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.x  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=20),
          axis.title.y = element_text(face="bold",
                                      size=30,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=.4,
                                      size=20))
  
  
  if (ShowSelected) {
    if (length(Selected) == 1) {
      plot <- plot + geom_point(data = subset(df, df[,Selected]),
                                aes(subset(df, df[,Selected])[,X],
                                    subset(df, df[,Selected])[,Y]),
                                na.rm = T,
                                size = 2,
                                colour = rgb(1,0,0))
    } else {
      plot <- plot + geom_point(data = subset(df, ifelse(apply(df[,Selected], 1,
                                                               function(x) sum(x, na.rm = T)) > 0, T, F)),
                                aes(subset(df, ifelse(apply(df[,Selected], 1,
                                                            function(x) sum(x, na.rm = T)) > 0, T, F))[,X],
                                    subset(df, ifelse(apply(df[,Selected], 1,
                                                            function(x) sum(x, na.rm = T)) > 0, T, F))[,Y]),
                                na.rm = T, size = 2, colour = rgb(1,0,0))
      
    }
    
  }
  
  
  
  if(ShowDensityNuc) {
    
    plot <- plot +
      
      stat_density_2d(data = subset(df, Nucleus.logical),
                      aes(alpha=..level..), colour = rgb(.8,0,0),
                      size = 2, na.rm = T, bins = 5, show.legend = F)
  }
  
  if(ShowDensityCyt) {
    
    plot <- plot +
      
      stat_density_2d(data = subset(df, Cytosol.logical),
                      aes(alpha=..level..), colour = rgb(0,0,.8), 
                      size = 2, na.rm = T, bins = 5, show.legend = F)
  }
  
  return(plot)
}












#### Scatter against Input LFQ ####

plot_Control_Cytosol <- MyScatterPlot(PG_summ_filtered,
                                      X = "LFQ.Control_Cytosol_Mean", Y = "LFQ.Control_Input_Mean",
                                      xmin = 24, xmax = 37,
                                      ymin = 24, ymax = 37,
                                      mainTitle = "Control",
                                      xTitle = "Cytosol (Log2 LFQ)",
                                      yTitle = "Input (Log2 LFQ)",
                                      ShowDensityNuc = T, ShowDensityCyt = T)



plot_Patient_Cytosol <- MyScatterPlot(PG_summ_filtered,
                                      X = "LFQ.Patient_Cytosol_Mean", Y = "LFQ.Patient_Input_Mean",
                                      xmin = 24, xmax = 37,
                                      ymin = 24, ymax = 37,
                                      mainTitle = "Patient",
                                      xTitle = "Cytosol (Log2 LFQ)",
                                      yTitle = "Input (Log2 LFQ)",
                                      ShowDensityNuc = T, ShowDensityCyt = T)



plot_Control_Nucleus <- MyScatterPlot(PG_summ_filtered,
                                      X = "LFQ.Control_Nucleus_Mean", Y = "LFQ.Control_Input_Mean",
                                      xmin = 24, xmax = 37,
                                      ymin = 24, ymax = 37,
                                      mainTitle = "Control",
                                      xTitle = "Nucleus (Log2 LFQ)",
                                      yTitle = "Input (Log2 LFQ)",
                                      ShowDensityNuc = T, ShowDensityCyt = T)



plot_Patient_Nucleus <- MyScatterPlot(PG_summ_filtered,
                                      X = "LFQ.Patient_Nucleus_Mean", Y = "LFQ.Patient_Input_Mean",
                                      xmin = 24, xmax = 37,
                                      ymin = 24, ymax = 37,
                                      mainTitle = "Patient",
                                      xTitle = "Nucleus (Log2 LFQ)",
                                      yTitle = "Input (Log2 LFQ)",
                                      ShowDensityNuc = T, ShowDensityCyt = T)




png(paste("Scatter_againsInput_LFQ.png", sep = ""), width = 1000, height = 1000, pointsize = 25)
grid.arrange(plot_Control_Cytosol, plot_Patient_Cytosol, plot_Control_Nucleus, plot_Patient_Nucleus,
             nrow = 2)
dev.off()
rm(plot_Control_Cytosol, plot_Patient_Cytosol, plot_Control_Nucleus, plot_Patient_Nucleus)

















#### Scatter Nucleus x Cytosol LFQ ####

plot_Control <- MyScatterPlot(PG_summ_filtered,
                              X = "LFQ.Control_Cytosol_Mean", Y = "LFQ.Control_Nucleus_Mean",
                              xmin = 23, xmax = 37,
                              ymin = 23, ymax = 37,
                              mainTitle = "Control",
                              xTitle = "Cytosol (Log2 LFQ)",
                              yTitle = "Nucleus (Log2 LFQ)",
                              ShowDensityNuc = T, ShowDensityCyt = T) + geom_abline(slope = 1)




plot_Patient <- MyScatterPlot(PG_summ_filtered,
                              X = "LFQ.Patient_Cytosol_Mean", Y = "LFQ.Patient_Nucleus_Mean",
                              xmin = 23, xmax = 37,
                              ymin = 23, ymax = 37,
                              mainTitle = "Patient",
                              xTitle = "Cytosol (Log2 LFQ)",
                              yTitle = "Nucleus (Log2 LFQ)",
                              ShowDensityNuc = T, ShowDensityCyt = T) + geom_abline(slope = 1)



png(paste("Scatter_NucleusXCytosol_LFQ.png", sep = ""), width = 1000, height = 500, pointsize = 25)
grid.arrange(plot_Control, plot_Patient,
             nrow = 1)
dev.off()
rm(plot_Control, plot_Patient)












#### Chrstimas Tree LFQ ####

plot_Cytosol <- MyScatterPlot(PG_summ_filtered,
                              X = "FC_LFQ_Cytosol_Mean", Y = "LFQ.Control_Cytosol_Mean",
                              xmin = -6, xmax = 6,
                              ymin = 22, ymax = 38,
                              mainTitle = "Cytosol",
                              xTitle = "Patient/Control (Log2_FC LFQ)",
                              yTitle = "Control (Log2 LFQ)") +
  
  geom_point(data = subset(PG_summ_filtered,
                           ifelse(apply(PG_summ_filtered[,c("Logic_Cytosol_increase_LFQ",
                                                            "Logic_Cytosol_decrease_LFQ")],
                                        1,function(x) sum(x, na.rm = T)) > 0, T, F)),
             na.rm = T,
             size = 4,
             colour = rgb(0,0,1))



plot_Input <- MyScatterPlot(PG_summ_filtered,
                            X = "FC_LFQ_Input_Mean", Y = "LFQ.Control_Input_Mean",
                            xmin = -6, xmax = 6,
                            ymin = 22, ymax = 38,
                            mainTitle = "Input",
                            xTitle = "Patient/Control (Log2_FC LFQ)",
                            yTitle = "Control (Log2 LFQ)") +
  
  geom_point(data = subset(PG_summ_filtered,
                           ifelse(apply(PG_summ_filtered[,c("Logic_Input_increase_LFQ",
                                                            "Logic_Input_decrease_LFQ")],
                                        1, function(x) sum(x, na.rm = T)) > 0, T, F)),
             na.rm = T,
             size = 4,
             colour = rgb(0,0,0))



plot_Nucleus <- MyScatterPlot(PG_summ_filtered,
                              X = "FC_LFQ_Nucleus_Mean", Y = "LFQ.Control_Nucleus_Mean",
                              xmin = -6, xmax = 6,
                              ymin = 22, ymax = 38,
                              mainTitle = "Nucleus",
                              xTitle = "Patient/Control (Log2_FC LFQ)",
                              yTitle = "Control (Log2 LFQ)") +
  
  geom_point(data = subset(PG_summ_filtered,
                           ifelse(apply(PG_summ_filtered[,c("Logic_Nucleus_increase_LFQ",
                                                            "Logic_Nucleus_decrease_LFQ")],
                                        1, function(x) sum(x, na.rm = T)) > 0, T, F)),
             na.rm = T,
             size = 4,
             colour = rgb(1,0,0))




png(paste("Xtree_LFQ.png", sep = ""), width = 1500, height = 500, pointsize = 25)
grid.arrange(plot_Cytosol, plot_Input, plot_Nucleus,
             nrow = 1)
dev.off()
rm(plot_Cytosol, plot_Input, plot_Nucleus)












#### Vulcano Plots ####

foo <- subset(PG_summ_filtered, !Logic_Input_increase_LFQ)
foo <- subset(foo, !Logic_Input_decrease_LFQ)


#Cytosol
plot_Cytosol <- MyScatterPlot(foo,
                              X = "FC_LFQ_Cytosol_Mean",
                              Y = "minus_LOg10_pval_Patient_CytosolXControl",
                              mainTitle = "Cytosol",
                              xTitle = "Patient/Control (Log2_FC LFQ)",
                              yTitle = "p-value (-Log10)", 
                              xmin = -6, xmax = 6,
                              ymin = 0, ymax = 8) +
  
  
  #  stat_function(fun = keilhauer,  xlim = c(-10, 10),
  #                color = rgb(.5,0,.5), size = 2, linetype = "dashed") +
  
  geom_point(data = subset(foo, ifelse(apply(foo[,c("Logic_Cytosol_increase_LFQ",
                                                    "Logic_Cytosol_decrease_LFQ")], 1,
                                             function(x) sum(x, na.rm = T)) > 0, T, F)),
             na.rm = T,
             size = 4,
             colour = rgb(0,0,1))# + 

#  geom_text_repel(data = subset(foo,
#                                Logic_Cytosol_increase_LFQ |
#                                  Logic_Cytosol_decrease_LFQ),
#                  aes(x = FC_LFQ_Cytosol_Mean, y = minus_LOg10_pval_Patient_CytosolXControl,
#                      label = Gene.names),
#                  size = 3,
#                  colour = rgb(0,0,1),
#                  fontface = 2)


#Input
plot_Input <- MyScatterPlot(foo,
                            X = "FC_LFQ_Input_Mean",
                            Y = "minus_LOg10_pval_Patient_InputXControl",
                            mainTitle = "Input",
                            xTitle = "Patient/Control (Log2_FC LFQ)",
                            yTitle = "p-value (-Log10)", 
                            xmin = -6, xmax = 6,
                            ymin = 0, ymax = 8) +
  
  
  #  stat_function(fun = keilhauer, xlim = c(-10, 10),
  #                color = rgb(.5,0,.5), size = 2, linetype = "dashed") +
  
  geom_point(data = subset(foo, ifelse(apply(foo[,c("Logic_Input_increase_LFQ",
                                                    "Logic_Input_decrease_LFQ")], 1,
                                             function(x) sum(x, na.rm = T)) > 0, T, F)),
             na.rm = T,
             size = 4,
             colour = rgb(0,0,0)) #+ 

#  geom_text_repel(data = subset(foo,
#                                Logic_Input_increase_LFQ |
#                                  Logic_Input_decrease_LFQ),
#                  aes(x = FC_LFQ_Input_Mean, y = minus_LOg10_pval_Patient_InputXControl,
#                      label = Gene.names),
#                  size = 3,
#                  colour = rgb(0,0,0),
#                  fontface = 2)



#Nucleus
plot_Nucleus <- MyScatterPlot(foo,
                              X = "FC_LFQ_Nucleus_Mean",
                              Y = "minus_LOg10_pval_Patient_NucleusXControl",
                              mainTitle = "Nucleus",
                              xTitle = "Patient/Control (Log2_FC LFQ)",
                              yTitle = "p-value (-Log10)", 
                              xmin = -6, xmax = 6,
                              ymin = 0, ymax = 8) +
  
  
  #  stat_function(fun = keilhauer,  xlim = c(-10, 10),
  #                color = rgb(.5,0,.5), size = 2, linetype = "dashed") +
  
  geom_point(data = subset(foo, ifelse(apply(foo[,c("Logic_Nucleus_increase_LFQ",
                                                    "Logic_Nucleus_decrease_LFQ")], 1,
                                             function(x) sum(x, na.rm = T)) > 0, T, F)),
             na.rm = T,
             size = 4,
             colour = rgb(1,0,0))# + 

#  geom_text_repel(data = subset(foo,
#                                Logic_Nucleus_increase_LFQ |
#                                  Logic_Nucleus_decrease_LFQ),
#                  aes(x = FC_LFQ_Nucleus_Mean, y = minus_LOg10_pval_Patient_NucleusXControl,
#                      label = Gene.names),
#                 size = 3,
#                  colour = rgb(1,0,0),
#                  fontface = 2)



png(paste("Vulcano_LFQ_no_Input.png", sep = ""), width = 1500, height = 500, pointsize = 25)
grid.arrange(plot_Cytosol, plot_Input, plot_Nucleus,
             nrow = 1)
dev.off()
rm(plot_Cytosol, plot_Input, plot_Nucleus)












#### List proteins Ethiraj ####
importins <- c("IPO4",
               "IPO5",
               "IPO7",
               "IPO8",
               "IPO9",
               "IPO11",
               "IPO13",
               "KPNA1",
               "KPNA2",
               "KPNA3",
               "KPNA4",
               "KPNA5",
               "KPNA6",
               "KPNA7",
               "KPNB1",
               "TNPO1",
               "TNPO2",
               "TNPO3")

exportins <- c("XPO1",
               "CSE1L",
               "XPOT",
               "XPO4",
               "XPO5",
               "XPO6",
               "XPO7")


regulators_cyto <- c("ARPC4",
                     "PAK4",
                     "CCNA1",
                     "ACTR2",
                     "PPP3CA",
                     "SSH2",
                     "IQGAP1",
                     "CDC42EP2",
                     "CDC42EP3",
                     "CFL1",
                     "CYFIP1",
                     "ACTR3",
                     "MID1",
                     "RAC1",
                     "CRK",
                     "CDC42",
                     "MAPRE1",
                     "ARHGAP6",
                     "AURKB",
                     "CLASP1",
                     "PAK1",
                     "CALD1",
                     "CDC42BPA",
                     "VASP",
                     "LIMK1",
                     "ARPC1B",
                     "MACF1",
                     "PHLDB2",
                     "NCK1",
                     "EZR",
                     "CDK5R1",
                     "MYLK",
                     "ARAP1",
                     "PPP3CB",
                     "CYFIP2",
                     "ARHGDIB",
                     "ARHGEF11",
                     "LIMK2",
                     "MAP3K11",
                     "ARPC5",
                     "CASK",
                     "RACGAP1",
                     "MARK2",
                     "MAP4",
                     "ARPC3",
                     "LLGL1",
                     "DIAPH1",
                     "TIAM1",
                     "GSN",
                     "SSH1",
                     "ARPC2",
                     "CLIP1",
                     "NCK2",
                     "MAPRE2",
                     "WAS",
                     "ARFIP2",
                     "WASF1",
                     "MYLK2",
                     "MSN",
                     "MAPT",
                     "CLASP2",
                     "FSCN2",
                     "CIT",
                     "BAIAP2",
                     "WASL",
                     "FNBP1L",
                     "MAPK13",
                     "PPP1R12A",
                     "PFN2",
                     "AURKC",
                     "CDK5",
                     "PPP1R12B",
                     "CLIP2",
                     "ROCK1",
                     "PIKFYVE",
                     "DSTN",
                     "CCNB2",
                     "AURKA",
                     "STMN1",
                     "CALM1",
                     "RHOA",
                     "CTTN",
                     "IQGAP2",
                     "RDX")




#### Cumulative plots #### 

#Mean
MyCumPlot <- function(df, X,
                      xlim.min = -2, xlim.max = 2, xlim.breaks = .5,
                      ylim.min = -1, ylim.max = 1, ylim.breaks = .1,
                      Title = NULL, xlabel = "Patient/Control (Log2_FC LFQ)", ylabel = "Cumulative density",
                      Plot.main = T, main.size = 4, main.color = rgb(0,0,0,1),
                      ShowSelected = T, selected.size = 2.5) {
  
  plot <- ggplot(df, 
                 aes(x = df[,X])) +
    
    geom_vline(xintercept = 0,
               linetype = 2,
               color = rgb(.25,.25,.25)) + 
    
    geom_hline(yintercept = .5,
               linetype = 2,
               color = rgb(.25,.25,.25)) + 
    
    coord_cartesian(xlim = c(xlim.min, xlim.max)) +
    scale_x_continuous(breaks=seq(xlim.min, xlim.max, xlim.breaks)) +
    scale_y_continuous(breaks=seq(ylim.min, ylim.max, ylim.breaks)) +
    ylab("Cumulative density") +
    xlab(xlabel) +
    labs(title = Title) +
    theme_bw() +
    theme(plot.title = element_text(face="bold",
                                    size = 35,
                                    hjust = 0.5,
                                    vjust = 0.4),
          axis.title.x = element_text(face="bold",
                                      size=30,
                                      hjust = 0.5,
                                      vjust = 0.4),
          axis.text.x  = element_text(face = "bold", color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=25),
          axis.title.y = element_text(face="bold",
                                      size=30,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(face = "bold", color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=25),
          panel.grid=element_blank())
  
  if (Plot.main) {
    plot <- plot +
      geom_step(stat="ecdf",
                size = main.size, 
                color = main.color,
                na.rm = T)
  }
  
  
  if (ShowSelected) {
    
    plot <- plot +
      
      
      #### Exportins 
      geom_step(data = subset(df, df[,"Gene.names"] %in% exportins),
                aes(x = subset(df, df[,"Gene.names"] %in% exportins)[,X]),
                stat="ecdf",
                size = selected.size, 
                color = rgb(.75,0,.75,1),
                na.rm = T) +
      
      annotate("text",
               x = -2,
               y = 1,
               hjust = 0,
               label = paste0("Exportins\n", "p-value <= ",
                              round(wilcox.test(df[,X],
                                                subset(df, df[,"Gene.names"] %in% exportins)[,X])$p.value,
                                    digits = 7)),
               fontface = "bold",
               colour = rgb(.75,0,.75,1),
               size = 6) +
      
      
      #### Importins 
      geom_step(data = subset(df, df[,"Gene.names"] %in% importins),
                aes(x = subset(df, df[,"Gene.names"] %in% importins)[,X]),
                stat="ecdf",
                size = selected.size,  
                color = rgb(0,.75,.75,1),
                na.rm = T) +
      
      annotate("text",
               x = -2,
               y = .8,
               hjust = 0,
               label = paste0("Importins\n", "p-value <= ",
                              round(wilcox.test(df[,X],
                                                subset(df, df[,"Gene.names"] %in% importins)[,X])$p.value,
                                    digits = 7)),
               fontface = "bold",
               colour = rgb(0,.75,.75,1),
               size = 6) +
      
      
      #### Cyto Regulators
      geom_step(data = subset(df, df[,"Gene.names"] %in% regulators_cyto),
                aes(x = subset(df, df[,"Gene.names"] %in% regulators_cyto)[,X]),
                stat="ecdf",
                size = selected.size,  
                color = rgb(.75,.75,0,1),
                na.rm = T) +
      
      annotate("text",
               x = -2,
               y = .6,
               hjust = 0,
               label = paste0("Cyto Regulators\n", "p-value <= ",
                              round(wilcox.test(df[,X],
                                                subset(df, df[,"Gene.names"] %in% regulators_cyto)[,X])$p.value,
                                    digits = 7)),
               fontface = "bold",
               colour = rgb(.75,.75,0,1),
               size = 6)
  }
}




Cytosol <- MyCumPlot(df = PG_summ_filtered, X = "FC_LFQ_Cytosol_Mean", Title = "Cytosol") 
Input <- MyCumPlot(df = PG_summ_filtered, X = "FC_LFQ_Input_Mean", Title = "Input") 
Nucleus <- MyCumPlot(df = PG_summ_filtered, X = "FC_LFQ_Nucleus_Mean", Title = "Nucleus") 




Vulcano_Cytosol <- MyScatterPlot(PG_summ_filtered,
                              X = "FC_LFQ_Cytosol_Mean",
                              Y = "minus_LOg10_pval_Patient_CytosolXControl",
                              mainTitle = "Cytosol",
                              xTitle = "Patient/Control (Log2_FC LFQ)",
                              yTitle = "p-value (-Log10)", 
                              xmin = -2, xmax = 2,
                              ymin = 0, ymax = 5) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% exportins),
             na.rm = T,
             size = 4,
             colour = rgb(.75,0,.75,.75)) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% importins),
             na.rm = T,
             size = 4,
             colour = rgb(0,.75,.75,.75)) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% regulators_cyto),
             na.rm = T,
             size = 4,
             colour = rgb(.75,.75,0,.75))






Vulcano_Input <- MyScatterPlot(PG_summ_filtered,
                                 X = "FC_LFQ_Input_Mean",
                                 Y = "minus_LOg10_pval_Patient_InputXControl",
                                 mainTitle = "Input",
                                 xTitle = "Patient/Control (Log2_FC LFQ)",
                                 yTitle = "p-value (-Log10)", 
                                 xmin = -2, xmax = 2,
                                 ymin = 0, ymax = 5) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% exportins),
             na.rm = T,
             size = 4,
             colour = rgb(.75,0,.75,.75)) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% importins),
             na.rm = T,
             size = 4,
             colour = rgb(0,.75,.75,.75)) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% regulators_cyto),
             na.rm = T,
             size = 4,
             colour = rgb(.75,.75,0,.75))




Vulcano_Nucleus <- MyScatterPlot(PG_summ_filtered,
                                 X = "FC_LFQ_Nucleus_Mean",
                                 Y = "minus_LOg10_pval_Patient_NucleusXControl",
                                 mainTitle = "Nucleus",
                                 xTitle = "Patient/Control (Log2_FC LFQ)",
                                 yTitle = "p-value (-Log10)", 
                                 xmin = -2, xmax = 2,
                                 ymin = 0, ymax = 5) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% exportins),
             na.rm = T,
             size = 4,
             colour = rgb(.75,0,.75,.75)) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% importins),
             na.rm = T,
             size = 4,
             colour = rgb(0,.75,.75,.75)) +
  
  geom_point(data = subset(PG_summ_filtered, PG_summ_filtered[,"Gene.names"] %in% regulators_cyto),
             na.rm = T,
             size = 4,
             colour = rgb(.75,.75,0,.75))





png("CumulativePlot.png", width = 1500, height = 1000, pointsize = 25)
grid.arrange(Vulcano_Cytosol, Vulcano_Input, Vulcano_Nucleus,
             Cytosol, Input, Nucleus, ncol=3)
dev.off()
rm(Cytosol, Input, Nucleus)








#### Writing table for sharing ####

Table_to_write <- subset(PG_summ_filtered, select = c("Majority.protein.IDs",
                                                      "Gene.names",
                                                      LFQ,
                                                      "LFQ.Control_Cytosol_Mean",
                                                      "LFQ.Patient_Cytosol_Mean",
                                                      "FC_LFQ_Cytosol_Mean",
                                                      "pval_Patient_CytosolXControl",
                                                      "padj_Patient_CytosolXControl",
                                                      "Logic_Cytosol_increase_LFQ",
                                                      "Logic_Cytosol_decrease_LFQ",
                                                      
                                                      "LFQ.Control_Nucleus_Mean",
                                                      "LFQ.Patient_Nucleus_Mean",
                                                      "FC_LFQ_Nucleus_Mean",
                                                      "pval_Patient_NucleusXControl",
                                                      "padj_Patient_NucleusXControl",
                                                      "Logic_Nucleus_increase_LFQ",
                                                      "Logic_Nucleus_decrease_LFQ",
                                                      
                                                      "LFQ.Control_Input_Mean",
                                                      "LFQ.Patient_Input_Mean",
                                                      "FC_LFQ_Input_Mean",
                                                      "pval_Patient_InputXControl",
                                                      "padj_Patient_InputXControl",
                                                      "Logic_Input_increase_LFQ",
                                                      "Logic_Input_decrease_LFQ"))

names(Table_to_write) <- c("Uniprot.IDs",
                           "Gene.names",
                           LFQ,
                           "Mean_LFQ.Control_Cytosol",
                           "Mean_LFQ.Patient_Cytosol",
                           "FC_LFQ_Cytosol_Mean_(Log2)",
                           "pval_Patient_CytosolXControl",
                           "FDR_Patient_CytosolXControl",
                           "Logic_Cytosol_increase",
                           "Logic_Cytosol_decrease",
                           
                           "Mean_LFQ.Control_Nucleus",
                           "Mean_LFQ.Patient_Nucleus",
                           "FC_LFQ_Nucleus_Mean_(Log2)",
                           "pval_Patient_NucleusXControl",
                           "FDR_Patient_NucleusXControl",
                           "Logic_Nucleus_increase",
                           "Logic_Nucleus_decrease",
                           
                           "Mean_LFQ.Control_Input",
                           "Mean_LFQ.Patient_Input",
                           "FC_LFQ_Input_Mean_(Log2)",
                           "pval_Patient_InputXControl",
                           "FDR_Patient_InputXControl",
                           "Logic_Input_increase",
                           "Logic_Input_decrease")



fwrite(Table_to_write, file = "ProteinGroups_analyze_CV.txt", sep = "\t", na = "", quote = F, row.names = F)
rm(Table_to_write)


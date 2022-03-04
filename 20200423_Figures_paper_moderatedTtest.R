setwd("K:/Collaborations/Ethiraj_Ravindran/20190410_Ethiraj_NucFrac_Rep1-6/output")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(limma)
library(ggrepel)
library(gplots)
library(data.table)







#### proteinGroups table load and preparation ####
PG_summ_Cytosol <- fread("PG_summ_Cytosol.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
PG_summ_Input <- fread("PG_summ_Input.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)
PG_summ_Nucleus <- fread("PG_summ_Nucleus.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)



### use Dplyr to remove duplicated rows. For some reason HLA-A rows are duplicated...
PG_summ_Cytosol <- PG_summ_Cytosol %>% distinct(Gene.names, Majority.protein.IDs, .keep_all = TRUE)
PG_summ_Input <- PG_summ_Input %>% distinct(Gene.names, Majority.protein.IDs, .keep_all = TRUE)
PG_summ_Nucleus <- PG_summ_Nucleus %>% distinct(Gene.names, Majority.protein.IDs, .keep_all = TRUE)


#### My plot ####
MyScatterPlot <- function(df, X, Y,
                          xmin = round(min(df[X]))-1, xmax = round(max(df[X]))+1,
                          ymin = round(min(df[Y]))-1, ymax = round(max(df[Y]))+1,
                          x_breaks = round((xmax+xmin)/20) + 1,
                          y_breaks = round((ymax+ymin)/20) + 1,
                          mainTitle = NULL, xTitle = NULL, yTitle = NULL,
                          ShowDensityNuc = F, ShowDensityCyt = F,
                          GradientColorMin = rgb(.5,.5,.5), GradientColorMax = rgb(0,0,0),
                          FilterDF.by = NULL, NegFilterDF.by = NULL,
                          ShowSelected = F, Selected = NULL,
                          PointSize = 2, PointAlpha = 0.5) {
  
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
               size = PointSize,
               alpha = PointAlpha,
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
                                    #face="bold", 
                                    vjust=0.5,
                                    hjust = 0.5,
                                    size=15),
          axis.title.x = element_text(#face="bold",
                                      size=15,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.x  = element_text(#face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=1,
                                      size=10),
          axis.title.y = element_text(#face="bold",
                                      size=15,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(#face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=.4,
                                      size=10),
          aspect.ratio = 1)
  
  
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



#### Fold changes comparison ####
PG_merge_Cytosol <- merge(PG_summ_Cytosol[, c("Gene.names", "Majority.protein.IDs", "id",
                                              "logFC", "Logic_increase","Logic_decrease")],
                          PG_summ_Input[, c("id", "logFC", "Logic_increase","Logic_decrease")], 
                          by = "id", suffixes = c("", "_Input"))
PG_merge_Cytosol <- subset(PG_merge_Cytosol, Logic_increase | Logic_decrease)

PG_merge_Nucleus <- merge(PG_summ_Nucleus[, c("Gene.names", "Majority.protein.IDs", "id",
                                              "logFC", "Logic_increase","Logic_decrease")],
                          PG_summ_Input[, c("id", "logFC", "Logic_increase","Logic_decrease")], 
                          by = "id", suffixes = c("", "_Input"))
PG_merge_Nucleus <- subset(PG_merge_Nucleus, Logic_increase | Logic_decrease)






# Cytosol
plot_Cytosol <- MyScatterPlot(PG_merge_Cytosol,
                              X = "logFC",
                              Y = "logFC_Input",
                              #mainTitle = "Cytosol",
                              xTitle = "Cytosol,\nPatient/Control (Log2_FC LFQ)",
                              yTitle = "Whole proteome,\nPatient/Control (Log2_FC LFQ)", 
                              xmin = -6, xmax = 6, x_breaks = 2,
                              ymin = -6, ymax = 6, y_breaks = 2,
                              GradientColorMin = rgb(0,0,1),
                              GradientColorMax = rgb(0,0,1),
                              PointSize = 3, PointAlpha = .75)


# Nucleus
plot_Nucleus <- MyScatterPlot(PG_merge_Nucleus,
                              X = "logFC",
                              Y = "logFC_Input",
                              #mainTitle = "Cytosol",
                              xTitle = "Nucleus,\nPatient/Control (Log2_FC LFQ)",
                              yTitle = "Whole proteome,\nPatient/Control (Log2_FC LFQ)", 
                              xmin = -6, xmax = 6, x_breaks = 2,
                              ymin = -6, ymax = 6, y_breaks = 2,
                              GradientColorMin = rgb(0,0,1),
                              GradientColorMax = rgb(0,0,1),
                              PointSize = 3, PointAlpha = .75)



pdf(paste("FC_against_input.pdf", sep = ""), width = 10, height = 5, useDingbats = F)
grid.arrange(plot_Cytosol, plot_Nucleus,
             nrow = 1)
dev.off()
rm(plot_Cytosol, plot_Nucleus)





#### Vulcano Plots ####

Input_increase <- subset(PG_summ_Input, Logic_increase)$Majority.protein.IDs
Input_decrease <- subset(PG_summ_Input, Logic_decrease)$Majority.protein.IDs



#Cytosol
PG_summ_Cytosol_subset <- subset(PG_summ_Cytosol, !((Majority.protein.IDs %in%
                                                       Input_increase) & Logic_increase) &
                                   !((Majority.protein.IDs %in%
                                        Input_decrease) & Logic_decrease))

plot_Cytosol <- MyScatterPlot(PG_summ_Cytosol_subset,
                              X = "logFC",
                              Y = "logadj.P.Val",
                              mainTitle = "Cytosol",
                              xTitle = "Patient/Control (Log2_FC LFQ)",
                              yTitle = "p-value (-Log10)", 
                              xmin = -6, xmax = 6,
                              ymin = 0, ymax = 6) +
  
  geom_point(data = subset(PG_summ_Cytosol_subset, Logic_increase | Logic_decrease),
             na.rm = T, shape = 21, size = 3, alpha = 1, colour = rgb(0,0,1))


#Input
plot_Input <- MyScatterPlot(PG_summ_Input,
                            X = "logFC",
                            Y = "logadj.P.Val",
                            mainTitle = "Input",
                            xTitle = "Patient/Control (Log2_FC LFQ)",
                            yTitle = "p-value (-Log10)", 
                            xmin = -6, xmax = 6,
                            ymin = 0, ymax = 6) +
  
  geom_point(data = subset(PG_summ_Input, Logic_increase | Logic_decrease),
             na.rm = T, shape = 21, size = 3, alpha = 1, colour = rgb(0,0,1))

#Nucleus
PG_summ_Nucleus_subset <- subset(PG_summ_Nucleus, !((Majority.protein.IDs %in%
                                                       Input_increase) & Logic_increase) &
                                   !((Majority.protein.IDs %in%
                                        Input_decrease) & Logic_decrease))

plot_Nucleus <- MyScatterPlot(PG_summ_Nucleus_subset,
                              X = "logFC",
                              Y = "logadj.P.Val",
                              mainTitle = "Nucleus",
                              xTitle = "Patient/Control (Log2_FC LFQ)",
                              yTitle = "p-value (-Log10)", 
                              xmin = -6, xmax = 6,
                              ymin = 0, ymax = 6) +
  
  geom_point(data = subset(PG_summ_Nucleus_subset, Logic_increase | Logic_decrease),
             na.rm = T, shape = 21, size = 3, alpha = 1, colour = rgb(0,0,1))



pdf(paste("Vulcano_LFQ_no_Input.pdf", sep = ""), width = 15, height = 5, useDingbats = F)
grid.arrange(plot_Cytosol, plot_Input, plot_Nucleus,
             nrow = 1)
dev.off()
rm(plot_Cytosol, plot_Input, plot_Nucleus)












#### List proteins Ethiraj ####
{importins <- c("IPO4",
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
}



#### Cumulative plots #### 

#Mean
MyCumPlot <- function(df, X,
                      xlim.min = -2, xlim.max = 2, xlim.breaks = .5,
                      ylim.min = -0, ylim.max = 1, ylim.breaks = .2,
                      Title = NULL, xlabel = "Patient/Control (Log2_FC LFQ)", ylabel = "Cumulative density",
                      Plot.main = T, main.size = 1, main.color = rgb(0,0,0,1),
                      CytoRegulators = T, selected.size = 1,
                      Importins = T, Exportins = T) {
  
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
    theme(plot.title = element_text(#face="bold",
                                    size = 15,
                                    hjust = 0.5,
                                    vjust = 0.4),
          axis.title.x = element_text(#face="bold",
                                      size=15,
                                      hjust = 0.5,
                                      vjust = 0.4),
          axis.text.x  = element_text(#face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=10),
          axis.title.y = element_text(#face="bold",
                                      size=15,
                                      hjust = 0.5,
                                      vjust = 1.5),
          axis.text.y  = element_text(#face = "bold",
                                      color = "black",
                                      angle=0, 
                                      vjust=0.5,
                                      hjust = 0.5,
                                      size=10),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1)
  
  if (Plot.main) {
    plot <- plot +
      geom_step(stat="ecdf",
                size = main.size, 
                color = main.color,
                na.rm = T)
  }
  
  
  if (Exportins) {
    
    plot <- plot +
      geom_step(data = subset(df, df[,"Gene.names"] %in% exportins),
                aes(x = subset(df, df[,"Gene.names"] %in% exportins)[,X]),
                stat="ecdf",
                size = selected.size, 
                color = rgb(.75,0,.75,1),
                na.rm = T) +
      
      annotate("text",
               x = xlim.min,
               y = 1,
               hjust = 0,
               label = paste0("Exportins\n", "p-value <= ",
                              round(wilcox.test(df[,X],
                                                subset(df, df[,"Gene.names"] %in% exportins)[,X])$p.value,
                                    digits = 3)),
               colour = rgb(.75,0,.75,1),
               size = 6)
  }
  
  
  if (Importins) {
    
    plot <- plot +
      geom_step(data = subset(df, df[,"Gene.names"] %in% importins),
                aes(x = subset(df, df[,"Gene.names"] %in% importins)[,X]),
                stat="ecdf",
                size = selected.size,  
                color = rgb(0,.75,.75,1),
                na.rm = T) +
      
      annotate("text",
               x = xlim.min,
               y = .8,
               hjust = 0,
               label = paste0("Importins\n", "p-value <= ",
                              round(wilcox.test(df[,X],
                                                subset(df, df[,"Gene.names"] %in% importins)[,X])$p.value,
                                    digits = 3)),
               colour = rgb(0,.75,.75,1),
               size = 6)
    
  }
  
  
  if (CytoRegulators) {
    
    plot <- plot +
      geom_step(data = subset(df, df[,"Gene.names"] %in% regulators_cyto),
                aes(x = subset(df, df[,"Gene.names"] %in% regulators_cyto)[,X]),
                stat="ecdf",
                size = selected.size,  
                color = rgb(.75,.75,0,1),
                na.rm = T) +
      
      annotate("text",
               x = xlim.min,
               y = .6,
               hjust = 0,
               label = paste0("Cyto Regulators\n", "p-value <= ",
                              round(wilcox.test(df[,X],
                                                subset(df, df[,"Gene.names"] %in% regulators_cyto)[,X])$p.value,
                                    digits = 3)),
               colour = rgb(.75,.75,0,1),
               size = 6)
  }
  return(plot)
}




Cytosol <- MyCumPlot(df = PG_summ_Cytosol, X = "logFC", Title = "Cytosol") 
Input <- MyCumPlot(df = PG_summ_Input, X = "logFC", Title = "Input") 
Nucleus <- MyCumPlot(df = PG_summ_Nucleus, X = "logFC", Title = "Nucleus") 







MyVulcanoPlot <- function(df,
                          X = "logFC",
                          Y = "logadj.P.Val",
                          xmin = -6, xmax = 6,
                          ymin = 0, ymax = 6,
                          x_breaks = 2,
                          y_breaks = 1,
                          mainTitle = NULL, 
                          xTitle = "IV.3/Control (Log2_FC LFQ)",
                          yTitle = "p-value (-Log10)",
                          Exp = F, Imp = F, Cyto = F, NUP85=T, DiffExp = T) {
  
  plot <- MyScatterPlot(df=df, X=X, Y=Y,
                        mainTitle = mainTitle, xTitle = xTitle, yTitle = yTitle,
                        xmin = xmin, xmax = xmax, x_breaks = x_breaks,
                        ymin = ymin, ymax = ymax, y_breaks = y_breaks,
                        ShowDensityNuc = F, ShowDensityCyt = F,
                        GradientColorMin = rgb(.5,.5,.5), GradientColorMax = rgb(0,0,0),
                        FilterDF.by = NULL, NegFilterDF.by = NULL,
                        ShowSelected = F, Selected = NULL, PointSize = 3, PointAlpha = .5) 
  
  if(DiffExp) {
    plot <- plot + 
      geom_point(data = subset(df, Logic_increase | Logic_decrease),
                 na.rm = T, size = 3, alpha = 0.75, colour = rgb(0,0,.75))
  }
  
  if(Exp) {
    plot <- plot + 
      geom_point(data = subset(df, df[,"Gene.names"] %in% exportins), #shape = 21,
                 na.rm = T, size = 3, alpha = 0.75, colour = rgb(.75,0,.75))
  }
  
  if(Imp) {
    plot <- plot + 
      geom_point(data = subset(df, df[,"Gene.names"] %in% importins), #shape = 21,
                 na.rm = T, size = 3, alpha = 0.75, colour = rgb(0,.75,.75))
  }
  
  if(Cyto) {
    plot <- plot + 
      geom_point(data = subset(df, df[,"Gene.names"] %in% regulators_cyto), #shape = 21,
                 na.rm = T, size = 3, colour = rgb(.75,.75,0))
  }
  
  if(NUP85) {
    plot <- plot + 
      geom_point(data = subset(df, df[,"Gene.names"] == "NUP85"), #shape = 21,
                 na.rm = T, size = 3, colour = rgb(0,.75,0)) +
      geom_text_repel(data = subset(df, df[,"Gene.names"] == "NUP85"),
                      aes(label=Gene.names),
                      color = rgb(0,.75,0),
                      size = 6,
                      fontface = "bold")
  }
  
  return(plot)
}

Vulcano_Cytosol <- MyVulcanoPlot(df = PG_summ_Cytosol,
                                 mainTitle = "Cytosol")

Vulcano_Input <- MyVulcanoPlot(PG_summ_Input,
                               mainTitle = "Input")

Vulcano_Nucleus <- MyVulcanoPlot(PG_summ_Nucleus,
                                 mainTitle = "Nucleus")



pdf("AlltogetherPlot.pdf", width = 15, height = 10, useDingbats = F)
grid.arrange(Vulcano_Cytosol, Vulcano_Input, Vulcano_Nucleus,
             Cytosol, Input, Nucleus, ncol=3)
dev.off()
rm(Cytosol, Input, Nucleus,
   Vulcano_Cytosol, Vulcano_Input, Vulcano_Nucleus)

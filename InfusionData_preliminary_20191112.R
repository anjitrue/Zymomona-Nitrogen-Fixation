library(ggplot2)
library(plotly)

#### WT_3c_FAIMS_ISPHpeptides_infusion_20191112 ####

WT_FAIMS_ISPH <- read.csv("G:/Projects/Proteomics/Zymomona/WT_3c_FAIMS_ISPHpeptides_infusion_20191112.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

WT_FAIMS_ISPH$Peptide <- as.factor(WT_FAIMS_ISPH$Peptide)

p <- ggplot(WT_FAIMS_ISPH, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
  geom_line(aes(color=Peptide), size=2) + 
  scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3'))  
p




#### ISPH_3a_FAIMS_ISPHpeptides_infusion_20191112 ####

ISPH_FAIMS_ISPH <- read.csv("G:/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_ISPHpeptides_infusion_20191112.csv", 
                          header = TRUE, sep = ",", stringsAsFactors = FALSE)


ISPH_plot <- function(x){
  x[,1] <- as.factor(x[,1])
  p <- ggplot(x, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
    geom_line(aes(color=Peptide), size=2) + 
    ylim(0,1.0e+06) +
    xlim(0.9, 2.1)+
    scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3'))+  
    theme_light()
  return(p)
}

ISPH_plot(ISPH_FAIMS_ISPH)

#### ISPH_3a_FAIMS_ISPHpeptides_infusion_240K ####

ISPH_FAIMS_240K <- read.csv("G:/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_ISPHpeptides_infusion_240K.csv", 
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_plot(ISPH_FAIMS_240K)

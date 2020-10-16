library(tidyverse)
library(R.matlab)
library(cowplot)
library(reshape2)
setwd('D:/我的文章/normalized_SIF_GPP/code/')

# import data
all_datas <- readMat('results/CVs_wheat.mat')

all_data <- as.data.frame(all_datas['CVs.far'])  
colnames(all_data) <- c('Raw-Data','Nadir-NIRv', 'Nadir-KD', 'Total-FPAR', expression(Total-i[0]))
all_data <- melt(all_data)

plot_wheat_far <- ggplot(all_data, mapping = aes(x = variable, y = value, fill = variable)) + ylab(expression(CV)) +  ggtitle('Wheat - Far-red SIF') + 
  #geom_violin(aes(fill = variable),trim=False, size = 2.5, show.legend = FALSE, adjust = .8) + 
  stat_boxplot(geom ='errorbar', width = 0.6, size = 1, show.legend = FALSE) +
  geom_boxplot(notch=T, size = 0.8, show.legend = FALSE, width = 0.6, outlier.shape = NA) + 
  geom_jitter( stroke = 0, alpha = .2, fill = 'black', size = 2, width = 0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", color = 'white', alpha = 1, size=2, show.legend = F) +  
  scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels = c('Raw-Data','Nadir-NIRv', 'Nadir-KD', 'Total-FPAR', expression(Total-i[0]))) + 
  #scale_color_brewer(palette="Dark2") + 
  theme_classic() + 
  theme(plot.title = element_text(face="bold", color="black",size=15, angle=0),
    axis.title.x=element_blank(),
        axis.text.x = element_text(color="black",size=10, angle=0),
        axis.text.y = element_text( color="black", size=10, angle=0),
        axis.title.y = element_text(face="bold", color="black",size=15),
        axis.line = element_line(size = 1)) +
  ylim(c(0,1))

all_data <- as.data.frame(all_datas['CVs.red'])  
colnames(all_data) <- c('Raw-Data','Nadir-Redv', 'Nadir-KD', 'Total-FPAR', expression('Total-i'[0]))
all_data <- melt(all_data)

plot_wheat_red <- ggplot(all_data, mapping = aes(x = variable, y = value, fill = variable)) + ylab(expression(CV)) +  ggtitle('Wheat - Red SIF') + 
  #geom_violin(aes(fill = variable),trim=False, size = 2.5, show.legend = FALSE, adjust = .8) + 
  stat_boxplot(geom ='errorbar', width = 0.6, size = 1, show.legend = FALSE) +
  geom_boxplot(notch=T, size = 0.8, show.legend = FALSE, width = 0.6, outlier.shape = NA) + 
  geom_jitter( stroke = 0, alpha = .2, fill = 'black', size = 2, width = 0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", color = 'white', alpha = 1, size=2, show.legend = F) +  
  scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels = c('Raw-Data','Nadir-Redv', 'Nadir-KD', 'Total-FPAR', expression(Total-i[0]))) + 
  #scale_color_brewer(palette="Dark2") + 
  theme_classic() + 
  theme(plot.title = element_text(face="bold", color="black",size=15, angle=0),
        axis.title.x=element_blank(),
        axis.text.x = element_text( color="black",size=10, angle=0),
        axis.text.y = element_text( color="black", size=10, angle=0),
        axis.title.y = element_text(face="bold", color="black",size=15),
        axis.line = element_line(size = 1)) +
  ylim(c(0,1))


# import data
all_datas <- readMat('results/CVs_corn.mat')

all_data <- as.data.frame(all_datas['CVs.far'])  
colnames(all_data) <- c('Raw-Data','Nadir-NIRv', 'Nadir-KD', 'Total-FPAR', expression('Total-i'[0]))
all_data <- melt(all_data)

plot_corn_far <- ggplot(all_data, mapping = aes(x = variable, y = value, fill = variable)) + ylab(expression(CV)) +  ggtitle('Corn - Far-red SIF') + 
  #geom_violin(aes(fill = variable),trim=False, size = 2.5, show.legend = FALSE, adjust = .8) + 
  stat_boxplot(geom ='errorbar', width = 0.6, size = 1, show.legend = FALSE) +
  geom_boxplot(notch=T, size = 0.8, show.legend = FALSE, width = 0.6, outlier.shape = NA) + 
  geom_jitter( stroke = 0, alpha = .2, fill = 'black', size = 2, width = 0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", color = 'white', alpha = 1, size=2, show.legend = F) +  
  scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels = c('Raw-Data','Nadir-NIRv', 'Nadir-KD', 'Total-FPAR', expression(Total-i[0]))) + 
  #scale_color_brewer(palette="Dark2") + 
  theme_classic() + 
  theme(plot.title = element_text(face="bold", color="black",size=15, angle=0),
        axis.title.x=element_blank(),
        axis.text.x = element_text( color="black",size=10, angle=0),
        axis.text.y = element_text( color="black", size=10, angle=0),
        axis.title.y = element_text(face="bold", color="black",size=15),
        axis.line = element_line(size = 1)) +
  ylim(c(0,1))

all_data <- as.data.frame(all_datas['CVs.red'])  
colnames(all_data) <- c('Raw-Data','Nadir-Redv', 'Nadir-KD', 'Total-FPAR', expression('Total-i'[0]))
all_data <- melt(all_data)

plot_corn_red <- ggplot(all_data, mapping = aes(x = variable, y = value, fill = variable)) + ylab(expression(CV)) +  ggtitle('Corn - Red SIF') + 
  #geom_violin(aes(fill = variable),trim=False, size = 2.5, show.legend = FALSE, adjust = .8) + 
  stat_boxplot(geom ='errorbar', width = 0.6, size = 1, show.legend = FALSE) +
  geom_boxplot(notch=T, size = 0.8, show.legend = FALSE, width = 0.6, outlier.shape = NA) + 
  geom_jitter( stroke = 0, alpha = .2, fill = 'black', size = 2, width = 0.1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="point", color = 'white', alpha = 1, size=2, show.legend = F) +  
  scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels = c('Raw-Data','Nadir-Redv', 'Nadir-KD', 'Total-FPAR', expression(Total-i[0]))) + 
  #scale_color_brewer(palette="Dark2") + 
  theme_classic() + 
  theme(plot.title = element_text(face="bold", color="black",size=15, angle=0),
        axis.title.x=element_blank(),
        axis.text.x = element_text( color="black",size=10, angle=0),
        axis.text.y = element_text( color="black", size=10, angle=0),
        axis.title.y = element_text(face="bold", color="black",size=15),
        axis.line = element_line(size = 1)) +
  ylim(c(0,1))
plot_figure <- plot_grid(   plot_wheat_far,plot_wheat_red,
                            plot_corn_far,plot_corn_red,
                            labels = c("a", "b","c","d"),
                            label_size = 20,
                            nrow = 2, ncol = 2, align = 'v')
plot_figure
ggsave("plot/figure_1_CVs.tiff", plot = plot_figure, width = 25, height = 15, units = "cm", dpi = 600, limitsize = FALSE, compression = "lzw")
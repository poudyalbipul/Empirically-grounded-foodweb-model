library(tidyverse)
library(readxl)
library(patchwork)
newdata <- read_excel("community_outcomes.xlsx")
df <- newdata[c("Diversity","Total_biomass", "Link_density",
                "Stability", "Mean_trophic_level")]

plotDB <- ggplot(newdata,aes(x= Diversity,y= Total_biomass)) +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  labs(title = "", x = "", y = "Biomass") +
  stat_function(fun = function(x) (x*0.4758)/(1+ x*0.0369),
                color = "red", size = 1)+
  scale_x_continuous(limits = c(0,80),
                     breaks = c(0,40,80),
                     labels = NULL)+
  scale_y_continuous(limits =c(4,20))+
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  annotate("text",
           x = 65, y = 18,   # adjust position as needed
           label = "\u03C1 = 0.65",
           size = 3)

plotDS <- ggplot(newdata,aes(x= Diversity,y= -Stability)) +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  labs(title = "", x = "Species richness", y = "Stability") +
  stat_function(fun = function(x) 0.017445 - 0.0022*sqrt(x),
                color = "red", size = 1)+
  scale_x_continuous(limits = c(0, 80),
                     breaks = c(0,40,80)
                     )+
  scale_y_continuous(limits =c(0,0.1),
                     breaks = c(0,0.05, 0.1))+
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  annotate("text",
           x = 65, y = 0.09,   # adjust position as needed
           label = "\u03C1 = -0.24",
           size = 3)

plotDL <- ggplot(newdata,aes(x= Diversity,y= Link_density)) +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  stat_function(fun = function(x) 0.05988*x,
                color = "red", size = 1)+
  labs(title = "", x = "", y = "Link density")+
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  scale_x_continuous(limits = c(0, 80),
                     breaks = c(0,40,80),
                     labels = NULL)+
  scale_y_continuous(limits =c(0,4),
                     breaks = c(0,2,4))+
  annotate("text",
           x = 65, y = 3,   # adjust position as needed
           label = "\u03C1 = 0.88",
           size = 3)


plotDT <- ggplot(newdata,aes(x= Diversity,y= Mean_trophic_level)) +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  stat_function(fun = function(x)  log(x)-1,
                color = "red", size = 1)+
  labs(title = "", x = "", y = "Mean trophic level") +
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  scale_x_continuous(limits = c(0, 80),
                     breaks = c(0,40,80),
                     labels = NULL)+
  scale_y_continuous(limits =c(1,4),
                     breaks = c(1,2,3,4))+
  annotate("text",
           x = 65, y = 3.6,   # adjust position as needed
           label = "\u03C1 = 0.73",
           size = 3)

plotTB <- ggplot(newdata, aes(x= Mean_trophic_level,y= Total_biomass)) +
  labs(title = "", x = "", y = "") +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  stat_function(fun = function(x) (2.17*exp(x))/(1+0.19977*exp(x)) ,
                color = "red", size = 1)+
  scale_y_continuous(limits = c(5, 20),
                     labels = NULL)+
  theme_bw()+
  scale_x_continuous(limits =c(1,4),
                     breaks = c(1,2,3,4),
                     labels = NULL)+
  annotate("text",
           x = 3.25, y = 18,   # adjust position as needed
           label = "\u03C1 = 0.61",
           size = 3)

plotLB <- ggplot(newdata,aes(x= Link_density, y= Total_biomass)) +
  labs(title = "", x = "", y = "") +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  stat_function(fun = function(x) 0.285*x/(0.06+x*0.0097),
                color = "red", size = 1)+
  scale_y_continuous(limits = c(5, 20),
                     labels = NULL)+
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  scale_x_continuous(limits =c(0,4),
                     breaks = c(0,2,4),
                     labels = NULL)+
  annotate("text",
           x = 3.25, y = 18,   # adjust position as needed
           label = "\u03C1 = 0.59",
           size = 3)

plotBS <- ggplot(newdata,aes(x= Total_biomass,y= -Stability)) +
  labs(title = "", x = "Biomass", y = "") +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  stat_function(fun = function(x) 0.017445- 0.0022*sqrt((x)/(0.4758-x*0.0369)),
                color = "red", size = 1)+
  scale_x_continuous(limits = c(5, 20)
                     )+
  scale_y_continuous(limits =c(0,0.1),
                     breaks = c(0, 0.05, 0.1),
                     labels = NULL)+
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  annotate("text",
           x = 16.25, y = 0.09,   # adjust position as needed
           label = "\u03C1 = -0.09",
           size = 3)

plotTS <- ggplot(newdata,aes(x= Mean_trophic_level,y= -Stability)) +
  labs(title = "", x = "Mean trophic level", y = "") +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  stat_function(fun = function(x) 0.015 - 0.0033*sqrt(exp(x)),
                color = "red", size = 1)+
  scale_x_continuous(limits = c(1, 4),
                     breaks = c(1,2,3,4)
                     )+
  scale_y_continuous(limits =c(0,0.1),
                     breaks = c(0,0.05, 0.1),
                     labels = NULL)+
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  annotate("text",
           x = 3.25, y = 0.09,   # adjust position as needed
           label = "\u03C1 = -0.18",
           size = 3)

plotLT <- ggplot(newdata,aes(x= Link_density,y= Mean_trophic_level)) +
  labs(title = "", x = "", y = "") +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  stat_function(fun = function(x) log(x/0.08)-1,
                color = "red", size = 1)+
  scale_y_continuous(limits = c(1, 4),
                     breaks = c(1,2,3,4),
                     labels = NULL)+
  scale_x_continuous(limits =c(0,4),
                     breaks = c(0,2,4),
                     labels = NULL)+
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  annotate("text",
           x = 3.25, y = 3.6,   # adjust position as needed
           label = "\u03C1 = 0.78",
           size = 3)

plotLS <- ggplot(newdata,aes(x= Link_density ,y= -Stability )) +
  labs(title = "", x = "Link density", y = "") +
  geom_point(size = 0.5, alpha = 0.1, color = "blue")+
  #stat_function(fun = function(x) (0.0177-x)^2 /(0.0097)^2,
  #              color = "red", size = 1, xlim = c(0,0.0177))+
  stat_function(fun = function(x) 0.0177 - 0.0097*(sqrt(x)) ,
                color = "red", size = 1)+
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = c(0, 0.05, 0.1),
                     labels = NULL)+
  scale_x_continuous(limits =c(0,4),
                     breaks = c(0,2,4)
                     )+ 
  theme_bw()+
  theme(axis.title = element_text(size = 16))+
  annotate("text",
           x = 3.25, y = 0.09,   # adjust position as needed
           label = "\u03C1 = -0.22",
           size = 3)

hist_div <- ggplot(df, aes(x = Diversity, y = after_stat(count / max(count)))) +
  geom_histogram(bins=30, fill = "blue", alpha = 0.6, color = "black") +
  labs(x = "", y = "Frequency") +
  # theme(axis.text.y = element_blank(),
  #      axis.ticks.y = element_blank())+
  scale_x_continuous(limits = c(0, 80),
                     breaks = c(0,40,80),
                     labels = NULL)+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme_bw()+
  theme(axis.title = element_text(size = 16))

hist_pop <- ggplot(df, aes(x = Total_biomass, y = after_stat(count / max(count)))) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.6, color = "black") +
  labs(x = "", y = "") +
  #theme(axis.text.y = element_blank(),
  #     axis.ticks.y = element_blank())+
  scale_x_continuous(limits = c(5, 20),
                     labels = NULL)+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme_bw()

hist_stab <- ggplot(df, aes(x = -Stability, y = after_stat(count / max(count)))) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.6, color = "black") +
  labs(x = "Stability", y = "") +
  #theme(axis.text.y = element_blank(),
  #     axis.ticks.y = element_blank())+
  scale_x_continuous(limits = c(0, 0.1),
                     breaks = c(0,0.05,0.1),
                     labels = NULL)+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme_bw()+
  theme(axis.title = element_text(size = 16))


hist_links <- ggplot(df, aes(x = Link_density, y = after_stat(count / max(count)))) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.6, color = "black") +
  labs(x = "", y = "") +
  #theme(axis.text.y = element_blank(),
  #     axis.ticks.y = element_blank())+
  scale_x_continuous(limits = c(0, 4),
                     breaks = c(0,2,4)
                     )+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme_bw()+
  theme(axis.title = element_text(size = 16))

hist_tl <- ggplot(df, aes(x = Mean_trophic_level, y = after_stat(count / max(count)))) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.6, color = "black") +
  labs(x = "", y = "") +
  #theme(axis.text.y = element_blank(),
  #     axis.ticks.y = element_blank())+
  scale_x_continuous(limits = c(1,4),
                     labels = NULL)+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme_bw()

all_plots <- list( hist_div, plot_spacer(), plot_spacer(), plot_spacer(),plot_spacer(),
                  plotDL, hist_links, plot_spacer(), plot_spacer(), plot_spacer(), 
                  plotDT, plotLT, hist_tl, plot_spacer(), plot_spacer(), 
                  plotDB, plotLB, plotTB, hist_pop, plot_spacer(), 
                  plotDS, plotLS, plotTS, plotBS, hist_stab) 

final_plot <- wrap_plots(all_plots, ncol = 5, nrow = 5)
final_plot

setwd("~/Documents/Research/loZIBpaper/R files")

library(tidyverse)
library(ggpubr)

sbirt <- readRDS("sbirt_clean.rds")




# visit1 <- na.omit(futimes) %>%
#   ggplot( aes(x=futimes[,2])) +
#   geom_histogram(alpha=0.6, position = 'identity',binwidth = 1) +
#   coord_cartesian(ylim=c(0,60), xlim=c(0,30)) +
#   theme(legend.position="none",
#         panel.background = element_blank(),
#         axis.title.x = element_text(size=30),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(size = 20),
#         axis.text.y = element_text(size = 20),
#         panel.grid.major=element_line(colour="grey90"),
#         title = element_text(size = 30)) +
#   #panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
#   labs(x = "Time Since Previous Visit", y = "Count", title = "3 Month") +
#   annotate(geom="text",x=15, y=45, label="Number of zero visits = 374", size=7)
# #  scale_fill_manual(values=group.colors) 
# #geom_smooth(aes(data.frame(x=c(0:10)),y=dpois(x, 2.5)), colour="red")
# #  geom_smooth(method="glm",method.args = list(family = "poisson"), formula = count(data$prcarevisit)~1)
# #scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
# 


sum(sbirt$visit==1 & sbirt$heavy==0, na.rm=T)
sum(sbirt$visit==2 & sbirt$heavy==0, na.rm=T)
sum(sbirt$visit==3 & sbirt$heavy==0, na.rm=T)
sum(sbirt$visit==4 & sbirt$heavy==0, na.rm=T)


# visit1

inflationplot_heavy_v1 <-   sbirt %>% filter(!is.na(heavy), visit==1) %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), alpha=.6, position = 'identity', binwidth = 1, color="white") +
  coord_cartesian(ylim = c(0, 85)) +
  theme(legend.position="none",
                 panel.background = element_blank(),
                 axis.title.x = element_text(size=12),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(size = 12),
                 axis.text.y = element_text(size = 12),
                 panel.grid.major=element_line(colour="grey90"),
                 title = element_text(size = 12),
                  plot.title = element_text(hjust = .5),
                  panel.border = element_rect(colour = NA, fill=NA, linewidth =.5)) +
     labs(x = "Days of Heavy Drinking", y = "Count", title = "Baseline") +
     annotate(geom="text",x=50, y=50, label="Zeros = 270", size=4)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

inflationplot_heavy_v1

#ggsave("inflationheavybase.pdf", plot=inflationplot_heavy_v1)



#Visit 2
inflationplot_heavy_v2 <-   sbirt %>% filter(!is.na(heavy), visit==2) %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), alpha=.6, position = 'identity', binwidth = 1, color="black") +
  coord_cartesian(ylim = c(0, 85)) +
  theme(legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 12),
        plot.title = element_text(hjust = .5),
        panel.border = element_rect(colour = NA, fill=NA, linewidth =.5)) +
  labs(x = "Days of Heavy Drinking", y = "Count", title = "3 Months") +
  annotate(geom="text",x=50, y=50, label="Zeros = 374", size=4)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

inflationplot_heavy_v2


#Visit 3
inflationplot_heavy_v3 <-   sbirt %>% filter(!is.na(heavy), visit==3) %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), alpha=.6, position = 'identity', binwidth = 1, color="black") +
  coord_cartesian(ylim = c(0, 85)) +
  theme(legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 12),
        plot.title = element_text(hjust = .5),
        panel.border = element_rect(colour = NA, fill=NA, linewidth =.5)) +
  labs(x = "Days of Heavy Drinking", y = "Count", title = "6 Months") +
  annotate(geom="text",x=50, y=50, label="Zeros = 396", size=4)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

inflationplot_heavy_v3



#Visit 4
inflationplot_heavy_v4 <-   sbirt %>% filter(!is.na(heavy), visit==4) %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), alpha=.6, position = 'identity', binwidth = 1, color="black") +
  coord_cartesian(ylim = c(0, 85)) +
  theme(legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 12),
        plot.title = element_text(hjust = .5),
        panel.border = element_rect(colour = NA, fill=NA, linewidth =.5)) +
  labs(x = "Days of Heavy Drinking", y = "Count", title = "12 Months") +
  annotate(geom="text",x=50, y=50, label="Zeros = 438", size=4)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

inflationplot_heavy_v4


library(grid)

inflation_all_plot <- ggarrange(
  inflationplot_heavy_v1 + rremove("ylab") + rremove("xlab"),
  inflationplot_heavy_v2 + rremove("ylab") + rremove("xlab"),
  inflationplot_heavy_v3 + rremove("ylab") + rremove("xlab"), 
  inflationplot_heavy_v4 + rremove("ylab") + rremove("xlab"),
  nrow = 1
) 

annotate_figure(inflation_all_plot, left = textGrob("Count", rot = 90, vjust = 1, gp = gpar(cex = 1.3, fontsize=10)),
                bottom = textGrob("Days of Heavy Drinking", gp = gpar(cex = 1.3, fontsize = 10)))

ggsave("/Users/bwrogers/Documents/Research/loZIBpaper/Latex Files/inflation_all_plot.pdf", width=6, height=3 )




# Visit 2
inflationplot_heavy_v2 <-   sbirt %>% filter(!is.na(heavy), visit==2) %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), position = 'identity', binwidth = 1, alpha=.6, color="black") +
  coord_cartesian(ylim = c(0, 65)) +
  theme(legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=30),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 30),
        panel.border = element_rect(colour = NA, fill=NA, size=.5)) +
  labs(x = "Days of Heavy Drinking", y = "Count", title = "Baseline") +
  annotate(geom="text",x=50, y=45, label="Number of zero visits = 374", size=10)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

inflationplot_heavy_v2

#visit 3

inflationplot_heavy_v3 <-   sbirt %>% filter(!is.na(heavy), visit==3) %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), position = 'identity', binwidth = 1, alpha=.6, color="black") +
  coord_cartesian(ylim = c(0, 65)) +
  theme(legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=30),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 30),
        panel.border = element_rect(colour = NA, fill=NA, size=.5)) +
  labs(x = "Days of Heavy Drinking", y = "Count", title = "Baseline") +
  annotate(geom="text",x=50, y=45, label="Number of zero visits = 396", size=7)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

inflationplot_heavy_v3


inflationplot_heavy_v4 <-   sbirt %>% filter(!is.na(heavy), visit==4) %>% 
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), position = 'identity', binwidth = 1, alpha=.6, color="black") +
  coord_cartesian(ylim = c(0, 65)) +
  theme(legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=30),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 30),
        panel.border = element_rect(colour = NA, fill=NA, size=.5)) +
  labs(x = "Days of Heavy Drinking", y = "Count", title = "Baseline") +
  annotate(geom="text",x=50, y=45, label="Number of zero visits = 438", size=7)+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

inflationplot_heavy_v4



#### All inflation plots with geom_label


heavy_labs <- tibble(visit_string = factor(c("Baseline", "3 Month", "6 Month", "12 Month"), levels=c("Baseline", "3 Month", "6 Month", "12 Month")),
                     x=rep(25,4),
                     y= rep(50,4),
                     label = c("Zeros = 270",
                               "Zeros = 374",
                               "Zeros = 396",
                               "Zeros = 438"
                               ))

inflationplot_heavy_all <-   sbirt %>% filter(!is.na(heavy)) %>% 
  mutate(visit_string = case_when(
    visit == 1 ~ "Baseline",
    visit == 2 ~ "3 Month",
    visit == 3 ~ "6 Month",
    visit == 4 ~ "12 Month"
  )) %>%
  mutate(visit_string = factor(visit_string, levels = c("Baseline", "3 Month", "6 Month", "12 Month"))) %>%
  ggplot() +
  geom_histogram(mapping=aes(x=heavy), alpha=.6, position = 'identity', binwidth = 1, color="white") +
  coord_cartesian(ylim = c(0, 85)) +
  theme(legend.position="none",
        panel.background = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major=element_line(colour="grey90"),
        title = element_text(size = 10),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5),
        panel.border = element_rect(colour = NA, fill=NA, linewidth =.5)) +
  labs(x = "Days of Heavy Drinking", y = "Count", title = "Days of Heavy Drinking", subtitle= "SBIRT Study") +
  #annotate(geom="text",x=50, y=50, label="Zeros = 270", size=5)+
  geom_label(data=heavy_labs, aes(label=label, x=x, y=y, size=10), label.size = NA) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  facet_wrap(~visit_string, nrow = 2) +
  theme(
    strip.text = element_text( size = rel(1.2)),
    strip.background = element_rect(fill = "white", colour = "white", size = 1)
  )

inflationplot_heavy_all


ggsave("/Users/bwrogers/Documents/Research/loZIBpaper/Latex Files/inflation_all_plot.pdf", width=5, height=5 )



#################### LOO Plot

loomodels <- c("UN", "AD", "UNcv", "ADcv", "ARcv", "AR", "RI", "1RI")



loo_results <- data.frame(Model= factor(loomodels, levels=rev(loomodels)),
                          LOO = c(7515.7, 7550.0, 7561.6,7564.8, 7584.0, 7616.8,  16829.6),
                          SE = c(147.0, 146.2, 150.4, 150.8, 151.5, 152.2, 815.1)
                          )

loo_plot_dat <- loo_results %>% 
  mutate(upper=LOO + 1.96 * SE, lower = LOO - 1.96*SE)


model_name_color <- c("blue", rep("orange",6))

loo_plot <- ggplot(loo_plot_dat, aes(x=LOO, y= Model)) +
    geom_errorbarh(aes(xmax=upper, xmin=lower), linewidth=1, height = .3, color = rev(model_name_color)) +
    geom_point(size=5, color=rev(model_name_color)) +
    labs(title =  "LOO Estimates for Each Model"
         #subtitle = "Lower LOO indicates better fit",
         #caption = "*Red indicates OLRE models"
         ) +
    theme_bw() +
    theme(axis.text.y = element_text(color= model_name_color, size=30),
          panel.background = element_blank(),
          axis.title.x = element_text(size=30),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 30),
          panel.grid.major=element_line(colour="grey90"),
          title = element_text(size = 30))
  
loo_plot

ggsave("looplot.pdf", plot=loo_plot )

#Re-save plots in Latex files folder
setwd("/Users/bwrogers/Documents/Research/loZIBpaper/Latex Files")

ggsave("looplot.pdf", plot=loo_plot )

ggsave("inflationheavybase.pdf", plot=inflationplot_heavy_v1)

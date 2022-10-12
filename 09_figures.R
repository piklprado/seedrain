library(ggplot2)
library(ggthemes)
library(gridExtra)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggrepel)


################################################################################
## Ggplot themes
## From https://rpubs.com/Koundy/71792
################################################################################
theme_Publication <- function(base_size=14, base_family="helvetica") {
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.2), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                ##legend.margin = unit(0, "cm"),
                legend.spacing = unit(0.2, "cm"),
                legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
                ))
    
}

scale_fill_Publication <- function(...){
    library(scales)
    discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f",
                                                              "#ef3b2c","#662506","#a6cee3",
                                                              "#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
    library(scales)
    discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f",
                                                                "#ef3b2c","#662506","#a6cee3",
                                                                "#fb9a99","#984ea3","#ffff33")), ...)

}

################################################################################
## Tree frequency distribution across months and traps (box plots)
################################################################################
## Traps
p1 <-
    tudo82spp %>%
    dplyr::select(species, nt_p1, nt_p2, nt_p3, nt_p_tot40) %>%
    pivot_longer(-species, names_to = "period", values_to = "n.traps") %>%
    ggplot(aes(period, n.traps)) +
    geom_boxplot() +
    theme_Publication() +
    xlab("") +
    ylab("Number of traps") +
    scale_y_continuous(limits = c(0,40), expand = c(0.01, 0.01)) +
    scale_x_discrete(labels=c("Pooled", "Year 1", "Year 2", "Year 3"))
## Months
p2 <-
    tudo82spp %>%
    dplyr::select(species, nm_p1, nm_p2, nm_p3, nm_p_tot12) %>%
    pivot_longer(-species, names_to = "period", values_to = "n.months") %>%
    ggplot(aes(period, n.months)) +
    geom_boxplot() +
    theme_Publication() +
    xlab("") +
    ylab("Number of months") +
    scale_y_continuous(limits = c(0,12), expand = c(0.01, 0.01)) +
    scale_x_discrete(labels=c("Pooled", "Year 1", "Year 2", "Year 3"))

## A single figure with both as panels
p12 <- arrangeGrob(p1, p2, ncol = 2, nrow = 1)
ggsave(filename= "figures/frequency_boxplots.pdf", plot = p12, width =8, height = 4, device=cairo_pdf)
ggsave(filename= "figures/frequency_boxplots.png", plot = p12, width =8, height = 4)

ggsave(filename= "figures/frequency_boxplots.tiff", plot = p12, width =8, height = 4, dpi=300)
ggsave(filename= "figures/frequency_boxplots.eps", plot = p12, width =8, height = 4, device=cairo_ps)


################################################################################
## Synchrony null models limits x observed values
################################################################################

p3 <-
    ggplot(summary.dmed, aes(year, mean.ov)) +
    geom_point(data = observed, aes(y = mean.ov, shape = "Observed"),
               size = 3, position=position_nudge(x = -0.3)) +
    geom_errorbar(aes(ymin = mean.ov.low, ymax = mean.ov.upp, color = by, lwd= by),
                  position = position_dodge(width=0.2)) +
    ylab("Mean temporal asynchrony (Bray Curtis)") +
    xlab("") +
    labs(color = "Null model:", shape = "", size = "Null model:", linetype = "Null model:") +
    scale_color_manual(values = c("black",  "darkgrey", "lightgrey")) +
    ##scale_linetype_manual(values = c("11","1","11")) +
    scale_size_manual(values = c(1.5,1.5,1.5)) +
    theme_Publication() 

ggsave(filename= "figures/null_models_synchrony.pdf", plot = p3, width = 9, height = 7.5, device=cairo_pdf)
ggsave(filename= "figures/null_models_synchrony.png", plot = p3, width = 9, height = 7.5)

ggsave(filename= "figures/null_models_synchrony.eps", plot = p3, width = 9, height = 7.5, device=cairo_ps)
ggsave(filename= "figures/null_models_synchrony.tiff", plot = p3, width = 9, height = 7.5, dpi = 300)

################################################################################
## SSL predicted by glmms and observed
################################################################################
## In separated panels
p4a <-
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc,
                                 paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq,
                               paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq))) %>%
    ggplot(aes(mass, ssl.mean)) +
    geom_ribbon(aes(y = pfit, ymin = plower, ymax = pupper), data = ssl.pred, fill="gray", alpha=0.75) +
    geom_point(size=2) +
    geom_linerange(aes(ymin=ssl.min, ymax = ssl.max)) +
    geom_line(aes(y = pfit), data = ssl.pred)+
    scale_x_log10() +
    ylab("SSL") +
    xlab("Seed mass (g)") +
    facet_grid(height.class ~ freq.class) +
    theme_Publication()
## In a single panel
ab.sp.rsl %<>%
    mutate(height.class = ifelse(height <= median.Hloc,
                                 paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq,
                               paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq)),
           classe = paste(height.class,freq.class, sep =" , "))
                      
p4b <-
    ssl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_point(data = ab.sp.rsl, aes(y=ssl.mean, color = classe), size=2) +
    geom_line(aes(y=pfit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = plower, ymax =pupper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw() +
    ylab("SSL")+
    xlab("Seed mass (g)") +
    theme_Publication() +
    theme(legend.position = c(0.75, 0.25), legend.title = element_blank(), legend.direction = "vertical")


ggsave(filename= "figures/SSL_pred_prob.pdf", plot = p4a, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/SSL_pred.prob.png", plot = p4a, width = 9, height = 9)
ggsave(filename= "figures/SSL_pred_prob_single.pdf", plot = p4b, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/SSL_pred.prob_single.png", plot = p4b, width = 9, height = 9)

ggsave(filename= "figures/SSL_pred_prob.eps", plot = p4a, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/SSL_pred.prob.tiff", plot = p4a, width = 9, height = 9, dpi = 300)
ggsave(filename= "figures/SSL_pred_prob_single.eps", plot = p4b, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/SSL_pred.prob_single.tiff", plot = p4b, width = 9, height = 9, dpi = 300)


################################################################################
## SSL Predicted by glms and observed
################################################################################
p5  <-
    ggplot(abundants3, aes(mass)) +
    geom_point(aes(y=ssl_tot40)) +
    geom_line(data = ab.ssl.tot40.newdata, aes(mass, fit.resp)) +
    geom_ribbon(data = ab.ssl.tot40.newdata, aes(mass, ymin=lower.resp, ymax =upper.resp), alpha =0.2) +
    facet_grid(Hloc.class ~ freq.class) +
    scale_x_log10() +
    ylab("SSL")+
    xlab("Seed mass (g)") +
    theme_Publication()

ggsave(filename= "figures/SSL_all_pred_prob_glm.pdf", plot = p5, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/SSL_all_pred.prob_glm.png", plot = p5, width = 9, height = 9)

ggsave(filename= "figures/SSL_all_pred_prob_glm.eps", plot = p5, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/SSL_all_pred.prob_glm.png", plot = p5, width = 9, height = 9, dpi = 300)


################################################################################
## TSL predicted by glmms and observed
################################################################################
## In separated panels
p6a <-
    ab.sp.rsl %>%
    mutate(height.class = ifelse(height <= median.Hloc,
                                 paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq,
                               paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq))) %>%
    ggplot(aes(mass, tsl.mean)) +
    geom_point(size=2) +
    geom_linerange(aes(ymin=tsl.min, ymax = tsl.max)) +
    geom_line(aes(y = pfit), data = tsl.pred)+
    geom_ribbon(aes(y = pfit, ymin = plower, ymax = pupper), data = tsl.pred, fill="gray", alpha=0.25) +
    facet_grid(height.class ~ freq.class) +
    scale_x_log10() +
    theme_bw() +
    ylab("TSL") +
    xlab("Seed mass (g)") +
    facet_grid(height.class ~ freq.class) +
    theme_Publication()

## In a single panel
ab.sp.rsl %<>%
    mutate(height.class = ifelse(height <= median.Hloc,
                                 paste0("Tree Height < ", median.Hloc, " m"),
                                 paste0("Tree Height > ", median.Hloc, " m")),
           freq.class = ifelse(freq <= median.freq,
                               paste0("Adult Frequency < ", median.freq),
                               paste0("Adult Frequency > ", median.freq)),
           classe = paste(height.class,freq.class, sep =" , "))
                      
p6b <-
    tsl.pred %>%
    mutate(classe = paste(height.class,freq.class, sep =" , ")) %>%
    ggplot(aes(x=mass)) +
    geom_point(data = ab.sp.rsl, aes(y=tsl.mean, color = classe), size=2) +
    geom_line(aes(y=pfit, color = classe), size=1.2) +
    geom_ribbon(aes(ymin = plower, ymax =pupper, fill = classe), alpha =0.1) +
    scale_x_log10() +
    theme_bw() +
    ylab("TSL")+
    xlab("Seed mass (g)") +
    theme_Publication() +
    theme(legend.position = c(0.75, 0.25), legend.title = element_blank(), legend.direction = "vertical")


ggsave(filename= "figures/TSL_pred_prob.pdf", plot = p6a, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/TSL_pred.prob.png", plot = p6a, width = 9, height = 9)
ggsave(filename= "figures/TSL_pred_prob_single.pdf", plot = p6b, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/TSL_pred.prob_single.png", plot = p6b, width = 9, height = 9)

ggsave(filename= "figures/TSL_pred_prob.eps", plot = p6a, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/TSL_pred.prob.tiff", plot = p6a, width = 9, height = 9, dpi = 300)
ggsave(filename= "figures/TSL_pred_prob_single.eps", plot = p6b, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/TSL_pred.prob_single.tiff", plot = p6b, width = 9, height = 9, dpi = 300)

################################################################################
## SSL Predicted by glms and observed
################################################################################
p7  <-
    ggplot(abundants3, aes(mass)) +
    geom_point(aes(y=tsl_tot12)) +
    geom_line(data = ab.tsl.tot12.newdata, aes(mass, fit.resp)) +
    geom_ribbon(data = ab.tsl.tot12.newdata, aes(mass, ymin=lower.resp, ymax =upper.resp), alpha =0.2) +
    facet_grid(Hloc.class ~ freq.class) +
    scale_x_log10() +
    ylab("TSL")+
    xlab("Seed mass (g)") +
    theme_Publication()

ggsave(filename= "figures/TSL_all_pred_prob_glm.pdf", plot = p7, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/TSL_all_pred.prob_glm.png", plot = p7, width = 9, height = 9)

ggsave(filename= "figures/TSL_all_pred_prob_glm.eps", plot = p7, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/TSL_all_pred.prob_glm.tiff", plot = p7, width = 9, height = 9, dpi = 300)

################################################################################
## TSL x SSL relationship
################################################################################
## Joining data for each year and of pooled years
tmp1 <-
    abundants2 %>%
    mutate(year = "Pooled", ssl = ssl_tot40, tsl = tsl_tot12) %>%
    dplyr::select(species, sp.i, ssl, tsl, year) %>%
    bind_rows(abundants[, c("species", "sp.i", "ssl", "tsl", "year")]) %>%
    mutate(year = ifelse(grepl("Pooled", year), year, paste("Year",year))) %>%
    arrange(year, ssl, tsl, sp.i)
tmp1$dupl <- duplicated(tmp1[, c("year", "ssl", "tsl")])
## For the figure caption of the 2nd plot below: Species that share the same point in each panel
tmp1 %>%
    group_by(year,ssl,tsl) %>%
    summarise(N = n(), sigla = paste(sp.i, collapse =",")) %>%
    filter(N>1) %>%
    data.frame()
## The plots
p8 <-
    tmp1 %>%
    ggplot(aes(ssl, tsl, label = sp.i)) +
    geom_point() +
    geom_text_repel(alpha =.3, size = 3, max.overlaps = 20) +
    xlab("SSL") +
    ylab("TSL") +
    xlim(0,1.05) +
    ylim(0,1.05) +
    theme_Publication()
## All points and all years (a bit messy)
p8a <- p8 + facet_wrap(~year, nrow=2)
## select a single point and label for overlaps
p8b <- p8 %+% filter(tmp1, !dupl) + facet_wrap(~year, nrow=2)
## Only pooled year, all points
p8c <- p8 %+% filter(tmp1, year == "Pooled")


ggsave(filename= "figures/TSLxSSL_all_spp.pdf", plot = p8a, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/TSLxSSL_all_spp.png", plot = p8a, width = 9, height = 9)
ggsave(filename= "figures/TSLxSSL_ommit_overlaps.pdf", plot = p8b, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/TSLxSSL_ommit_overlaps.png", width = 9, height = 9, plot = p8b)
ggsave(filename= "figures/TSLxSSL_pooled.pdf", plot = p8c, width = 9, height = 9, device=cairo_pdf)
ggsave(filename= "figures/TSLxSSL_pooled.png", width = 9, height = 9, plot = p8c)

ggsave(filename= "figures/TSLxSSL_all_spp.eps", plot = p8a, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/TSLxSSL_all_spp.tiff", plot = p8a, width = 9, height = 9, dpi = 300)
ggsave(filename= "figures/TSLxSSL_ommit_overlaps.eps", plot = p8b, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/TSLxSSL_ommit_overlaps.tiff", width = 9, height = 9, plot = p8b, dpi = 300)
ggsave(filename= "figures/TSLxSSL_pooled.eps", plot = p8c, width = 9, height = 9, device=cairo_ps)
ggsave(filename= "figures/TSLxSSL_pooled.tiff", width = 9, height = 9, plot = p8c, dpi = 300)

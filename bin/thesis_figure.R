#!/us/bin/env Rscript
# Generate figure on HOG pathway work for thesis
library(tidyverse)
library(ggpubr)
library(multipanelfigure)
library(Biostrings)
data("BLOSUM62")
blosum <- as_tibble(BLOSUM62, rownames = "wt") %>%
  pivot_longer(-wt, names_to = "mut", values_to = "blosum62")
  
blank_plot <- function(text = NULL){
  p <- ggplot(tibble(x=c(0, 1)), aes(x=x, y=x)) +
    geom_blank() +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
  
  if (!is.null(text)){
    p <- p + annotate(geom = 'text', x = 0.5, y = 0.5, label = text)
  }
}

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

### Panel  - P(Aff) ###
jelier <- read_tsv("data/jelier.impact") %>%
  select(gene, position = pos_aa, wt = ref_aa, mut = alt_aa, effect, foldx = foldx_ddG, sift = sift_score) %>%
  left_join(blosum, by = c("wt", "mut")) %>%
  mutate(sift_bin = cut(sift, 20, labels=FALSE),
         foldx_bin = cut(foldx, 20, labels=FALSE)) %>%
  group_by(sift_bin) %>%
  mutate(sift_prob = sum(effect == 0) / n()) %>%
  group_by(foldx_bin) %>%
  mutate(foldx_prob = sum(effect == 0) / n()) %>%
  ungroup()

min_sift <- min(jelier$sift[jelier$sift > 0], na.rm = TRUE)
f_sift <- function(x,a,b){1/(1 + exp(a*log(x + min_sift) + b))}
fit_sift <- nls(sift_prob ~ f_sift(sift,a,b), data = drop_na(jelier, sift), start = list(a=-1,b=1))

f_foldx <- function(x,a,b){1/(1 + exp(a*x + b))}
fit_foldx <- nls(foldx_prob ~ f_foldx(foldx,a,b), data = drop_na(jelier, foldx), start = list(a=0, b=0))

jelier_preds <- bind_rows(SIFT4G = tibble(x = seq(0, 1, 0.01), y = f_sift(x, coef(fit_sift)[1], coef(fit_sift)[2])),
                          FoldX = tibble(x = seq(-10, 40, 0.1), y = f_foldx(x, coef(fit_foldx)[1], coef(fit_foldx)[2])),
                          .id = "model")

jelier_probs <- bind_rows(SIFT4G = select(jelier, x = sift, y = sift_prob),
                          FoldX = select(jelier, x = foldx, y = foldx_prob),
                          .id = "model")
  
lab <- as_labeller(c(SIFT4G='"SIFT4G Score"', FoldX='"FoldX"~Delta*Delta*"G (kj~mol"[-1]*")"'), default = label_parsed)
p_paff <- ggplot(mapping = aes(x = x, y = y, colour = model)) +
  facet_wrap(~model, ncol = 1, scales = "free_x", labeller = lab, strip.position = "bottom") +
  geom_point(data = jelier_probs, shape = 20, show.legend = FALSE) +
  geom_line(data = jelier_preds, colour = "black") +
  scale_colour_manual(values = c(SIFT4G = "cornflowerblue", FoldX = "firebrick2")) +
  labs(x = "", y = "P(Neut)") +
  theme(strip.placement = "outside")


### Panel  - Gene Phenotype Associations ###
# Summarise associations between gene paff and phenoypes
p_associations <- blank_plot("P(Aff) Associations")

### Panel - HOG Pathway ###
p_hog <- ggplot(path_elements, aes(x = x, y = y, colour = colour)) +
  geom_point(shape = 15) +
  lims(x = c(0, 0.8), y = c(0, 1)) +
  annotation_raster(readPNG("meta/hog_pathway.png"), xmin = 0, xmax = 0.8, ymin = 0, ymax = 1, interpolate = TRUE) +
  scale_colour_manual(name = "", values = c(Adaptor="#6666FF", GTPase="#FF66FF", Kinase="#FF6666",
                                             Phosphatase="#66FF66", Sensor="#FFFF66", TF="#33FFCC")) +
  coord_fixed(expand = FALSE) +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10, 0,-10,-40))

### Panel - Pathway Models ###
# Summarise performance of different models - plot Rsquared, cor and prob?
p_pathway <- blank_plot("Pathway Models")

### Figure Assembly ###
size <- theme(text = element_text(size = 8))
p1 <- p_paff + labs(tag = 'A') + size
p2 <- p_associations + labs(tag = 'B') + size
p3 <- p_hog + labs(tag = 'C') + size
p4 <- p_pathway + labs(tag = 'D') + size

figure <- multi_panel_figure(width = 360, height = 240, columns = 3, rows = 2,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2:3) %>%
  fill_panel(p3, row = 2, column = 1) %>%
  fill_panel(p4, row = 2, column = 2:3)

ggsave('figures/thesis_figure.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/thesis_figure.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

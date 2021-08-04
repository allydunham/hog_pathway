#!/us/bin/env Rscript
# Generate figure on HOG pathway work for thesis
library(tidyverse)
library(ggpubr)
library(multipanelfigure)
library(ggtext)
library(png)

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

# Make pretty p value categories for plots
pretty_p_values <- function(p, breaks = c(0.001, 0.01, 0.05), markdown_exp=FALSE, prefix_p=FALSE){
  breaks <- sort(breaks, decreasing = TRUE)
  break_str <- as.character(breaks)
  gt <- '>'
  lt <- '<'
  
  if (markdown_exp){
    break_str <- str_replace(break_str, "^1?e(-?[0-9]*\\.?[0-9]*)", "10<sup>\\1</sup>")
    gt <- '&gt;'
    lt <- '&lt;'
  }
  
  p_out <- rep(str_c(gt, ' ', break_str[1]), length(p))
  for (i in 1:length(breaks)){
    p_out[p < breaks[i]] <- str_c(lt, ' ', break_str[i])
  }
  break_levels <- c(str_c(gt, ' ', break_str[1]), unlist(map(break_str, ~str_c(lt, ' ', .))))
  
  if (prefix_p){
    p_out <- str_c('p ', p_out)
    break_levels <- str_c('p ', break_levels)
  }
  
  return(factor(p_out, levels = break_levels))
}

#### Load Data ####
data("BLOSUM62", package = "Biostrings")
blosum <- as_tibble(BLOSUM62, rownames = "wt") %>%
  pivot_longer(-wt, names_to = "mut", values_to = "blosum62")

strains <- readRDS('data/Rdata/strain_meta.rds')
filtered_strains <- filter(strains, Ploidy <= 2, Aneuploidies == 'euploid') %>%
  pull(`Standardized name`) %>%
  setdiff(. , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))

genes <- readRDS('data/Rdata/gene_meta_all.rds')
essential <- readRDS('data/Rdata/essential_genes.rds')

# Take 72H where both available
conditions <- c(
  "2,4-Dichlorophenoxyacetic acid (48H)", "20ºC (72H)", "39ºC (72H)", "42ºC (72H)", "5-FU (48H)", 
  "6-AU (48H)", "6-AU + 39ºC (72H)", "aa starvation (48H)", "Acetic acid (48H)",
  "Amphotericin B (48H)", "Amphotericin B + anaerobic (48H)", "Anaerobic growth (48H)", "Cadmium chloride (48H)",
  "Caffeine 15mM (72H)", "Caffeine 20mM (72H)", "Caspofungin (48H)", "Clioquinol (72H)", 
  "Clozapine (48H)", "Cyclohexamide (72H)", "DMSO 1%  (48H)", "Glucose 20% (48H)", 
  "Glycerol 2%  (72H)", "NaCl 0.4M (72H)", "NaCl 0.4M + 39ºC (72H)", "NaCl 0.6M (72H)", 
  "NaCl 0.6M + 39ºC (72H)", "NiSO4 (48H)", "Nitrogen starvation (48H)", "Nystatin (48H)", "Paraquat (48H)", 
  "SC (48H)", "SC + hepes (48H)", "Sorbitol 1M (48H)", "YPD (24H)"
)

growth <- read_tsv('data/raw/yeasts_liti_fixed.tsv', col_names = TRUE, col_types = cols(strain = col_character())) %>%
  filter(info %in% filtered_strains, condition %in% conditions) %>% 
  filter(subset == 'liti') %>%
  rename(strain_id = strain, strain=info)

# Growth cor (comment out growth filtering)
# read_tsv('data/raw/yeasts_liti_fixed.tsv', col_names = TRUE, col_types = cols(strain = col_character())) %>%
#   filter(info %in% filtered_strains) %>% 
#   filter(subset == 'liti') %>%
#   rename(strain_id = strain, strain=info) %>%
#   select(strain, condition, score) %>%
#   extract(condition, c("condition", "time"), "([^\\(\\)]*) \\((48H|72H)\\)") %>%
#   drop_na() %>%
#   pivot_wider(names_from = time, values_from = score) %>%
#   drop_na() %>%
#   group_by(condition) %>%
#   group_modify(~broom::tidy(cor.test(.$`48H`, .$`72H`)))

ko_thresh <- 0.5
probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% growth$strain) %>%
  mutate(ko = p_aff > ko_thresh)

associations <- read_tsv("data/gene_paff_growth_assocations.tsv") %>%
  left_join(select(genes, gene = id, name), by = "gene") %>%
  left_join(select(essential, gene = locus, essential), by = "gene") %>%
  filter(condition %in% conditions)

path <- readRDS('data/Rdata/hog_path_probs.rds') %>%
  filter(strain %in% unique(growth$strain))

ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

#### Panel - P(Aff) ####
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
  
lab <- as_labeller(c(SIFT4G='"SIFT4G Score"', FoldX='"FoldX"~Delta*Delta*"G (kj"~"mol"[-1]*")"'), default = label_parsed)
p_paff <- ggplot(mapping = aes(x = x, y = y, colour = model)) +
  facet_wrap(~model, ncol = 1, scales = "free_x", labeller = lab, strip.position = "bottom") +
  geom_point(data = jelier_probs, shape = 20, show.legend = FALSE) +
  geom_line(data = jelier_preds, colour = "black") +
  scale_colour_manual(values = c(SIFT4G = "cornflowerblue", FoldX = "firebrick2")) +
  labs(x = "", y = "P(Neut)") +
  theme(strip.placement = "outside")

# BLOSUM Cor
# cor.test(jelier$blosum62, f_sift(jelier$sift, coef(fit_sift)[1], coef(fit_sift)[2]))
# cor.test(jelier$blosum62, f_foldx(jelier$foldx, coef(fit_foldx)[1], coef(fit_foldx)[2]))

#### Panel - Gene Phenotype Association Summary ####
association_summary <- group_by(associations, condition) %>%
  summarise(Essential = sum(p_adj < 0.01 & essential == "E", na.rm = TRUE),
            `Non-Essential` = sum(p_adj < 0.01 & essential == "NE", na.rm = TRUE)) %>%
  pivot_longer(-condition, names_to = "type", values_to = "count") %>%
  drop_na()

p_associations <- ggplot(association_summary, aes(x = reorder(condition, count), y = count, fill = type)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(name = "", values = c(Essential = "#4daf4a", `Non-Essential` = "#377eb8")) +
  scale_y_continuous(expand = expansion(0)) +
  labs(y = "Significant Genes (P<sub>adj</sub> < 0.01)", x = "") +
  theme(panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_markdown(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,-10,0),
        legend.key.size = unit(2, "mm"))

#### Panel - Example Gene Phenotype Associations ####
association_examples <- left_join(probs, select(genes, gene = id, name), by = "gene") %>%
  filter(name %in% c("HOG1", "PBS2", "SER3", "CDC42")) %>%
  left_join(filter(select(growth, condition, strain, score), condition == "NaCl 0.6M (72H)"), by = "strain") %>%
  mutate(ko_desc = ifelse(ko, "P(Aff) > 0.5", "P(Aff) ≤ 0.5"),
         name = factor(name, levels = c("HOG1", "PBS2", "CDC42", "SER3")))

p_association_examples <- ggplot(association_examples, aes(x = ko_desc, y = score, fill = ko_desc)) +
  facet_wrap(~name, nrow = 1) +
  geom_boxplot(outlier.shape = 20) +
  stat_compare_means(method = "wilcox", comparisons = list(c("P(Aff) > 0.5", "P(Aff) ≤ 0.5"))) +
  stat_summary(geom = "text", fun.data = function(x) {data.frame(y = -4.4, label = length(x))}) +
  scale_fill_manual(name = "", values = c(`P(Aff) > 0.5` = "firebrick2", `P(Aff) ≤ 0.5` = "cornflowerblue")) +
  labs(x = "", y = "S-Score") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.box.margin = margin(-20, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0))

### Panel - HOG Pathway ###
path_elements <- tibble(x = 2, y = 2, colour = c("Adaptor", "GTPase", "Kinase", "Phosphatase", "Sensor", "TF"))
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

#### Panel - HOG Pathway Models ####
# Summarise performance of different models - plot Rsquared, cor and prob?
assoc_gene_counts <- filter(associations, p_adj < 0.05, mean_ko < mean_active, essential == "NE") %>%
  select(condition, gene) %>%
  left_join(probs, by = "gene") %>%
  group_by(strain, condition) %>%
  summarise(sig_assoc_count = sum(ko), .groups = "drop")

ko_gene_counts <- filter(ko_growth, qvalue < 0.05) %>%
  select(condition, gene) %>%
  distinct() %>%
  left_join(probs, by = "gene") %>%
  drop_na() %>%
  group_by(strain, condition) %>%
  summarise(sig_ko_count = sum(ko, na.rm = TRUE), .groups = "drop")

path_models <- left_join(select(growth, condition, sscore = score, strain), select(path, strain, hog_probability), "strain") %>%
  left_join(assoc_gene_counts, by = c("strain", "condition")) %>%
  left_join(ko_gene_counts, by = c("strain", "condition")) %>%
  pivot_longer(c(-strain, -sscore, -condition), names_to = "model", values_to = "value") %>%
  drop_na() %>%
  group_by(model, condition) %>%
  group_modify(~broom::glance(lm(sscore ~ value, data = .))) %>%
  ungroup() %>%
  mutate(model = c(hog_probability = "P(HOG Pathway)", sig_assoc_count = "P(Aff) Association\nGene Sets",
                   sig_ko_count = "KO Growth\nGene Set")[model],
         model = factor(model, levels = c("P(HOG Pathway)", "P(Aff) Association\nGene Sets", "KO Growth\nGene Set")),
         p_cat = pretty_p_values(p.value, breaks = c(1e-10, 1e-5, 0.0001, 0.001, 0.01, 0.05, 0.1)))

condition_levels <- filter(path_models, model == "P(HOG Pathway)") %>% arrange(r.squared) %>% pull(condition)

p_pathway <- ggplot(path_models, aes(x = factor(condition, levels = condition_levels), y = r.squared, fill = p_cat)) +
  facet_wrap(~model, nrow = 1) +
  geom_col() +
  coord_flip() +
  scale_fill_brewer(name = "p-value", palette = "YlGnBu") +
  labs(x = "", y = expression(r^2)) +
  theme(panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        panel.grid.major.y = element_blank())

#### Figure Assembly ####
size <- theme(text = element_text(size = 12))
p1 <- p_paff + labs(tag = 'A') + size
p2 <- p_associations + labs(tag = 'B') + size
p3 <- p_association_examples + labs(tag = 'C') + size
p4 <- p_hog + labs(tag = 'D') + size
p5 <- p_pathway + labs(tag = 'E') + size

figure <- multi_panel_figure(width = 360, height = 240, columns = 3, rows = 2,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 1, column = 3) %>%
  fill_panel(p4, row = 2, column = 1) %>%
  fill_panel(p5, row = 2, column = 2:3)

ggsave('figures/thesis_figure_hog.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_hog.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

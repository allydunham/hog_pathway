#!/us/bin/env Rscript
# Generate figure on deletion phenotypes for thesis
library(tidyverse)
library(ggpubr)
library(multipanelfigure)
library(ggtext)
library(ggdendro)

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

blank_plot <- function(text = ''){
  ggplot(tibble(x=c(0, 1)), aes(x=x, y=x)) +
    geom_blank() +
    annotate(geom = 'text', x = 0.5, y = 0.5, label = text) +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

#### Load Data ####
gene_meta <- readRDS('data/Rdata/gene_meta_all.rds') %>%
  mutate(name = ifelse(is.na(name), id, name)) %>%
  rename(gene = id)

ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

# geno <- readRDS('data/Rdata/genotypes_all_genes.rds') %>%
#   select(UWOP = SACE_GAT, Y55 = ADA, YPS = AKN, S288C = SACE_GAV) %>%
#   mutate(across(everything(), ~ifelse(. > 0, 1, 0))) %>%
#   as.matrix() %>%
#   t()
# geno_dist <- dist(geno, method = "manhattan")

#### Panel - Strain Relationships - genetic and phenotypic (confusion) ####
# dendrgrams from genetic distance and KO profiles
# get_dend_order <- function(tbl, ...) {
#   mat <- tblhelpr::tibble_to_matrix(arrange(tbl, strain), -strain, row_names = "strain")
#   d <- dist(mat, method = "manhattan")
#   hc <- hclust(d)
#   out <- matrix(c(t(hc$merge)), nrow = 1)
#   return(as_tibble(out, .name_repair = ~str_c(str_c("p", c(1,1,2,2,3,3)), "_", c(1,2,1,2,1,2))))
# }
# 
# growth_wide <- mutate(ko_growth, sig = sign(score) * (qvalue < 0.01)) %>%
#   select(condition, strain, name, sig) %>%
#   pivot_wider(names_from = name, values_from = sig) %>%
#   filter(!condition == "Maltose 2% (72H)")
# 
# growth_dends <- group_by(growth_wide, condition) %>%
#   group_modify(., get_dend_order) %>%
#   group_by(p1_1, p1_2, p2_1, p2_2, p3_1, p3_2) %>%
#   summarise(n = n(), conditions = str_c(condition, collapse  = ", "), .groups = "drop")
# 
# get_dend <- function(con, ...) {
#   tbl <- filter(growth_wide, condition == con)
#   mat <- tblhelpr::tibble_to_matrix(arrange(tbl, strain), -strain, -condition, row_names = "strain")
#   d <- dist(mat, method = "manhattan")
#   hc <- hclust(d)
#   hc$height <- 1:3
#   dendro_data(hc)
# }
# 
# ggdendrogram(get_dend("39ºC (48H)"))
# 
# 
# pivot_longer(growth_dists, c(-n, -conditions), names_to = "pair", values_to = "ind") %>%
#   separate(pair, into = c("pair", "num")) %>%
#   pivot_wider(names_from = num, values_from = ind, names_prefix = "i") %>%
#   mutate()
# 
# ggplot(growth_dists, aes(x = pair, y = distance)) +
#   geom_boxplot()

p_strains <- blank_plot("Strain dendrograms\n(genetic/growth)")

#### Panel - Relative Sensitivity per strain ####
# Which strain has most sig knockouts per condition
condition_sensitivity <- filter(ko_growth, qvalue < 0.05) %>%
  select(strain, condition, name, score) %>%
  group_by(condition) %>%
  summarise(S288C = sum(!is.na(score) & strain == "S288C"),
            UWOP = sum(!is.na(score) & strain == "UWOP"),
            Y55 = sum(!is.na(score) & strain == "Y55"),
            YPS = sum(!is.na(score) & strain == "YPS")) %>%
  mutate(condition = factor(condition, levels = condition[order(.$S288C/(.$S288C + .$UWOP + .$Y55 + .$YPS))], ordered = TRUE)) %>%
  pivot_longer(names_to = 'strain', values_to = 'sig_count', S288C:YPS)

p_sensitivity <- ggplot(condition_sensitivity, aes(x=condition, y=sig_count, fill=strain)) +
  geom_col(position='fill') +
  scale_fill_brewer(name = "Strain", type = "qual", palette = "Set1") +
  scale_y_continuous(expand = expansion(0)) +
  labs(y = 'Relative conditionally significant gene count', x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major.y = element_blank(),
        axis.ticks.x = element_blank())
  
#### Panel - Mal gene example heatmap ####
# MAL genes noticed first
mal_ko <- filter(ko_growth, condition == "Maltose 2% (48H)", str_starts(name, "MAL")) %>%
  mutate(sig = ifelse(qvalue < 0.01, "*", ""))

p_mal <- ggplot(mal_ko, aes(x = strain, y = name, fill = score, label = sig)) +
  geom_raster() +
  geom_text() +
  scale_fill_gradientn(colours = c("#b10026", "#fc4e2a", "#feb24c", "#f7f7f7", "#4393c3", "#2d004b"),
                       values = scales::rescale(c(-7.5, -5, -2.5, 0, 2.5, 5)), name = "S-Score") +
  labs(subtitle = "Maltose 2% (48H)") +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank())

#### Panel - Other example heatmap ####
# A second set of example genes showing switching behavior
nacl_ko <- filter(ko_growth, condition == "NaCl 0.6M (72H)") %>%
  group_by(name) %>%
  filter(sum(qvalue < 0.01) > 1) %>%
  ungroup() %>%
  arrange(name, strain) %>%
  mutate(sig = ifelse(qvalue < 0.01, "*", ""))

p_nacl <- ggplot(nacl_ko, aes(x = strain, y = name, fill = score, label = sig)) +
  geom_raster() +
  geom_text() +
  scale_fill_gradientn(colours = c("#b10026", "#fc4e2a", "#feb24c", "#f7f7f7", "#4393c3", "#2d004b"),
                       values = scales::rescale(c(-7.5, -5, -2.5, 0, 2.5, 5)), name = "S-Score") +
  labs(subtitle = "NaCl 0.6M (72H)") +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.y = element_blank())

#### Panel - Proportion of shared phenotypes ####
double_prop <- mutate(ko_growth, sig = sign(score) * (qvalue < 0.01)) %>%
  select(strain, name, condition, sig) %>%
  pivot_wider(names_from = strain, values_from = sig) %>%
  mutate(S288C_UWOP = S288C == UWOP,
         S288C_Y55 = S288C == Y55,
         S288C_YPS = S288C == YPS,
         UWOP_Y55 = UWOP == Y55,
         UWOP_YPS = UWOP == YPS,
         YPS_Y55 = YPS == Y55) %>%
  select(S288C_UWOP:YPS_Y55) %>%
  pivot_longer(everything(), names_to = "pair", values_to = "shared") %>%
  group_by(pair) %>%
  summarise(prop = sum(shared, na.rm = TRUE) / sum(!is.na(shared)))

p_prop <- blank_plot("Proportion of shared phenotypes")

#### Panel - Condition pair correlation ####
# Condition pair correlations within strains
p_pair_cor <- blank_plot("Condition pair correlations")

#### Figure Assembly ####
size <- theme(text = element_text(size = 8))
p1 <- p_strains + labs(tag = 'A') + size
p2 <- p_sensitivity + labs(tag = 'B') + size
p3 <- p_mal + labs(tag = 'C') + size
p4 <- p_nacl + labs(tag = 'D') + size
p5 <- p_prop + labs(tag = 'E') + size
p6 <- p_pair_cor + labs(tag = 'F') + size

figure <- multi_panel_figure(width = 360, height = 240, columns = 3, rows = 2,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 1, column = 2) %>%
  fill_panel(p3, row = 1, column = 3) %>%
  fill_panel(p4, row = 2, column = 1) %>%
  fill_panel(p5, row = 2, column = 2) %>%
  fill_panel(p6, row = 2, column = 3)

ggsave('figures/thesis_figure_deletions.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
ggsave('figures/thesis_figure_deletions.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')

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

# Growth cor (comment out growth filtering)
# select(ko_growth, strain, gene, condition, score) %>%
#   extract(condition, c("condition", "time"), "([^\\(\\)]*) \\((48H|72H)\\)") %>%
#   drop_na() %>%
#   pivot_wider(names_from = time, values_from = score) %>%
#   drop_na() %>%
#   group_by(condition) %>%
#   group_modify(~broom::tidy(cor.test(.$`48H`, .$`72H`)))

ko_comparisons <- read_tsv("data/raw/ko_comparisons.tsv", skip = 2) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>%
  left_join(select(gene_meta, gene, name), by = "gene")

gene_exclusivity <- readxl::read_xlsx("data/raw/gene_exclusivity.xlsx", sheet = 2, col_names = c("gene", "exclusiveness"), skip = 1) %>%
  left_join(select(gene_meta, gene, name), by = "gene")

mash_dist <- tribble(
  ~strain, ~S288c,    ~UWOP,    ~Y55,     ~YPS,
  "S288C",  0.000000, 0.005539, 0.005905, 0.006047,
  "UWOP",   0.005539, 0.000000, 0.005615, 0.004640,
  "Y55",    0.005905, 0.005615, 0.000000, 0.006115,
  "YPS", 0.006047, 0.004640, 0.006115, 0.000000,
) %>%
  tblhelpr::tibble_to_matrix(-strain, row_names = "strain") %>%
  as.dist()

#### Panel - Strain Relationships - genetic and phenotypic ####
geno_hc <- hclust(mash_dist, method = "single")
geno_data <- dendro_data(geno_hc)

margin_size <- 15
p_strains <- ggplot() +
  geom_segment(data = geno_data$segments, aes(x = x, xend = xend, y = y, yend = ifelse(yend == 0, 0.004, yend))) +
  geom_point(data = geno_data$labels, aes(x = x, y = y + 0.004, colour = label), show.legend = FALSE, size = 3) +
  scale_x_continuous(breaks = 1:4, labels = geno_data$labels$label, sec.axis = dup_axis(name = "")) +
  scale_y_continuous(expand = expansion(0), sec.axis = dup_axis(name = "")) +
  scale_colour_brewer(name = "Strain", type = "qual", palette = "Set1") +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Mash Distance") +
  theme(axis.title.x.top = element_text(margin = margin(margin_size,0,0,0,unit = "mm")),
        axis.title.y.right = element_text(margin = margin(0,margin_size,0,0,unit = "mm")),
        axis.title.x.bottom = element_text(margin = margin(0,0,margin_size,0,unit = "mm")),
        axis.title.y.left = element_text(margin = margin(0,0,0,margin_size,unit = "mm")),
        axis.ticks = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text.y.right = element_blank())

#### Panel - Relative Sensitivity per strain ####
# Which strain has most sig knockouts per condition
condition_sensitivity <- filter(ko_growth, qvalue < 0.01) %>%
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
  scale_y_continuous(expand = expansion(c(0.01, 0.05))) +
  coord_flip(clip = "off") +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(y = 'Relative conditionally significant\ngene count', x = "") +
  theme(panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,-10,0))
  
#### Panel - Mal gene example heatmap ####
# MAL genes noticed first
mal_ko <- filter(ko_growth, condition == "Maltose 2% (48H)", str_starts(name, "MAL")) %>%
  mutate(sig = ifelse(qvalue < 0.01, "*", ""))

mal_margin <- 5
p_mal <- ggplot(mal_ko, aes(x = strain, y = name, fill = score, label = sig)) +
  geom_raster() +
  geom_text() +
  scale_fill_gradientn(colours = c("#b10026", "#fc4e2a", "#feb24c", "#f7f7f7", "#4393c3", "#053061"),
                       values = scales::rescale(c(-7.5, -5, -2.5, 0, 2.5, 5)), name = "S-Score") +
  labs(subtitle = "Maltose 2% (48H)", x = "", y = "") +
  theme(axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(margin = margin(mal_margin, 0, 0, 0, unit = "mm")),
        legend.box.margin = margin(0, mal_margin, 0, 0, unit = "mm"),
        axis.title.x = element_text(margin = margin(0, 0, mal_margin, 0, unit = "mm")),
        axis.title.y = element_text(margin = margin(0, 0, 0, mal_margin, unit = "mm")))

#### Panel - Other example heatmap ####
# A second set of example genes showing switching behavior
nacl_ko <- filter(ko_growth, condition == "NaCl 0.6M (72H)") %>%
  group_by(name) %>%
  filter(sum(qvalue < 0.01) > 1) %>%
  ungroup() %>%
  arrange(name, strain) %>%
  mutate(sig = ifelse(qvalue < 0.01, "*", ""))

nacl_margin <- 5
p_nacl <- ggplot(nacl_ko, aes(x = strain, y = name, fill = score, label = sig)) +
  geom_raster() +
  geom_text() +
  scale_fill_gradientn(colours = c("#b10026", "#fc4e2a", "#feb24c", "#f7f7f7", "#4393c3", "#053061"),
                       values = scales::rescale(c(-7.5, -5, -2.5, 0, 2.5, 5)), name = "S-Score") +
  labs(subtitle = "NaCl 0.6M (72H)", x = "", y = "") +
  theme(axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(margin = margin(nacl_margin, 0, 0, 0, unit = "mm")),
        legend.box.margin = margin(0, nacl_margin, 0, 0, unit = "mm"),
        axis.title.x = element_text(margin = margin(0, 0, nacl_margin, 0, unit = "mm")),
        axis.title.y = element_text(margin = margin(0, 0, 0, nacl_margin, unit = "mm")))

#### Panel - Proportion of shared phenotypes ####
calc_prop <- function(x) {
  x <- x[!is.na(x)]
  sum(x) / length(x)
}

get_props <- function(tbl, ...) {
  tbl <- select(tbl, condition, name, strain = strain2, shared) %>%
    pivot_wider(names_from = strain, values_from = shared)
  
  quad <- calc_prop(tbl[,3] & tbl[,4] & tbl[,5])
  triple <- c(calc_prop(tbl[,3] & tbl[,4]), calc_prop(tbl[,3] & tbl[,5]), calc_prop(tbl[,4] & tbl[,5]))
  double <- c(calc_prop(tbl[,3]), calc_prop(tbl[,4]) ,calc_prop(tbl[,5]))
  
  tibble(strains = 2:4, mean = c(mean(double), mean(triple), quad), sd = c(sd(double), sd(triple), NA))
}

props <- left_join(ko_comparisons, select(ko_growth, strain1 = strain, condition, gene, qvalue1 = qvalue),
          by = c("strain1", "condition", "gene")) %>%
  left_join(select(ko_growth, strain2 = strain, condition, gene, qvalue2 = qvalue),
            by = c("strain2", "condition", "gene")) %>%
  bind_rows(., set_names(., c("condition", "gene", "strain2", "strain1", "scores2", "scores1", "qvalue", "name", "qvalue2", "qvalue1"))) %>%
  filter((qvalue1 < 0.01 & qvalue < 0.01) | (qvalue1 < 0.01 & qvalue2 < 0.01 & qvalue >= 0.01)) %>%
  mutate(shared = qvalue >= 0.01) %>%
  group_by(strain = strain1) %>%
  group_modify(get_props)

p_prop <- ggplot(props, aes(x = strains, y = mean, ymin = mean - sd, ymax = mean + sd, fill = strain)) +
  geom_col(position = "dodge", width = 0.75) + 
  geom_errorbar(position = position_dodge(0.75), width = 0.25) +
  scale_fill_brewer(name = "", type = "qual", palette = "Set1") +
  scale_x_continuous(breaks = c(2, 3, 4)) +
  lims(y = c(0, 1)) +
  labs(x = "Number of Strains", y = "Proportion of Shared Significant Phenotypes")

#### Panel - Gene phenotypes ####
gene_exclusivity <- mutate(gene_exclusivity, rank = 1:n()) %>%
  drop_na()

top_genes <- filter(gene_exclusivity, rank <= 10) %>%
  mutate(x = 300, y = 0.833 - (rank - 1) * 0.05)

p_exclusiveness <- ggplot(gene_exclusivity, aes(x = rank, y = exclusiveness)) +
  geom_point(shape = 20) +
  geom_text(data = top_genes, mapping = aes(x = x, y = y, label = name), hjust = -0.05, vjust = 0.5) +
  geom_segment(data = top_genes, mapping = aes(x = x, y = y, xend = rank, yend = exclusiveness)) +
  labs(x = "Rank", y = "Exclusive phenotype propensity")

#### Figure Assembly ####
size <- theme(text = element_text(size = 12))
p1 <- p_strains + labs(tag = 'A') + size
p2 <- p_sensitivity + labs(tag = 'B') + size + theme(axis.text.y = element_text(size = 8))
p3 <- p_mal + labs(tag = 'C') + size
p4 <- p_nacl + labs(tag = 'D') + size
p5 <- p_prop + labs(tag = 'E') + size
p6 <- p_exclusiveness + labs(tag = 'F') + size

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

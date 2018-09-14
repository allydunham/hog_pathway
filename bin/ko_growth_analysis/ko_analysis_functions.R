#### Script containing the functions used for KO growth analysis

## Extract a matrix from ko_growth tibble for a strain giving a variable spread over gene or condition
# Does not apply sorting because ko_scores tsv comes sorted
split_strains <- function(str, var, tbl, col='gene', row='condition'){
  # Function to filter ko_growth to a given strain with a matrix of condition against 'var' (score or qvalue)
  ko <- filter(tbl, strain == str) %>%
    select(!!row, !!col, !!var) %>%
    spread(key = col, value = !!var)
  return(ko)
}

## calculate condition ko growth profile dendogram from a tibble of strain scores
get_dend <- function(x){
  mat <- as.matrix(select(x, -condition))
  rownames(mat) <- x$condition
  return(as.dendrogram(hclust(dist(mat))))
}

## Calculate correlation between pairs of genes or conditions (based on profile according to the other)
# requires a tibble with the dependant variable in rows labeled by column 'var'
get_cor <- function(x, cor_meth = 'pearson', cor_use = 'pairwise', upper_tri=TRUE, var='condition'){
  mat <- as.matrix(select(x, -!!var))
  rownames(mat) <- pull(x, !!var)
  cors <- cor(t(mat), method = cor_meth, use = cor_use)
  if (upper_tri){
    cors[lower.tri(cors, diag = TRUE)] <- NA
  }
  var1 <- paste0(var,'1')
  var2 <- paste0(var,'2')
  cors %<>% as_tibble(rownames=var1) %>%
    gather(key = !!var2, value = 'cor', -!!var1) %>%
    drop_na()
  return(cors)
}

# Determine which level list a gene is in
# sig_list should have the most desirable category last if there are overlaps
get_sig <- function(item, sig_list){
  for (i in length(sig_list):1){
    if (item %in% sig_list[[i]]){
      return(i)
    }
  }
  return(0)
}

# Determine genes significant in a condition
get_sig_genes <- function(con, growth_tbl, threshold = 0.01){
  gene <- growth_tbl %>%
    filter(condition == con, qvalue < threshold) %>%
    pull(name) %>%
    table()
  return(names(gene[gene > 1]))
}

## Generic function to give heatmaps for gene s-scores in a set of cons
plot_con_gene_heatmaps <- function(tbl, genes, cons=NULL, strains=NULL, primary_strain='S288C', facet_cols=2, facet_rows=2, sig_level=0.01){
  if (is.null(cons)){
    cons <- unique(tbl$condition)
  }
  if (is.null(strains)){
    strains <- unique(tbl$strain)
  }
  if (is.null(primary_strain) | !primary_strain %in% strains){
    primary_strain <- strains[1]
  }
  
  
  tbl <- filter(tbl, name %in% genes, condition %in% cons)
  sig <- mutate(tbl, sig = ifelse(qvalue < sig_level, '*', '')) %>%
    select(strain, name, condition, sig)
  
  tbl <- select(tbl, strain, condition, name, score) %>%
    spread(key = condition, value = score) %>%
    gather(key = 'condition', value = 'score', -strain, -name) %>%
    left_join(., sig, by = c('strain', 'name', 'condition')) %>%
    mutate(sig = ifelse(is.na(sig), '', sig))
  
  limits <- c(min(tbl$score, na.rm = TRUE), max(tbl$score, na.rm = TRUE))
  
  if (limits[2] - limits[1] < 10){
    clrs <- c('firebrick2', 'white', 'cornflowerblue')
    vals <- scales::rescale(c(limits[1], 0, limits[2]))
  } else {
    clrs <- c('firebrick2', 'darkgoldenrod1', 'white', 'cornflowerblue', 'darkorchid')
    vals <- scales::rescale(c(limits[1], limits[1]/2, 0, limits[2]/2, limits[2]))
  }
  
  #### Plot heatmap faceted by strain
  ## Determine ordering based on primary strain
  mat <- split_strains(str=primary_strain, tbl=tbl, row='name', col='condition', var='score') %>%
    tbl_to_matrix(., row = 'name') %>%
    set_na(., 0)
  
  gene_dend <- as.dendrogram(hclust(dist(mat)))
  con_dend <- as.dendrogram(hclust(dist(t(mat))))
  gene_order <- rownames(mat)[order.dendrogram(gene_dend)]
  con_order <- colnames(mat)[order.dendrogram(con_dend)]
  
  # Add any gene/cons not in the primary strain in a random order
  gene_order <- c(gene_order, genes[!genes %in% gene_order])
  con_order <- c(con_order, cons[!cons %in% con_order])
  
  tbl %<>% mutate(name = factor(name, levels = gene_order, ordered = TRUE),
                  condition = factor(condition, levels = con_order, ordered = TRUE))
  
  p_prim_gene_dend <- ggdendrogram(gene_dend, rotate=TRUE)
  p_prim_con_dend <- ggdendrogram(con_dend, rotate=TRUE)
  
  ## plot heatmap
  p_strain_heatmaps <- ggplot(tbl, aes(x=condition, y=name, fill=score)) +
    geom_raster() +
    geom_text(aes(label=sig)) +
    facet_wrap(~strain, ncol=facet_cols, nrow=facet_rows) +
    xlab('Condition') +
    ylab('Gene') +
    scale_fill_gradientn(colours = clrs, limits=limits, na.value = 'black', values = vals) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    guides(fill=guide_colourbar(title = 'S-Score'))
  
  #### Plot general combined heatmap
  ## Determine ordering
  mat <- tbl %>% select(-sig) %>%
    spread(key = condition, value = score) %>%
    unite(col='name', strain, name) %>%
    tbl_to_matrix(., row = 'name') %>%
    set_na(., 0)
  
  gene_dend <- as.dendrogram(hclust(dist(mat)))
  con_dend <- as.dendrogram(hclust(dist(t(mat))))
  gene_order <- rownames(mat)[order.dendrogram(gene_dend)]
  con_order <- colnames(mat)[order.dendrogram(con_dend)]
  
  # Add any gene/cons not in the primary strain in a random order
  gene_order <- c(gene_order, genes[!genes %in% gene_order])
  con_order <- c(con_order, cons[!cons %in% con_order])
  
  tbl %<>% unite(col='unit', strain, name, remove = FALSE) %>%
    mutate(unit = factor(unit, levels = gene_order, ordered = TRUE),
           condition = factor(condition, levels = con_order, ordered = TRUE)) 
  
  p_all_gene_dend <- ggdendrogram(gene_dend, rotate=TRUE, labels = TRUE)
  p_all_con_dend <- ggdendrogram(con_dend, labels = TRUE)
  
  ## Plot heatmap
  p_all_heatmap <- ggplot(tbl, aes(x=condition, y=unit, fill=score)) +
    geom_raster() +
    geom_text(aes(label=sig)) +
    xlab('Condition') +
    ylab('Gene') +
    scale_fill_gradientn(colours = clrs, na.value = 'black', values = vals) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    guides(fill=guide_colourbar(title = 'S-Score'))
  
  return(list(strain_heatmap=p_strain_heatmaps,
              strain_gene_dend=p_prim_gene_dend,
              strain_condition_dend=p_prim_con_dend,
              all_heatmap=p_all_heatmap,
              all_gene_dend=p_all_gene_dend,
              all_condition_dend=p_all_con_dend))
}

# Function to plot RNA seq expression levels relative to S288C in a gene set. Gene IDs take precedent over gene names
plot_rna_seq_fold_change <- function(tbl, gene_ids=NULL, gene_names=NULL, cons=NULL, strains=NULL, sig_level=0.01){
  if (!is.null(strains)){
    tbl <- filter(tbl, strain %in% strains)
  }
  if (!is.null(gene_ids)){
    tbl <- filter(tbl, gene %in% gene_ids)
    var = 'gene'
  } else if (!is.null(gene_names)){
    tbl <- filter(tbl, name %in% gene_names)
    var = 'name'
  } else {
    stop("One of 'genes' or 'names' must be provided")
  }
  
  tbl <- mutate(tbl, sig = ifelse(padj < sig_level, '*', ''))
  
  limits <- c(min(tbl$log2FoldChange, na.rm = TRUE), max(tbl$log2FoldChange, na.rm = TRUE))
  max_abs <- max(abs(limits))
  
  if (max_abs < 2){
    clrs <- c('firebrick2', 'white', 'cornflowerblue')
  } else {
    clrs <- c('firebrick2', 'darkgoldenrod1', 'white', 'cornflowerblue', 'darkorchid')
  }
  
  p <- ggplot(tbl, aes_string(y=var, x='strain', fill='log2FoldChange')) +
    geom_raster() +
    geom_text(aes(label=sig)) +
    xlab('Strain') +
    ylab('Gene') + 
    scale_fill_gradientn(colours = clrs, na.value = 'black', limits=c(-max_abs, max_abs)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    guides(fill=guide_colourbar(title = 'Log2 Fold Change'))
  
  return(p)
}



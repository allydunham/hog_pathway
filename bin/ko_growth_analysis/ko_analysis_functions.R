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

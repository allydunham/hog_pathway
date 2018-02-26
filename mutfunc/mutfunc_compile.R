## Script compiling full mutfunc Rds objects into tsv tables
setwd('/Users/ally/Projects/hog/mutfunc/')
require(data.table)

## Import orf/uniprot id conversions
uniprot2orf <- readRDS(file.path('mutfunc/uniprot2orf.rds'))
orf2uniprot <- names(uniprot2orf)
names(orf2uniprot) <- uniprot2orf

#### Import foldx data
fx_exp = fread( 'zcat < mutfunc/exp.tab.gz' )

# Add ORF names - plus filter to those in the uniprot/orf ID set
ind = fx_exp$uniprot_id %in% names(uniprot2orf)
fx_exp = fx_exp[ind,]
fx_exp$orf = uniprot2orf[ fx_exp$uniprot_id ]

## Add the same two ORF ids as in the aa muts table
fx_exp$pos_id = with(fx_exp, sprintf('%s %s', orf, uniprot_pos))
fx_exp$id = with(fx_exp, sprintf('%s %s%s%s', orf, aa_wt, uniprot_pos, aa_mt))

# Keep highest - remove any lower impact duplicates
fx_exp = fx_exp[order(ddG, decreasing = T),]
fx_exp = subset(fx_exp, !duplicated(id))

fx_exp$pdb_id = with(fx_exp, sprintf('%s_%s_%s_%s', toupper(pdb_id), chain, pdb_pos, aa_wt) )
fx_exp = fx_exp[,c('id', 'pos_id','pdb_id', 'ddG', 'ddG_sd'), with=F]
fx_exp$evidence = 'exp'

## Foldx modeled data
fx_mod = fread(paste('zcat <', 'mutfunc/mod.tab.gz'))

# Add orf names
ind = fx_mod$uniprot_id %in% names(uniprot2orf)
fx_mod = fx_mod[ind,]
fx_mod$orf = uniprot2orf[ fx_mod$uniprot_id ]

# Add mut id
fx_mod$id = with(fx_mod, sprintf('%s %s%s%s', orf, aa_wt, uniprot_pos, aa_mt))
fx_mod$pos_id = with(fx_mod, sprintf('%s %s', orf, uniprot_pos))

# Keep highest
fx_mod = fx_mod[order(ddG, decreasing = T),]
fx_mod = subset(fx_mod, !duplicated(id))

# Filter calls with experimental evidence
fx_mod = subset(fx_mod, !pos_id %in% fx_exp$pos_id)

fx_mod$pdb_id = sapply( strsplit(fx_mod$pdb_id, '_'), function(f) f[1] )
fx_mod$pdb_id = with(fx_mod, sprintf('%s_%s_%s_%s', toupper(pdb_id), chain, pdb_pos, aa_wt) )
fx_mod = fx_mod[,c('id', 'pos_id','pdb_id', 'ddG', 'ddG_sd'), with=F]

fx_mod$ddG = as.numeric(fx_mod$ddG)
fx_mod$ddG_sd = as.numeric(fx_mod$ddG_sd)
fx_mod$evidence = 'mod'

fx = as.data.frame( rbind(fx_exp, fx_mod) )

#### Foldx Interactions
fx_int = fread(paste('zcat <', 'mutfunc/int_clean.tab.gz'))

# Add orf names - filter by allowed names again
ind = fx_int$prot_a %in% names(uniprot2orf) &
  fx_int$prot_b %in% names(uniprot2orf)
fx_int = fx_int[ind,]

# Switch a/b - for interactions in the B chain perspective?
to_switch = fx_int$chain == 'B'
new_a = fx_int$prot_b[ to_switch ]
new_b = fx_int$prot_a[ to_switch ]
fx_int$prot_a[ to_switch ] = new_a
fx_int$prot_b[ to_switch ] = new_b

## Get Orf names
fx_int$orf_a = uniprot2orf[ fx_int$prot_a ]
fx_int$orf_b = uniprot2orf[ fx_int$prot_b ]

# Add mut id - same as before
fx_int$id = with(fx_int, sprintf('%s %s%s%s', orf_a, aa_wt, pdb_pos, aa_mt))
fx_int$pos_id = with(fx_int, sprintf('%s %s', orf_a, pdb_pos))

fx_int = fx_int[,c('id', 'pos_id','orf_a', 'orf_b', 'pdb','chain','evidence', 'ddG', 'ddG_sd'),with=F]

# Keep highest
fx_int = fx_int[order(ddG, decreasing = T)]
fx_int = subset(fx_int, !duplicated(id))

fx_int = as.data.frame(fx_int)

#### Elm data
elm = fread( file.path('mutfunc/elm_mut.tsv') )
elm$id = with(elm, sprintf('%s %s%s%s', orf, wt, pos, mt))
elm$pos_id = with(elm, sprintf('%s %s', orf, pos))

#### Other ptm data
ptms = readRDS(file.path('mutfunc/dbptm_yeast_parsed.rds'))
AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

ptms_dat = lapply(1:nrow(ptms), function(i){
  message(i)
  row = ptms[i,]
  cbind(mt=AA, row)
})

ptms_dat = rbindlist(ptms_dat)
ptms_dat = subset(ptms_dat, modified_residue != mt)

ptms_dat$id = with(ptms_dat, sprintf('%s %s%s%s', orf, modified_residue, position, mt))

#### Ps data
ps = readRDS(file.path('mutfunc/mimp_results_yeast_lr1_all_possible_muts.rds'))
ps$id = with(ps, sprintf('%s %s', gene, mut))

#### TF binding site muts
tfbs = readRDS(file.path('mutfunc/muts_chip_tfko_nofilter.rds'))


#### SIFT scores
cons = readRDS('mutfunc/sift_all_yeast.rds')
setDT(cons)
ind = cons$acc %in% names(uniprot2orf) # filter to desired muts
cons = cons[ind,]

cons$orf = uniprot2orf[ cons$acc ] # add orf ID
cons$pos_id = with(cons, sprintf('%s %s', orf, pos))

cons$id = with(cons, sprintf('%s %s%s%s', orf, ref, pos, alt))

#### Write full tsvs ####
write.table(fx,"../data/mutfunc/foldx.tsv",sep = '\t',quote = FALSE,row.names = FALSE)
write.table(fx_int,"../data/mutfunc/foldx-int.tsv",sep = '\t',quote = FALSE,row.names = FALSE)
write.table(cons,"../data/mutfunc/sift.tsv",sep = '\t',quote = FALSE,row.names = FALSE)
write.table(ps,"../data/mutfunc/pho.tsv",sep = '\t',quote = FALSE,row.names = FALSE)
write.table(ptms_dat,"../data/mutfunc/ptms.tsv",sep = '\t',quote = FALSE,row.names = FALSE)
write.table(elm,"../data/mutfunc/elms.tsv",sep = '\t',quote = FALSE,row.names = FALSE)
write.table(tfbs,"../data/mutfunc/tfbs.tsv",sep = '\t',quote = FALSE,row.names = FALSE)



  
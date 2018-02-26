#!/usr/bin/env Rscript
require(data.table)

args = commandArgs(trailingOnly=TRUE)

MUTFUNC_DIR = 'mutfunc'
## Load list of orf IDs by uniprot ID
## then reverse to generate uniprotIDs by orf ID
uniprot2orf = readRDS(file.path(MUTFUNC_DIR, 'uniprot2orf.rds'))
sp = split( names(uniprot2orf), uniprot2orf)
orf2uniprot = unlist(  lapply(sp, unname) )

## Load coding variants table (gives pos, codon change, aa change etc.)
df = read.table(args[1], sep='\t', header=T)
df$id =  with(df, sprintf('%s %s%s%s', gene, ref_aa, pos_aa, alt_aa))
df$pos_id = with(df, sprintf('%s %s', gene, pos_aa))

# genomic muts - list of genomic mutations in the form chrom:pos_ref/mut
genomic_muts = read.table(args[2], sep='\t', header=T)
genomic_muts = genomic_muts$mut_id

##################################
# foldx_pdb 
##################################
message('foldx_exp') # read in foldx data for mutants, with all possible aa changes and consequent sturctural impact
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

# Get minimal set - reduce to only those positions in your muts table
fx_exp_min = subset(fx_exp, pos_id %in% df$pos_id)
fx_exp_min$pdb_id = with(fx_exp_min, sprintf('%s_%s_%s_%s', toupper(pdb_id), chain, pdb_pos, aa_wt) )
fx_exp_min = fx_exp_min[,c('id', 'pos_id','pdb_id', 'ddG', 'ddG_sd'), with=F]
fx_exp_min$evidence = 'exp'

##################################
# foldx_mod
################################## - same process for fold x modeled impacts? vs experimental
message('foldx_mod')
mod_path = 'mutfunc/mod.tab.gz'
fx_mod = fread(paste('zcat <', mod_path))

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

# Get minimal set
fx_mod_min = subset(fx_mod, pos_id %in% df$pos_id)
fx_mod_min = subset(fx_mod_min, !pos_id %in% fx_exp$pos_id)

fx_mod_min$pdb_id = sapply( strsplit(fx_mod_min$pdb_id, '_'), function(f) f[1] )
fx_mod_min$pdb_id = with(fx_mod_min, sprintf('%s_%s_%s_%s', toupper(pdb_id), chain, pdb_pos, aa_wt) )
fx_mod_min = fx_mod_min[,c('id', 'pos_id','pdb_id', 'ddG', 'ddG_sd'), with=F]

fx_mod_min$ddG = as.numeric(fx_mod_min$ddG)
fx_mod_min$ddG_sd = as.numeric(fx_mod_min$ddG_sd)
fx_mod_min$evidence = 'mod'

fx_mod_exp = as.data.frame( rbind(fx_exp_min, fx_mod_min) )

write.table(fx_mod_exp, 'mutfunc_foldx_exp_mod.tsv', sep='\t', row.names=F, quote=F)

##################################
# foldx_int
################################## - process fold X interaction impacts
message('foldx_int')
int_path = 'mutfunc/int_clean.tab.gz'
fx_int = fread(paste('zcat <', int_path))

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

## subset those interactions of interest
fx_int_min = subset(fx_int, pos_id %in% df$pos_id)
fx_int_min = fx_int_min[,c('id', 'pos_id','orf_a', 'orf_b', 'pdb','chain','evidence', 'ddG', 'ddG_sd'),with=F]

fx_int_min = fx_int_min[order(ddG, decreasing = T)]
fx_int_min = subset(fx_int_min, !duplicated(id))

fx_int_min = as.data.frame(fx_int_min)
#rownames(fx_int_min) = fx_int_min$id
#fx_int_min = fx_int_min[,names(fx_int_min) != 'id']

write.table(fx_int_min, 'mutfunc_foldx_int.tsv', sep='\t', row.names=F, quote=F)

##################################
# elm
################################## - add consistent ids to elm (?) data
message('elm')
elm = fread( file.path(MUTFUNC_DIR, 'elm_mut.tsv') )
elm$id = with(elm, sprintf('%s %s%s%s', orf, wt, pos, mt))
elm$pos_id = with(elm, sprintf('%s %s', orf, pos))

elm = as.data.frame( subset(elm, pos_id %in% df$pos_id) )
# rownames(elm) = elm$id
# elm = elm[,!names(elm) %in% c('pos_id', 'id')]

#saveRDS(elm, 'mutfunc_elm.rds')
write.table(elm, 'mutfunc_elm.tsv', sep='\t', row.names=F, quote=F)

##################################
# ptms
##################################
message('ptms')
ptms = readRDS(file.path(MUTFUNC_DIR, 'dbptm_yeast_parsed.rds'))
AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

ptms_dat = lapply(1:nrow(ptms), function(i){
message(i)
row = ptms[i,]
cbind(mt=AA, row)
})

ptms_dat = rbindlist(ptms_dat)
ptms_dat = subset(ptms_dat, modified_residue != mt)

ptms_dat$id = with(ptms_dat, sprintf('%s %s%s%s', orf, modified_residue, position, mt))

## Generates a df of all possible mutations at PTM sites and their consistent ID

#saveRDS(ptms_dat, 'mutfunc_ptm.rds')
write.table(ptms_dat, 'mutfunc_ptm.tsv', sep='\t', row.names=F, quote=F)

##################################
# phos
################################## - add ids to phosphorylation mut sites from MIMP
message('phos')
ps = readRDS(file.path(MUTFUNC_DIR, 'mimp_results_yeast_lr1_all_possible_muts.rds'))
ps$id = with(ps, sprintf('%s %s', gene, mut))

#df_phos = merge(df, ps, by='id')
#saveRDS(ps, 'mutfunc_phos.rds')
write.table(ps, 'mutfunc_phos.tsv', sep='\t', row.names=F, quote=F)

##################################
# tfbs
################################## - subset TF binding site mutants to those considered
message('tfbs')
tfbs = readRDS(file.path(MUTFUNC_DIR, 'muts_chip_tfko_nofilter.rds'))
tfbs_muts = subset(tfbs, tfbs$mut_id %in% genomic_muts)
#saveRDS(tfbs_muts, 'mutfunc_tfbs.rds')
write.table(tfbs_muts, 'mutfunc_tfbs.tsv', sep='\t', row.names=F, quote=F)

##################################
# sift
##################################
message('sift')
cons = readRDS('mutfunc/sift_all_yeast.rds')
setDT(cons)
ind = cons$acc %in% names(uniprot2orf) # filter to desired muts
cons = cons[ind,]

cons$orf = uniprot2orf[ cons$acc ] # add orf ID
cons$pos_id = with(cons, sprintf('%s %s', orf, pos))

cons_min = as.data.frame( subset(cons, pos_id %in% df$pos_id) ) # sebset to only desired muts
cons_min$id = with(cons_min, sprintf('%s %s%s%s', orf, ref, pos, alt))

#saveRDS( cons_min[,c('id', 'pos_id', 'score')], 'mutfunc_cons.rds')
write.table(cons_min[,c('id', 'pos_id', 'score')], 'mutfunc_cons.tsv', sep='\t', row.names=F, quote=F)

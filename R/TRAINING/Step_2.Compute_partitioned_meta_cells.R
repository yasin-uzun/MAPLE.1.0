library(data.table)
library(GenomicRanges)
library(rdist)
library(reshape2)

source('met_matrix_processing.R')
source('classes.R')

#dataset = 'scmandt_muscle'
#samples = c('GSM')

#dataset = 'scmandt'
#samples = c('Serum')

#dataset = 'scnmt'
#sample_name = 'combined'

dataset = 'scnmt_gastr'
sample_name = 'combined'



run_parallel = T

part_rate = 0.2
#part_rate = 0


#bin_size = 200
bin_size = 500
num_neighbors = 20


data_dir = paste0('public_met/',dataset ,'/data/')
part_dir = paste0('public_met/',dataset ,'/data/part/')

dir.create(part_dir)

binned_met_dir = paste0('public_met/',dataset ,'/data/binned_met/')

if(run_parallel)
{
  library(doSNOW)
  cl <- makeCluster(36, outfile="", type = 'SOCK')
  registerDoSNOW(cl)
}



list_meta_objects = list()
meta_list_file = paste0(part_dir, '/meta_list.bin_size_', bin_size, '.',
                        sample_name, '.part_rate_', part_rate,'.num_nei_',num_neighbors,'.rds')

metric = 'met_count'
binned_met_shared_rds = paste0(binned_met_dir, '/binned_met.bin_size_',bin_size, '.shared.', metric, '.',sample_name,'.rds')
metric = 'demet_count'
binned_demet_shared_rds = paste0(binned_met_dir, '/binned_met.bin_size_',bin_size, '.shared.', metric, '.',sample_name,'.rds')


exp_shared_file = paste0(data_dir,'/expression_matrix/count_matrix.shared.',sample_name,'.rds')


df_exp = readRDS(exp_shared_file)
dim(df_exp)
df_exp[1:5, 1:5]


df_met_shared = readRDS(binned_met_shared_rds)
dim(df_met_shared)
df_met_shared[1:5, 1:5]

df_demet_shared = readRDS(binned_demet_shared_rds)
dim(df_demet_shared)
df_demet_shared[1:5, 1:5]

length(unique(df_met_shared$gene_symbol))

print(meta_list_file)
if(part_rate > 0)
{
  list_part_cells_file = paste0(part_dir, '/' , part_rate ,'/list_part_cells.',sample_name,'.rds')

  list_part_cells = readRDS(list_part_cells_file)

  all_cells = unique(df_met_shared$cell)


  for(i in 1:length(list_part_cells))
  {

    cells_leaveout = as.character(list_part_cells[[i]] )
    cells_keep = setdiff(all_cells, cells_leaveout)


    df_met_part_keep     = df_met_shared[df_met_shared$cell %in% cells_keep, ]
    df_demet_part_keep   = df_demet_shared[df_demet_shared$cell %in% cells_keep, ]

    df_met_part_leaveout     = df_met_shared[df_met_shared$cell %in% cells_leaveout, ]
    df_demet_part_leaveout   = df_demet_shared[df_demet_shared$cell %in% cells_leaveout, ]


    df_exp_part_keep = df_exp[, cells_keep]
    df_exp_part_leaveout = df_exp[, cells_leaveout]

    print(df_exp_part_keep[1:5,       1:5])
    print(df_exp_part_leaveout[1:5,       1:5])
    
    dim(df_met_part_keep)
    dim(df_demet_part_keep)

    dim(df_met_part_leaveout)
    dim(df_demet_part_leaveout)

    meta_keep = compute_meta_cells(df_met_shared = df_met_part_keep, df_demet_shared = df_demet_part_keep, df_exp = df_exp_part_keep, exp_aggr_type = 'sum')
    meta_leaveout = compute_meta_cells(df_met_shared = df_met_part_leaveout, df_demet_shared = df_demet_part_leaveout, df_exp = df_exp_part_leaveout, 
                                       exp_aggr_type = 'sum')

    list_meta_objects[[i]] = list('meta_keep' = meta_keep, 'meta_leaveout' = meta_leaveout)
  }


}else
{

  meta_object = compute_meta_cells(df_met_shared = df_met_shared, df_demet_shared = df_demet_shared, df_exp = df_exp, exp_aggr_type = 'sum')

  list_meta_objects[[1]] = list('meta_keep' = meta_object)

}

saveRDS(list_meta_objects, file = meta_list_file)
print(meta_list_file)



stopCluster(cl)


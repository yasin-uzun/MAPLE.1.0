library(data.table)
library(GenomicRanges)
library(rdist)
library(reshape2)

source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/met_matrix_processing.v11.R')
source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/classes.v01.R')

#dataset = 'scmandt_muscle'
#samples = c('GSM')

#dataset = 'scmandt'
#samples = c('Serum')

dataset = 'scnmt'
#samples = c('EB')
#samples = c('combined')
sample_name = 'combined'

dataset = 'sctrio'
#samples = c('EB')
#samples = c('combined')
sample_name = 'CRC01'

dataset = 'scnmt_gastr'
#sample_names = c('combined')
sample_name = 'combined'


dataset = 'snmc_mouse'
sample_name = 'combined'

dataset = 'snmc_human'
sample_name = 'combined'

run_parallel = F

#part_rate = 0.2 
part_rate = 0


#bin_size = 200
bin_size = 500
num_neighbors = 20
distance_dim_count = 10


data_dir = paste0('/mnt/isilon/tan_lab/uzuny/projects/public_met/',dataset ,'/data/')
part_dir = paste0('/mnt/isilon/tan_lab/uzuny/projects/public_met/',dataset ,'/data/part/')



binned_met_dir = paste0('/mnt/isilon/tan_lab/uzuny/projects/public_met/',dataset ,'/data/binned_met/')

if(run_parallel)
{
library(doSNOW)
cl <- makeCluster(20, outfile="", type = 'SOCK')
registerDoSNOW(cl)
}



#/mnt/isilon/tan_lab/uzuny/projects/public_met/snmc_mouse/data/binned_met//binned_demet.bin_size_500.demet_count.combined.rds
#/mnt/isilon/tan_lab/uzuny/projects/public_met/snmc_mouse/data/binned_met//binned_met.bin_size_500.demet_count.combined.rds

#for(sample_name in samples)
#{
  
  list_meta_objects = list()
  meta_list_file = paste0(part_dir, '/meta_list.bin_size_', bin_size, '.',sample_name, '.part_rate_', part_rate,'.num_nei_',num_neighbors,'.rds')
  
  metric = 'met_count'
  binned_met_shared_rds = paste0(binned_met_dir, '/binned_met.bin_size_',bin_size, '.shared.', metric, '.',sample_name,'.rds')
  metric = 'demet_count'
  binned_demet_shared_rds = paste0(binned_met_dir, '/binned_met.bin_size_',bin_size, '.shared.', metric, '.',sample_name,'.rds')
  
  if(dataset == 'sctrio')
  {
    exp_shared_file = paste0(data_dir,'/expression_matrix/expr_matrix.shared.',sample_name,'.rds')
    
  }else if(dataset == 'snmc_human' || dataset == 'snmc_mouse' )
  {
    exp_shared_file = paste0(data_dir,'/expression_matrix/counts_matrix.full_coding.rds')
    metric = 'met_count'
    binned_met_shared_rds = paste0(binned_met_dir, '/binned_met.bin_size_',bin_size, '.', metric, '.',sample_name,'.rds')
    metric = 'demet_count'
    binned_demet_shared_rds = paste0(binned_met_dir, '/binned_demet.bin_size_',bin_size, '.', metric, '.',sample_name,'.rds')
    
  }else
  {
    exp_shared_file = paste0(data_dir,'/expression_matrix/count_matrix.shared.',sample_name,'.rds')
  }
  
  df_exp = readRDS(exp_shared_file)
  dim(df_exp)
  df_exp[1:5, 1:5]

  #df_met_shared = df_binned_met
  #df_demet_shared = df_binned_demet
  
  
  df_met_shared = readRDS(binned_met_shared_rds)
  dim(df_met_shared)
  df_met_shared[1:5, 1:5]
  
  df_demet_shared = readRDS(binned_demet_shared_rds)
  dim(df_demet_shared)
  df_demet_shared[1:5, 1:5]
  
  length(unique(df_met_shared$gene_symbol))
  
  
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
      
      
      dim(df_met_part_keep)
      dim(df_demet_part_keep)
      
      dim(df_met_part_leaveout)
      dim(df_demet_part_leaveout)
      
      meta_keep = compute_meta_cells(df_met_shared = df_met_part_keep, df_demet_shared = df_demet_part_keep, df_exp = df_exp_part_keep, exp_aggr_type = 'sum')
      meta_leaveout = compute_meta_cells(df_met_shared = df_met_part_leaveout, df_demet_shared = df_demet_part_leaveout, df_exp = df_exp_part_leaveout, exp_aggr_type = 'sum')
      
      list_meta_objects[[i]] = list('meta_keep' = meta_keep, 'meta_leaveout' = meta_leaveout)
    }
    
  
  }else
  {
    
    if(dataset == 'sctrio')
    {
      meta_object = compute_meta_cells(df_met_shared = df_met_shared, df_demet_shared = df_demet_shared, df_exp = df_exp, exp_aggr_type = 'mean')
      
    }else if(dataset == 'snmc_mouse' | dataset == 'snmc_humans'){
      meta_object = compute_meta_cells(df_met_shared = df_met_shared, df_demet_shared = df_demet_shared, df_exp = NULL, exp_aggr_type = 'sum', run_parallel = F)
    }else
    {
      meta_object = compute_meta_cells(df_met_shared = df_met_shared, df_demet_shared = df_demet_shared, df_exp = df_exp, exp_aggr_type = 'sum')
    }
    
    list_meta_objects[[1]] = list('meta_keep' = meta_object)
    
  }
  
  saveRDS(list_meta_objects, file = meta_list_file)
  

  
#}


stopCluster(cl)


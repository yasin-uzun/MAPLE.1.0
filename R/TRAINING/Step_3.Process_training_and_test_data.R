library(data.table)
#library(keras)
#library(glmnet)

library(data.table)
library(GenomicRanges)

library(reshape2)

source('classes.R')
source('met_matrix_processing.R')
source('features_and_response_processing.R')

organism = 'mouse'
norm_scope = 'matrix'

dataset = 'scmandt_muscle'
sample_name = 'GSM'

#dataset = 'scmandt'
#sample_name = 'Serum'

#dataset = 'scnmt'
#sample_name = 'EB'


datasets = c('scmandt_muscle', 'scmandt', 'scnmt')
sample_names = c('GSM', 'Serum', 'EB')

datasets = c('scnmt')
sample_names = c('combined')




dataset = 'scnmt'
sample_name = 'combined'

dataset = 'scnmt_gastr'
sample_names = c('combined')
sample_name = 'combined'


bin_size = 500
#num_neighbors = 20
num_neighbors = 5
distance_dim_count = 10


if(organism == 'human'){
  cpg_content_file = 'human/hg38/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'
}

if(organism == 'mouse'){
  cpg_content_file = 'mouse/mm10/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'
}


df_cpg_content = readRDS(cpg_content_file)
head(df_cpg_content)
dim(df_cpg_content)
#df_cpg_content['H2BFS', ]



#part_object_file = paste0(part_dir, '/' , part_rate , '/obj_part.rds')

#library(doSNOW)
#cl <- makeCluster(20, outfile="", type = 'SOCK')
#registerDoSNOW(cl)

#part_rate = 0.2
part_rate = 0

datasets = c('scmandt_muscle', 'scmandt', 'scnmt', 'scnmt_gastr')
sample_names = c('GSM', 'Serum', 'combined', 'combined')


num_neighbors = 20
max_na_bins = 0
max_na_bins = 19


for(i in 1:length(datasets))
{


  dataset = datasets[i]
  sample_name = sample_names[i]


  print(dataset)

  data_dir = paste0('public_met/',dataset ,'/data/')
  part_dir = paste0('public_met/',dataset ,'/data/part/')


  meta_list_file = paste0(part_dir, '/meta_list.bin_size_', bin_size, '.',
                          sample_name, '.part_rate_', part_rate,
                          '.num_nei_',num_neighbors,'.rds')

 

  meta_list = readRDS(meta_list_file)


  fr_lists = list()
  

  if(part_rate == 0)
  {
    part_id = 1
    meta_part = meta_list[[part_id]]
    dim(meta_part$meta_keep@df_rate_meta_mean)
    meta_part$meta_keep@df_rate_meta_mean[1:5, 1:5]
    dim(meta_part$meta_keep@df_met_counts)
    
    part_id = 1
    meta_part = meta_list[[part_id]]
    
    fr_list_entire = get_fr_list(meta_data = meta_part$meta_keep, max_na_bins = 0, norm_scope = norm_scope)
    fr_lists[[part_id]] = list(fr_list_entire = fr_list_entire)
    
  
    head(fr_list_entire$features_mean)
    dim(fr_list_entire$features_matrix)
    fr_list_entire$features_matrix[1:5, 1:5]
    dim(fr_list_entire$features_array)
    head(fr_list_entire$response)
    length(fr_list_entire$response)
    sum(!is.na(fr_list_entire$response))
    sum(is.na(fr_list_entire$response))
    
    fr_lists[[part_id]] = list(fr_list_entire = fr_list_entire)

  }else
  {

    for(part_id in 1:5)
    {
      cat('Processing partition ', part_id, ' ... \n ')
      meta_part = meta_list[[part_id]]

      colnames(meta_part$meta_keep@df_exp_meta) = gsub('E4.5.5.5', 'E4.5-5.5', colnames(meta_part$meta_keep@df_exp_meta))
      colnames(meta_part$meta_leaveout@df_exp_meta) = gsub('E4.5.5.5', 'E4.5-5.5', colnames(meta_part$meta_leaveout@df_exp_meta))
      
      fr_list_train = get_fr_list(meta_data = meta_part$meta_keep, max_na_bins = 0, norm_scope = norm_scope)
      fr_list_test = get_fr_list(meta_data = meta_part$meta_leaveout, max_na_bins = 0, norm_scope = norm_scope)
      #fr_list_test = get_fr_list(meta_data = meta_part$meta_leaveout, max_na_bins = 5, norm_scope = norm_scope)
      print(head(fr_list_test$response))
      fr_lists[[part_id]] = list(fr_list_train = fr_list_train, fr_list_test = fr_list_test)
      
      dim(meta_part$meta_keep@df_exp_meta)
      dim(meta_part$meta_leaveout@df_exp_meta)
      
    }#for(part_id in 1:5)

  }#else

  fr_lists_file = paste0(part_dir, '/fr_lists.bin_size_',
                         bin_size, '.', sample_name,
                         '.part_rate_', part_rate,
                         '.norm_scope_', norm_scope,
                         '.max_na_bins_',max_na_bins,
                         '.num_nei_',num_neighbors,'.rds')
  
  #print(fr_lists_file)
  
  saveRDS(fr_lists, fr_lists_file)
  cat("fr_list_file: ", fr_lists_file, "\n")

}




temp0 = fr_lists[[1]]
head(temp0$fr_list_entire$response)

temp1 = fr_list_test
head(temp1$response)
length(temp1$response)


temp2 = fr_list_train
head(temp2$response)
length(temp2$response)


  


#stopCluster(cl)


 
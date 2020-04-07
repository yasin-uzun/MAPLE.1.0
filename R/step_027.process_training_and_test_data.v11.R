library(data.table)
#library(keras)
#library(glmnet)

library(data.table)
library(GenomicRanges)

library(reshape2)

source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/classes.v01.R')
source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/met_matrix_processing.v13.R')
source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/features_and_response_processing.v20.R')
#source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/training.v06.R')


organism = 'mouse'


#cpg_content_file = '/mnt/isilon/tan_lab/uzuny/genomes/human/hg38/annotated_regions//regions.genes.tss_ud_5K_cpg_ratio.bin_size_500.rds'
#cpg_content_file = '/mnt/isilon/tan_lab/uzuny/genomes/mouse/mm10/annotated_regions/regions.genes.tss_ud_5K.cpg_ratio.bin_size_500.rds'
#cpg_content_file = '/mnt/isilon/tan_lab/uzuny/genomes/mouse/mm10/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'


#dataset = 'scmandt_muscle'
#sample_name = 'GSM'

#dataset = 'scmandt'
#sample_name = 'Serum'

#dataset = 'scnmt'
#sample_name = 'EB'


datasets = c('scmandt_muscle', 'scmandt', 'scnmt')
sample_names = c('GSM', 'Serum', 'EB')

datasets = c('scnmt')
sample_names = c('combined')

norm_scope = 'matrix'
#norm_scope = 'gene'


dataset = 'scnmt'
sample_name = 'combined'



dataset = 'scnmt_gastr'
sample_names = c('combined')
sample_name = 'combined'

dataset = 'sctrio'
sample_name = 'CRC01'
#organism = 'human'

dataset = 'snmc_mouse'
sample_name = 'combined'

dataset = 'snmc_human'
sample_name = 'combined'
organism = 'human'


bin_size = 500
num_neighbors = 20
distance_dim_count = 10


if(organism == 'human'){
  cpg_content_file = '/mnt/isilon/tan_lab/uzuny/genomes/human/hg38/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'
}

if(organism == 'mouse'){
  cpg_content_file = '/mnt/isilon/tan_lab/uzuny/genomes/mouse/mm10/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'
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


#for(i in 1:length(datasets))
#{
  
  
  #dataset = datasets[i]
  #sample_name = sample_names[i]
  
  
  print(dataset)
  
  data_dir = paste0('/mnt/isilon/tan_lab/uzuny/projects/public_met/',dataset ,'/data/')
  part_dir = paste0('/mnt/isilon/tan_lab/uzuny/projects/public_met/',dataset ,'/data/part/')
  
  
  meta_list_file = paste0(part_dir, '/meta_list.bin_size_', bin_size, '.',
                          sample_name, '.part_rate_', part_rate,
                          '.num_nei_',num_neighbors,'.rds')
  
  fr_lists_file = paste0(part_dir, '/fr_lists.bin_size_',
                         bin_size, '.', sample_name,
                         '.part_rate_', part_rate,
                         '.norm_scope_', norm_scope,
                         '.num_nei_',num_neighbors,'.rds')
  
  
  meta_list = readRDS(meta_list_file)
  
  
  fr_lists = list()
  
  if(part_rate == 0)
  {
    part_id = 1
    meta_part = meta_list[[part_id]]
    
    fr_list_entire = get_fr_list(meta_data = meta_part$meta_keep, max_na_bins = 0, norm_scope = norm_scope)
    fr_lists[[part_id]] = list(fr_list_entire = fr_list_entire)
    
  }else
  {
    
    for(part_id in 1:5)
    {
      cat('Processing partition ', part_id, ' ... \n ')
      meta_part = meta_list[[part_id]]
      
      fr_list_train = get_fr_list(meta_part$meta_keep, max_na_bins = 0, norm_scope = norm_scope)
      fr_list_test = get_fr_list(meta_part$meta_leaveout, max_na_bins = 5, norm_scope = norm_scope)
      
      fr_lists[[part_id]] = list(fr_list_train = fr_list_train, fr_list_test = fr_list_train)
      
    }#for(part_id in 1:5)
    
  }#else
  
  
  saveRDS(fr_lists, fr_lists_file)


#}


#stopCluster(cl)


  
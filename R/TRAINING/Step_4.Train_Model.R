

library(reticulate)
py_config()


library(data.table)
library(keras)
library(glmnet)

library(data.table)
library(GenomicRanges)

library(reshape2)


#library(rsample)      # data splitting 
library(randomForest) # basic implementation
#library(ranger)       # a faster implementation of randomForest
#library(caret)        # an aggregator package for performing many machine learning models
#library(h2o)          # an extremely fast java-based platform


num_cpus = 20
library(doSNOW)
cl <- makeCluster(num_cpus, outfile="", type = 'SOCK')
registerDoSNOW(cl)



norm_scope = 'matrix'
#norm_scope = 'gene'



source('classes.R')
source('met_matrix_processing.R')
source('features_and_response_processing.R')
source('training.R')



bin_size = 500
num_neighbors = 20
distance_dim_count = 10

#train_subsample_size = 10000
train_subsample_size = 100000

flag_use_subsampled_train_set = T

# Directory to store the output
model_dir = 'models/'


#datasets = c('scmandt_muscle', 'scmandt', 'scnmt',  'scnmt')
#sample_names = c('GSM', 'Serum', 'EB', 'combined')

#datasets = c('scnmt')
#sample_names = c('combined')

#dataset = 'scnmt_gastr'
#sample_name = 'combined'

#datasets = c('scnmt_gastr')
#sample_names = c('combined')

dataset = 'scnmt'
sample_name = 'combined'


organism = 'mouse'



if(organism == 'mouse')
{
   cpg_content_file = 'mm10/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'
}else
{
  cpg_content_file = 'hg38/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'
}

df_cpg_content = readRDS(cpg_content_file)
head(df_cpg_content)
dim(df_cpg_content)




for(i in 1:length(datasets))
{
 
  dataset = datasets[i]
  sample_name = sample_names[i]
   
  # Input directory
  data_dir = paste0('public_met/',dataset ,'/data/')
  
  
  # For data partitioning (train/test)
  part_dir = paste0('public_met/',dataset ,'/data/part/')
  
  
  #part_rate = 0.2 
  part_rate = 0 
  part_id = 1
  
  
  #Training data###
  
  fr_lists_file = paste0(part_dir, '/fr_lists.bin_size_',   bin_size, '.', 
                          sample_name,
                         '.part_rate_', part_rate,
                         '.norm_scope_', norm_scope,
                         '.num_nei_',num_neighbors,'.rds')
  
  

  fr_lists = readRDS(fr_lists_file)
  
  
  fr_part = fr_lists[[part_id]]
  
  
  fr_list_train = fr_part$fr_list_entire
  

  if(flag_use_subsampled_train_set)
  {
    fr_list_train_to_be_used <- subsample_fr(fr_list_train, subsample_size = train_subsample_size, seed = 123)
    dim(fr_list_train_to_be_used$features_matrix)
    
  }else
  {
    fr_list_train_to_be_used <- fr_list_train
    train_subsample_size = 0
  }
  
  
  cnn_model_dir = paste0(model_dir, '/models/',dataset,'/cnn_model/')
  dir.create(cnn_model_dir, showWarnings = T, recursive = T)
  
  

  # ####Train CNN###

  patience = 10
  dropout_rate = 0.2
  hidden_act = 'relu'
  output_act = 'linear'
  loss = 'mse'

  training_result = train_cnn_model(fr_list_train_to_be_used$features_array,
                                   fr_list_train_to_be_used$response,
                                   patience = patience,
                                   hidden_act = hidden_act,
                                   output_act = output_act,
                                   loss = loss,
                                   
                                   dropout_rate = dropout_rate
                                   )
  
  cnn_model = training_result$model
  print('')



  cnn_model_file = paste0(cnn_model_dir,
                          'cnn_model.patience_',  patience,
                          '.hidden_act_', hidden_act,
                          '.output_act_',output_act,
                          '.loss_', loss,
                          '.tr_size_', train_subsample_size,
                          '.norm_scope_', norm_scope,
                          '.sample_name_', sample_name, 
                          
                          '.hd5')

  cnn_model %>% save_model_hdf5(cnn_model_file)


  
  
  ####Train Elastic###
  alpha = 0.5
  elastic_model = cv.glmnet(as.matrix(fr_list_train_to_be_used$features_matrix), as.numeric(fr_list_train_to_be_used$response), alpha = alpha)
  
  elastic_model_dir = paste0(model_dir, '/models/',dataset,'/elastic_model/')
  dir.create(elastic_model_dir, showWarnings = T, recursive = T)
  
  
  elastic_model_file = paste0(elastic_model_dir, 
                              'elastic_model.alpha_', alpha , 
                              '.tr_size_', train_subsample_size, 
                              '.norm_scope_', norm_scope, 
                              '.sample_name_', sample_name, 
                              '.rds')
  saveRDS(elastic_model, file = elastic_model_file)
  
  ####Train RF###
  
  
  sampsize = 0.8
  nodesize = 1
   
  rf_model <- foreach(ntree=rep(500/num_cpus, num_cpus), .combine=randomForest::combine,
                          .multicombine=TRUE, .packages='randomForest') %dopar% {
                            randomForest(y = as.numeric(fr_list_train_to_be_used$response),
                                         x = as.matrix(fr_list_train_to_be_used$features_matrix), 
                                         sampsize = sampsize * nrow(fr_list_train_to_be_used$features_matrix),
                                         nodesize = nodesize,
                                         ntree=ntree)
                          }
  
  rf_model_dir = paste0(model_dir, '/models/',dataset,'/rf_model/')
  dir.create(rf_model_dir, showWarnings = T, recursive = T)
  
  rf_model_file = paste0(rf_model_dir, 'rf_model.sampsize_', sampsize ,
                         '.nodesize_', nodesize , 
                         '.tr_size_', train_subsample_size,
                         '.norm_scope_', norm_scope,
                         '.sample_name_', sample_name, 
                         '.rds')
  
  saveRDS(rf_model, file = rf_model_file)
  

}






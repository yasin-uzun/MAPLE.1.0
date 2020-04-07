rm(list = ls())

# 
library(reticulate)
use_python("/opt/anaconda3/bin/python")
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

#stopCluster(cl)
library(doSNOW)
cl <- makeCluster(18, outfile="", type = 'SOCK')
registerDoSNOW(cl)


source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/classes.v01.R')
source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/met_matrix_processing.v05.R')
source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/features_and_response_processing.v20.R')
source('/mnt/isilon/tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/v18/generic_scripts/training.v06.R')

organism = 'mouse'
#





norm_scope = 'matrix'
#norm_scope = 'gene'


#num_cpus = 20
#library(doSNOW)
#cl <- makeCluster(num_cpus, outfile="", type = 'SOCK')
#registerDoSNOW(cl)


bin_size = 500
num_neighbors = 20
distance_dim_count = 10

#train_subsample_size = 10000
#test_subsample_size = 10000

train_subsample_size = 100000
test_subsample_size = 100000
#test_subsample_size = 0

flag_use_subsampled_test_set = F
#flag_use_subsampled_test_set = T

#test_dataset = 'scmandt_muscle'
#test_sample_name = 'GSM'

#test_dataset = 'scmandt'
#test_sample_name = 'Serum'

#test_dataset = 'scnmt'
#test_sample_name = 'EB'

package_data_dir = '/mnt/isilon//tan_lab/uzuny/scripts/dna_met/gene_expression_predictor/package_data/'

#datasets = c('scmandt_muscle', 'scmandt', 'scnmt',  'scnmt', 'scnmt_gastr')
#sample_names = c('GSM', 'Serum', 'EB', 'combined', 'combined')

datasets = c('scmandt_muscle', 'scmandt', 'scnmt',  'scnmt_gastr', 'snmc_mouse', 'snmc_human')
sample_names = c('GSM', 'Serum', 'combined', 'combined', 'combined', 'combined')

#datasets = c('snmc_human')
#sample_names = c('combined')


cor_mat_cnn = cor_mat_elastic = cor_mat_rf =  cor_mat_ensem = matrix(nrow = length(datasets), ncol = length(datasets), data = NA)
error_mat_cnn = error_mat_elastic = error_mat_rf = error_mat_ensem = matrix(nrow = length(datasets), ncol = length(datasets), data = NA)

cor_vec_mean = error_vec_mean = c()
cor_vec_null = error_vec_null = c()

rownames(cor_mat_cnn) = rownames(cor_mat_elastic) = rownames(cor_mat_rf) = rownames(cor_mat_ensem) = paste(datasets, sample_names)
colnames(cor_mat_cnn) = colnames(cor_mat_elastic) = colnames(cor_mat_rf) = colnames(cor_mat_ensem) = paste(datasets, sample_names)

rownames(error_mat_cnn) = rownames(error_mat_elastic) = rownames(error_mat_rf) = rownames(error_mat_ensem) = paste(datasets, sample_names)
colnames(error_mat_cnn) = colnames(error_mat_elastic) = colnames(error_mat_rf) = colnames(error_mat_ensem) = paste(datasets, sample_names)


list_response_test = list()
list_predict_mean = list()
list_predict_null = list()
list_predict_cnn = list()
list_predict_elastic = list()
list_predict_rf = list()
list_predict_ensem = list()

#i=5 #snmc_mouse
i=6 #snmc_human

#for(i in 1:length(datasets))
#{


#part_rate = 0.2 
part_rate = 0 
part_id = 1


test_dataset = datasets[i]
test_sample_name = sample_names[i]

cat('test dataset ', test_dataset, '\n')

if(test_dataset == 'snmc_human'){  organism = 'human'  }


if(organism == 'human'){
  cpg_content_file = '/mnt/isilon/tan_lab/uzuny/genomes/human/hg38/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'
}

if(organism == 'mouse'){
  cpg_content_file = '/mnt/isilon/tan_lab/uzuny/genomes/mouse/mm10/annotated_regions/regions.genes.tss_ud_5K.protein_coding.cpg_ratio.bin_size_500.rds'
}


df_cpg_content = readRDS(cpg_content_file)
head(df_cpg_content)
dim(df_cpg_content)

data_dir = paste0('/mnt/isilon/tan_lab/uzuny/projects/public_met/',test_dataset ,'/data/')
part_dir = paste0(data_dir ,'/part/')

#Test data###

fr_lists_file = paste0(part_dir, '/fr_lists.bin_size_',
                       bin_size, '.', test_sample_name,
                       '.part_rate_', part_rate,
                       '.norm_scope_', norm_scope,
                       '.num_nei_',num_neighbors,'.rds')




fr_lists = readRDS(fr_lists_file)

fr_part = fr_lists[[part_id]]

fr_list_test = fr_part$fr_list_entire


if(flag_use_subsampled_test_set)
{
  fr_list_test_to_be_used <- subsample_fr(fr_list_test, subsample_size = test_subsample_size, seed = 456)
}else
{
  fr_list_test_to_be_used <- fr_list_test
  #train_subsample_size = 0
  test_subsample_size = 0
}

sub_mat_list = list()
sub_array_list = list()

num_elements = nrow(fr_list_test_to_be_used$features_matrix)

length(fr_list_test_to_be_used$response)
sub_mat_size = 1000000
  
if(num_elements > sub_mat_size)
{
  sub_mat_count = ceiling(num_elements / sub_mat_size)
  for(sm in 1:sub_mat_count)
  {
    print(sm)
    start = (sm-1)*sub_mat_size + 1
    end = min(sm*sub_mat_size,  num_elements)
    sub_mat_list[[sm]]  = fr_list_test_to_be_used$features_matrix[start:end, ]
    sub_array_list[[sm]]  = fr_list_test_to_be_used$features_array[start:end, ,]
    
  }
}








list_response_test[[paste(test_dataset, test_sample_name)]] = fr_list_test_to_be_used$response


flag_true_exp_exists = length(fr_list_test_to_be_used$response) > 0 

#Mean Model                                        
predict_mean = 1 - fr_list_test_to_be_used$features_mean
names(predict_mean) = rownames(fr_list_test_to_be_used$features_matrix)

head(predict_mean)
length(predict_mean)

sum(is.na(predict_mean))

if( flag_true_exp_exists )
{
  spear_mean = cor(predict_mean, fr_list_test_to_be_used$response, method = 'spearman')
  spear_mean = round(spear_mean, 3)
  
  error_mean = mean(abs(predict_mean - fr_list_test_to_be_used$response))
  
  print(paste('Global test correlation (Mean):', spear_mean) )
  
  cor_vec_mean[paste(test_dataset, test_sample_name)] = spear_mean
  error_vec_mean[paste(test_dataset, test_sample_name)] = error_mean
  
  list_predict_mean[[paste(test_dataset, test_sample_name)]] = predict_mean
  
  
  #saveRDS(list_predict_mean, file = 'list_predict_mean.rds')
  
  #Null Model
  set.seed(123)
  predict_null = runif(length(fr_list_test_to_be_used$response))
  #predict_null = rnorm(length(fr_list_test_to_be_used$response))
  
  names(predict_null) = rownames(fr_list_test_to_be_used$features_matrix)
  
  spear_null = cor(predict_null, fr_list_test_to_be_used$response, method = 'spearman')
  spear_null = round(spear_null, 3)
  
  error_null = mean(abs(predict_null - fr_list_test_to_be_used$response))
  
  print(paste('Global test correlation (null):', spear_null) )
  
  cor_vec_null[paste(test_dataset, test_sample_name)] = spear_null
  error_vec_null[paste(test_dataset, test_sample_name)] = error_null
  
  list_predict_null[[paste(test_dataset, test_sample_name)]] = predict_null
  
}





for(j in 1:length(datasets))
{
  if(i == j) {next}
  
  model_dataset = datasets[j]
  model_sample_name = sample_names[j]
  
  cat('train dataset: ', model_dataset, ' train_sample_name: ', model_sample_name,  '\n')
  
  
  ####Test CNN###
  #patience = 20
  #dropout_rate = 0.2
  #hidden_act = 'relu'
  #output_act = 'tanh'
  #loss = 'mse'
  
  
  patience = 10
  dropout_rate = 0.2
  hidden_act = 'relu'
  output_act = 'linear'
  loss = 'mse'
  
  
  #patience = 10
  #dropout_rate = 0.5
  #hidden_act = 'relu'
  #output_act = 'sigmoid'
  #loss = 'mse'
  
  
  #patience = 10
  #dropout_rate = 0.3
  #hidden_act = 'relu'
  #output_act = 'sigmoid'
  #loss = 'mse'
  
  
  #patience = 10
  #dropout_rate = 0.3
  #hidden_act = 'relu'
  #output_act = 'tanh'
  #loss = 'mse'
  
  cnn_model_dir = paste0(package_data_dir, '/models/',model_dataset,'/cnn_model/')
  
  
  cnn_model_file = paste0(cnn_model_dir,
                          'cnn_model.patience_',  patience,
                          '.hidden_act_', hidden_act,
                          '.output_act_',output_act,
                          '.loss_', loss,
                          '.tr_size_', train_subsample_size,
                          '.norm_scope_', norm_scope,
                          '.sample_name_', model_sample_name, 
                          '.hd5')
  
  
  cnn_model = load_model_hdf5(cnn_model_file)
  
  #somehow parallelization does not work for CNN
  # if(num_elements > 1000000)
  # {
  #   
  #   predict_temp_list = parLapply(cl, sub_array_list, cnn_predict, cnn_model = cnn_model )
  #   
  #   #predict_elastic = do.call(cBind, predict_temp_list)
  #   head(predict_temp_list[[1]])
  #   predict_cnn = Reduce(c, predict_temp_list)
  #   head(predict_cnn)
  #   length(predict_cnn)
  #   
  # }else
  # {
  #   predict_cnn_res <- cnn_model %>% predict(fr_list_test_to_be_used$features_array)
  #   predict_cnn = predict_cnn_res[,1]
  # }
  
  predict_cnn_res <- cnn_model %>% predict(fr_list_test_to_be_used$features_array)
  predict_cnn = predict_cnn_res[,1]
  
  predict_cnn = minmax_norm(predict_cnn)
  
  names(predict_cnn) = rownames(fr_list_test_to_be_used$features_matrix)
  head(predict_cnn)
  length(predict_cnn)
  sum(is.na(predict_cnn))
  head(which(is.na(predict_cnn)))
  predict_cnn[176]
  fr_list_test_to_be_used$features_matrix[176, ]
  sum(is.na(predict_cnn))
  
  #predict_cnn = predict_cnn[!is.na(predict_cnn)]
  
  length(predict_cnn)
  max(predict_cnn)
  min(predict_cnn)
  
  list_predict_cnn[[paste(model_dataset, model_sample_name, test_dataset, test_sample_name)]] = predict_cnn
  
  #hist(predict_cnn)
  if(flag_true_exp_exists )
  {
    spear_test_cnn = cor(fr_list_test_to_be_used$response, predict_cnn, method = 'spearman')
    spear_test_cnn = round(spear_test_cnn, 3)
    
    error_test_cnn = mean(abs(fr_list_test_to_be_used$response - predict_cnn), na.rm = T)
    error_test_cnn = round(error_test_cnn, 3)
    
    print(paste('Global test correlation (CNN):', spear_test_cnn) )
    
    
    cor_mat_cnn[paste(model_dataset, model_sample_name), paste(test_dataset, test_sample_name)]  = spear_test_cnn
    error_mat_cnn[paste(model_dataset, model_sample_name), paste(test_dataset, test_sample_name)]  = error_test_cnn
    
  }
   
  #Split large matrix
  
  

  
  ####Test Elastic###
  alpha = 0.5
  
  elastic_model_dir = paste0(package_data_dir, '/models/',model_dataset,'/elastic_model/')
  
  elastic_model_file = paste0(elastic_model_dir, 'elastic_model.alpha_', alpha , '.tr_size_', train_subsample_size, 
                              '.norm_scope_', norm_scope, '.sample_name_', model_sample_name,  '.rds')
  elastic_model = readRDS(file = elastic_model_file)
  
  dim(fr_list_test_to_be_used$features_matrix)
  
  
  if(num_elements > 1000000)
  {
    
    predict_temp_list = parLapply(cl, sub_mat_list, elastic_predict, elastic_model = elastic_model )
    
    #predict_elastic = do.call(cBind, predict_temp_list)
    head(predict_temp_list[[1]])
    predict_elastic = Reduce(c, predict_temp_list)
    head(predict_elastic)
    length(predict_elastic)
    
  }else
  {
    
    predict_elastic = predict(elastic_model, s = "lambda.min", 
                              fr_list_test_to_be_used$features_matrix )[,1]
    
  }
  
  
  predict_elastic[predict_elastic < 0] = 0
  sum(is.na(predict_elastic))
  predict_elastic = predict_elastic[!is.na(predict_elastic)]
  
  predict_elastic = predict_elastic / max(predict_elastic, na.rm = T)
  names(predict_elastic) = rownames(fr_list_test_to_be_used$features_matrix_test)
  list_predict_elastic[[paste(model_dataset, model_sample_name, test_dataset, test_sample_name)]] = predict_elastic
  
  
  if(flag_true_exp_exists )
  {
    spear_elastic = cor(fr_list_test_to_be_used$response, predict_elastic, method = 'spearman')
    spear_elastic = round(spear_elastic, 3)
    
    error_test_elastic = mean(abs(fr_list_test_to_be_used$response - predict_elastic), na.rm = T )
    error_test_elastic = round(error_test_elastic, 3)    
    
    print(paste('Global test correlation (Elastic):', spear_elastic) )
    
    cor_mat_elastic[paste(model_dataset, model_sample_name), paste(test_dataset, test_sample_name)]  = spear_elastic
    error_mat_elastic[paste(model_dataset, model_sample_name), paste(test_dataset, test_sample_name)]  = error_test_elastic
  }
  
  
  ####Test RF###
  
  
  sampsize = 0.8
  nodesize = 1
  
  rf_model_dir = paste0(package_data_dir, '/models/',model_dataset,'/rf_model/')
  
  rf_model_file = paste0(rf_model_dir, 'rf_model.sampsize_', sampsize ,'.nodesize_', nodesize , '.tr_size_', train_subsample_size,
                         '.norm_scope_', norm_scope, '.sample_name_', model_sample_name,  '.rds')
  
  rf_model = readRDS(rf_model_file)
  
  if(num_elements > 1000000)
  {

    predict_temp_list = parLapply(cl, sub_mat_list, random_forest_predict, rf_model = rf_model )
    #predict_elastic = do.call(cBind, predict_temp_list)
    head(predict_temp_list[[1]])
    predict_rf = Reduce(c, predict_temp_list)
    head(predict_rf)
    length(predict_rf)
    
  }else
  {
    
    predict_rf = predict(rf_model, newdata = fr_list_test_to_be_used$features_matrix) 
    
  }
  names(predict_rf) = rownames(fr_list_test_to_be_used$features_matrix)                     
  head(predict_rf)
  predict_rf = predict_rf / max(predict_rf, na.rm = T)
  
  list_predict_rf[[paste(model_dataset, model_sample_name, test_dataset, test_sample_name)]] = predict_rf
  
  if(flag_true_exp_exists )
  {
    spear_rf = cor(fr_list_test_to_be_used$response, predict_rf, method = 'spearman')
    spear_rf = as.numeric(round(spear_rf, 3))
    
    error_test_rf = mean(abs(fr_list_test_to_be_used$response - predict_rf), na.rm = T )
    error_test_rf = round(error_test_rf, 3)  
    
    print(paste('Global test correlation (RF):', spear_rf) )
    
    cor_mat_rf[paste(model_dataset, model_sample_name), paste(test_dataset, test_sample_name)]  = spear_rf
    error_mat_rf[paste(model_dataset, model_sample_name), paste(test_dataset, test_sample_name)]  = error_test_rf
  }
  
  
  
  ###Ensemb###
  predict_ensem = (predict_cnn + predict_rf + predict_elastic) / 3
  names(predict_ensem) = rownames(fr_list_test_to_be_used$features_matrix)                     
  predict_ensem = predict_ensem / max(predict_ensem, na.rm = T)
  list_predict_ensem[[paste(model_dataset, model_sample_name, test_dataset, test_sample_name)]] = predict_ensem
  
  
  if(flag_true_exp_exists )
  {
    spear_ensem = cor(fr_list_test_to_be_used$response, predict_ensem, method = 'spearman')
    spear_ensem = as.numeric(round(spear_ensem, 3))
    
    error_test_ensem = mean(abs(fr_list_test_to_be_used$response - predict_ensem) , na.rm = T )
    error_test_ensem = round(error_test_ensem, 3)  
    
    print(paste('Global test correlation (Ensem):', spear_ensem) )
    
    cor_mat_ensem[paste(model_dataset, model_sample_name), paste(test_dataset, test_sample_name)]  = spear_ensem
    error_mat_ensem[paste(model_dataset, model_sample_name), paste(test_dataset, test_sample_name)]  = error_test_ensem
    
  }
  
  
  list_predicts = list(list_predict_cnn, list_predict_elastic, list_predict_rf, list_predict_ensem, list_predict_mean, list_predict_null)
  
  list_predicts_file = paste0(package_data_dir, '/models/', 'predictions', '.tr_size_', train_subsample_size,  '.ts_size_', test_subsample_size, 
                              '.norm_scope_', norm_scope,'.v10.',test_dataset,'.rds')
  
  saveRDS(list_predicts, file = list_predicts_file)
  
  
  
}



#}

diag(cor_mat_cnn) = diag(cor_mat_rf) = diag(cor_mat_elastic) = diag(cor_mat_ensem) = NA

cor_mat_cnn
cor_mat_elastic
cor_mat_rf
cor_mat_ensem

cor_vec_mean
cor_vec_null

error_mat_cnn
error_mat_elastic
error_mat_rf
error_mat_ensem

error_vec_mean



list_cor_mats = list(cor_mat_cnn, cor_mat_elastic, cor_mat_rf, cor_mat_ensem, cor_vec_mean, cor_vec_null)
list_error_mats = list(error_mat_cnn, error_mat_elastic, error_mat_rf, error_mat_ensem, error_vec_mean, error_vec_null)
list_predicts = list(list_predict_cnn, list_predict_elastic, list_predict_rf, list_predict_ensem, list_predict_mean, list_predict_null)


names(list_cor_mats) = names(list_error_mats) = names(list_predicts) = c('CNN', 'Elastic', 'RF', 'Ensemble', 'Mean', 'Null')

list_cor_mats_file = paste0(package_data_dir, '/models/', 'cor_mats', '.tr_size_', train_subsample_size,  '.ts_size_', test_subsample_size,  '.norm_scope_', norm_scope,'.v02.rds')
saveRDS(list_cor_mats, file = list_cor_mats_file)

list_error_mats_file = paste0(package_data_dir, '/models/', 'error_mats', '.tr_size_', train_subsample_size,  '.ts_size_', test_subsample_size,  '.norm_scope_', norm_scope,'.v02.rds')
saveRDS(list_error_mats, file = list_error_mats_file)


#list_predicts_file = paste0(package_data_dir, '/models/', 'predictions', '.tr_size_', train_subsample_size,  '.ts_size_', test_subsample_size, 
#                            '.norm_scope_', norm_scope,'.v02.rds')

list_predicts_file = paste0(package_data_dir, '/models/', 'predictions', '.tr_size_', train_subsample_size,  '.ts_size_', test_subsample_size, 
                            '.norm_scope_', norm_scope,'.v10.',test_dataset,'.rds')

saveRDS(list_predicts, file = list_predicts_file)


list_response_test_file = paste0(package_data_dir, '/models/', 'response_test', '.ts_size_', test_subsample_size,  '.norm_scope_', norm_scope,'.v02.rds')
saveRDS(list_response_test, file = list_response_test_file)




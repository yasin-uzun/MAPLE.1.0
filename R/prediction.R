library(data.table)
library(keras)
library(glmnet)



cnn_predict <- function(fr_list, model_file)
{

  library(keras)
  cnn_model = load_model_hdf5(model_file)
  #result = predict(cnn_model, newdata)
  temp_result <- cnn_model %>% predict(fr_list$features_array)
  result <- temp_result[,1]
  length(result)
  head(result)
  names(result) = rownames(fr_list$features_matrix)

  return(result)
}

elastic_predict <- function(fr_list, model_file)
{
  library(glmnet)
  elastic_model = readRDS(file = model_file)
  temp_result = predict(elastic_model, s = "lambda.min",
                   fr_list$features_matrix )
  result <- temp_result[,1]
  length(result)
  head(result)
  return(result)
}


rf_predict <- function(fr_list, model_file)
{
  library(randomForest)
  rf_model = readRDS(file = model_file)
  result = predict(rf_model, newdata = fr_list$features_matrix )
  length(result)
  head(result)
  return(result)

}


ensemble_predict <- function(prediction_list)
{
   list_rows    = lapply(prediction_list, names)
   shared_rows = Reduce(intersect, list_rows)
   length(shared_rows)
   head(shared_rows)

   prediction_mat = sapply(prediction_list, function (x) x[shared_rows])
   dim(prediction_mat)
   head(prediction_mat)

   result = rowMeans(prediction_mat, na.rm = T)

   head(result)
   length(result)

   return(result)
}


convert_preds_to_matrix <- function(predictions)
{

  df_predict <- sapply(names(predictions), function(x) strsplit(x, "_x_")[[1]], USE.NAMES=FALSE)
  dim(df_predict)
  df_predict[1:2, 1:3]
  df_predict = data.frame(t(df_predict), predictions)
  head(df_predict)
  colnames(df_predict) = c('cell', 'gene_symbol', 'activity')
  dna_mat = reshape2::dcast(df_predict, gene_symbol~cell)
  dna_mat[1:5, 1:5]
  dim(dna_mat)


  dna_mat[1:5, 1:5]
  rownames(dna_mat) = dna_mat$gene_symbol
  dna_mat$gene_symbol = NULL

  sum(is.na(dna_mat))

  dna_mat <- apply(dna_mat, MARGIN = 2, FUN = function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })

  sum(is.na(dna_mat))


  dna_mat[1:5, 1:5]
  dim(dna_mat)
  sum(is.na(dna_mat[1:5, 1:5]))

  head(sort(rowSums(dna_mat)))
  sds = apply(dna_mat, MARGIN = 1, FUN = sd)
  head(sds)
  head(sort(sds))
  dna_mat = dna_mat[sds > 0,]
  dim(dna_mat)

  #rownames(dna_mat) = gsub(' ', '_', rownames(dna_mat))
  return(dna_mat)

}


#
# subsample_fr <- function(fr_list, subsample_size = 10000, seed = 123)
# {
#
#   response = fr_list$response
#   features_array = fr_list$features_array
#   features_matrix = fr_list$features_matrix
#   features_mean = fr_list$features_mean
#   features_mean_5 = fr_list$features_mean_5
#
#   set.seed(seed)
#   subsample_index = sample(subsample_size)
#
#   features_matrix_subsampled = features_matrix[subsample_index, ]
#   response_subsampled = response[subsample_index ]
#   features_array_subsampled = features_array[subsample_index, ,]
#   features_mean_subsampled = features_mean[subsample_index]
#   features_mean_5_subsampled = features_mean_5[subsample_index]
#
#   fr_list$features_matrix = features_matrix_subsampled
#   fr_list$response = response_subsampled
#   fr_list$features_array = features_array_subsampled
#   fr_list$features_mean = features_mean_subsampled
#   fr_list$features_mean_5 = features_mean_5_subsampled
#
#   return(fr_list)
#
# }
#
#
#


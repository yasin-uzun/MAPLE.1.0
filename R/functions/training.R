library(data.table)
library(keras)
library(glmnet)



train_cnn_model <- function(array_train, response_train, 
                            patience = 10, 
                            hidden_act = 'relu', 
                            output_act = 'sigmoid', 
                            loss = 'mse',
                            dropout_rate = 0,  
                            epochs = 500)
{
  
  result_list = list()
  
  set.seed(123)
  
  print_dot_callback <- callback_lambda(
    on_epoch_end = function(epoch, logs) {
      if (epoch %% 100 == 0) cat("\n")
      cat(".")
    }
  ) 
  
  early_stop <- callback_early_stopping(monitor = "val_loss", patience = patience)
  
  set.seed(123)
  
 
  
  use_session_with_seed(123)
  
  if(dropout_rate == 0)
  {
    model <- keras_model_sequential() %>%
      layer_conv_1d(filters = 50, kernel_size = 5, activation = hidden_act, input_shape = dim(array_train[1,,]) ) %>%
      layer_max_pooling_1d(pool_size = 4)  %>%
      layer_flatten()  %>%
      layer_dense(units = 100, activation = hidden_act) %>% 
      layer_dense(units = 50, activation = hidden_act) %>% 
      layer_dense(units = 25, activation = hidden_act) %>% 
      layer_dense(units = 12, activation = hidden_act) %>% 
      layer_dense(units = 6, activation = hidden_act) %>% 
      layer_dense(units = 3, activation = hidden_act) %>% 
      layer_dense(units = 1, activation = output_act)     
  }else
  {

    model <- keras_model_sequential() %>%
      layer_conv_1d(filters = 50, kernel_size = 5, activation = hidden_act, input_shape = dim(array_train[1,,]) ) %>%
      layer_max_pooling_1d(pool_size = 4)  %>%
      layer_flatten()  %>%
      layer_dropout(rate = dropout_rate) %>% 
      layer_dense(units = 100, activation = hidden_act) %>% 
      layer_dropout(rate = dropout_rate) %>% 
      layer_dense(units = 50, activation = hidden_act) %>% 
      layer_dropout(rate = dropout_rate) %>% 
      layer_dense(units = 25, activation = hidden_act) %>% 
      layer_dropout(rate = dropout_rate) %>% 
      layer_dense(units = 12, activation = hidden_act) %>% 
      layer_dropout(rate = dropout_rate) %>% 
      layer_dense(units = 6, activation = hidden_act) %>% 
      layer_dropout(rate = dropout_rate) %>% 
      layer_dense(units = 3, activation = hidden_act) %>% 
      layer_dense(units = 1, activation = output_act) 
  }
  
  
  
  model %>% compile(
    loss = loss,
    optimizer = optimizer_rmsprop(),
    metrics = list("mean_absolute_error")
  )
  
  set.seed(123)
  
  history <- model %>% fit(
    array_train,
    response_train,
    epochs = epochs,
    validation_split = 0.2,
    verbose = 0,
    callbacks = list(early_stop, print_dot_callback)
  )
  
  return(list('model' = model, 'history' = history))
  
}


cnn_predict <- function(newdata, cnn_model)
{
  
  library(keras)
  result = predict(cnn_model, newdata)
  return(result[,1] )
}

elastic_predict <- function(newdata, elastic_model)
{
  library(glmnet)
  return(predict(elastic_model, s = "lambda.min", 
                 newdata )[,1])
}


random_forest_predict <- function(newdata, rf_model)
{
  library(randomForest)
  return(predict(rf_model, newdata = newdata))
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



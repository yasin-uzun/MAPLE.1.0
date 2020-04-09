library(reshape2)
library(data.table)
library(abind)

#input_binned_met = input_binned_met_test
#df_cpg_content = df_cpg_content
#input_mean_met = input_mean_met_test
#input_exp = input_exp_test
#gene_set = gene_set_test

#float_values = response_train

split_values_into_classes <- function(float_values, numclasses = 5)
{
  sorted = sort(float_values)
  cut_idx = round(1:num_classes / num_classes * length(float_values))
  cut_points = sorted[cut_idx]

  classes = rep(0, length(float_values))

  for(i in 1:num_classes)
  {
    classes[float_values>cut_points[i]] = i
    names(classes) = names(float_values)
  }
  return(factor(classes))
}



compute_feature_matrix <- function(input_binned_met, df_cpg_content, gene_set = NULL, max_na_bins = 5)
{
  gene_set = intersect(gene_set, rownames(df_cpg_content))
  length(gene_set)

  input_binned_met_gs = input_binned_met[input_binned_met$gene_symbol %in% gene_set, ]
  head(input_binned_met_gs)
  dim(input_binned_met_gs)
  sum(is.na(input_binned_met_gs))

  length(unique(input_binned_met_gs$gene_symbol))


  count_na_met_by_rows = rowSums(is.na(input_binned_met_gs))
  #nonna_met_features = input_binned_met_gs[count_na_met_by_rows == 0, ]
  nonna_met_features = input_binned_met_gs[count_na_met_by_rows <= max_na_bins, ]

  nonna_met_features[is.na(nonna_met_features)] = mean(as.matrix(nonna_met_features[, 3:ncol(nonna_met_features)]), na.rm = T)

  dim(nonna_met_features)
  dim(df_cpg_content)

  features = as.matrix(cbind(nonna_met_features[3:ncol(nonna_met_features)], df_cpg_content[as.character(nonna_met_features$gene_symbol),] ))
  head(features)

  indx <- which(is.na(features), arr.ind=TRUE)
  head(indx)
  features[1153,]
  nonna_met_features[1153,]

  dim(features)
  head(features)

  rownames(features) = paste0(nonna_met_features$cell, '_x_', nonna_met_features$gene_symbol)

  return(features)
}


feature_matrix_to_3d_array <- function(features)
{
  print('Transforming matrix into array')
  num_tot_features = ncol(features)
  num_bins = num_tot_features/ 2
  head(features)
  dim(features)
  #array_3d = array(dim = c(nrow(features), num_bins, 2))
  #for(i in 1:nrow(features))
  #{
  #  array_3d[i,,] = cbind(features[i, 1:num_bins], features[i, (num_bins+1):num_tot_features])
  #}
  array_3d = abind( features[,1:num_bins], features[,(num_bins+1):num_tot_features], along=3 )

  return(array_3d)
}



compute_mean_met_input <- function(input_mean_met, gene_set)
{
  input_mean_met = input_mean_met[gene_set, ]
  input_mean_met_melted = melt(as.matrix(input_mean_met))
  head(input_mean_met_melted)
  colnames(input_mean_met_melted) = c('gene_symbol', 'cell', 'mean_met')

  rownames(input_mean_met_melted) = paste0(input_mean_met_melted$cell, '_x_', input_mean_met_melted$gene_symbol)
  head(input_mean_met_melted)
  mean_met_vec_input = as.numeric(input_mean_met_melted$mean_met)
  names(mean_met_vec_input) = rownames(input_mean_met_melted)
  head(mean_met_vec_input)
  return(mean_met_vec_input)

}


compute_response_vector <- function(input_exp)
{
  input_exp_melted = melt(as.matrix(input_exp))
  colnames(input_exp_melted) = c('gene_symbol', 'cell', 'expr')
  rownames(input_exp_melted) = paste0(input_exp_melted$cell, '_x_', input_exp_melted$gene_symbol)
  head(input_exp_melted)
  exp_vec_input = as.numeric(input_exp_melted$expr)
  names(exp_vec_input) = rownames(input_exp_melted)
  head(exp_vec_input)
  return(exp_vec_input)

}




compute_features_and_responses <- function(input_binned_met, df_cpg_content, input_mean_met = NULL, input_mean_met_5 = NULL,
                                           input_exp = NULL, gene_set = NULL, max_na_bins = 5 )
{
  print('Computing features matrix')
  features_matrix = compute_feature_matrix(input_binned_met, df_cpg_content, gene_set, max_na_bins = max_na_bins)
  print('Transforming features matrix into array')
  features_array = feature_matrix_to_3d_array(features_matrix)

  features_mean = compute_mean_met_input(input_mean_met, gene_set)
  features_mean = features_mean[rownames(features_matrix)]

  if(is.null(input_mean_met_5))
  {
    features_mean_5 = NULL
  }else{
    features_mean_5 = compute_mean_met_input(input_mean_met_5, gene_set)
    features_mean_5 = features_mean_5[rownames(features_matrix)]
  }

  sum(is.na(features_matrix))
  head(features_mean)

  #names(features_mean) = rownames(features_matrix)

  dimnames(features_array)[[1]] = rownames(features_matrix)


  #df_cell_names <- sapply(temp123, function(x) strsplit(x, "_x_")[[1]], USE.NAMES=FALSE)
  #tail(sort(unique(df_cell_names[2,])))
  if(is.null(input_exp) || nrow(input_exp) == 0 )
  {
    response = NULL
  }else
  {
    response = compute_response_vector(input_exp)
    response = response[rownames(features_matrix)]
  }
  result_list = list('features_mean' = features_mean,
                     'features_mean_5' = features_mean_5,
                     'features_matrix' = features_matrix,
                     'features_array' = features_array,
                     'response' = response
  )

  return(result_list)
}

minmax_norm <- function(x)
{
  return(   (x - min(x, na.rm = T) )  /  ( max(x, na.rm = T) - min(x, na.rm = T)  ) )
}

normalize_rna_counts <- function(count_matrix, norm_scope = 'matrix')
{

  df_cpm = sweep(count_matrix, 2, colSums(count_matrix, na.rm = T),`/`) * 1000000
  #df_cpm[1:5, 1:5]
  dim(df_cpm)

  log_cpm = as.matrix(log10(df_cpm + 1))

  if(norm_scope == 'matrix')
  {
    min1 = min(as.numeric(log_cpm), na.rm = T)
    max1 = max(as.numeric(log_cpm), na.rm = T)

    log_cpm_minmax = (log_cpm - min1) / (max1- min1)
  }else if(norm_scope == 'gene')
  {
    mins = apply(log_cpm, MARGIN = 1, FUN = 'min')
    maxs = apply(log_cpm, MARGIN = 1, FUN = 'max')

    log_cpm_minmax = t(apply(log_cpm, 1, function(x)(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))))

    log_cpm_minmax[is.na(log_cpm_minmax)] = 0
  }else
  {
    stop('RNA expression normalization can be done either at matrix level or gene level.')
  }
  return(log_cpm_minmax)
}



standardize <- function(x)
{
  return(   (x - mean(x, na.rm = T) )  / sd(x, na.rm = T)  )
}


standardize_met_data <- function(met_matrix, norm_scope = 'matrix')
{
  demet_matrix = as.matrix(1 - as.matrix(met_matrix[, 3:ncol(met_matrix)]))

  demet_matrix <- apply(demet_matrix, MARGIN = 2, FUN = function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })

  if(norm_scope == 'matrix')
  {
    mean_glob = mean(as.numeric(demet_matrix), na.rm = T)
    sd_glob = mean(as.numeric(demet_matrix), na.rm = T)

    demet_std = (demet_matrix - mean_glob) / sd_glob
    std_demet_df = data.frame(met_matrix[, 1:2], demet_std)

  }else if(norm_scope == 'gene')
  {
    demet_df = data.frame(met_matrix[, 1:2], demet_matrix)

    demet_std <- ave(as.matrix(demet_df[, 3:ncol(demet_df)]), demet_df$gene_symbol, FUN=standardize    )
    std_demet_df = data.frame(met_matrix[, 1:2], demet_std)
    #std_demet_df = data.frame(met_matrix[, 1:2])

  }else
  {
    stop('Methylation normalization can be done either at matrix level or gene level.')
  }

  return(std_demet_df)

}






get_fr_list <- function(meta_data, cpg_content_file, max_na_bins = 0, norm_scope = 'matrix')
{

  df_cpg_content = readRDS(cpg_content_file)
  head(df_cpg_content)
  dim(df_cpg_content)

  input_binned_met = meta_data@df_rate_meta
  dim(input_binned_met)
  head(input_binned_met)

  input_exp = meta_data@df_exp_meta
  sum(rowSums(input_exp) == 0)
  dim(input_exp)

  sum(rowSums(input_exp) > 0)


  #par(mfrow = c(2,1))
  input_exp_norm = normalize_rna_counts(input_exp, norm_scope = norm_scope)
  #hist(input_exp_matrix_normalized[input_exp_matrix_normalized > 0],  main = 'Global Normalization', xlab = 'Normalized expression')

  #input_exp_gene_normalized = normalize_rna_counts(input_exp, norm_scope = 'gene')
  #hist(input_exp_gene_normalized[input_exp_gene_normalized > 0], main = 'Gene Normalization', xlab = 'Normalized expression')


  input_mean_met = meta_data@df_rate_meta_mean
  input_mean_met_5 = meta_data@df_rate_meta_mean_5

  dim(input_mean_met)
  head(input_mean_met)

  print('Processing met. data')
  std_demet = standardize_met_data(met_matrix = input_binned_met, norm_scope = norm_scope)

  #gene_set = rownames(input_exp_norm)
  gene_set = unique(std_demet$gene_symbol)
  head(gene_set)
  gene_set = intersect(gene_set, rownames(df_cpg_content))

  print('Computing features')

  fr_list = compute_features_and_responses(input_binned_met = std_demet,
                                           df_cpg_content = df_cpg_content,
                                           input_mean_met = input_mean_met,
                                           input_mean_met_5 = input_mean_met_5,
                                           input_exp = input_exp_norm,
                                           gene_set = gene_set,
                                           max_na_bins = 0)
  print('Features computed.')

  return(fr_list)

}




subsample_fr <- function(fr_list, subsample_size = 10000, seed = 123)
{


  set.seed(seed)
  subsample_index = sample(length(fr_list$response), subsample_size)
  length(subsample_index)

  fr_list$features_mean = fr_list$features_mean[subsample_index ]
  fr_list$features_mean_5 = fr_list$features_mean_5[subsample_index ]

  fr_list$response = fr_list$response[subsample_index]
  #fr_list$response = fr_list$response[subsample_index]

  fr_list$features_matrix = fr_list$features_matrix[subsample_index, ]
  #fr_list$features_matrix = fr_list$features_matrix[subsample_index, ]

  fr_list$features_array = fr_list$features_array[subsample_index,, ]
  #fr_list$features_array = fr_list$features_matrix[subsample_index,, ]


  return(fr_list)

}








merge_fr_lists_for_multiple_datasets <- function(fr_lists_multiple_datasets)
{
  features_mean_list = lapply( fr_lists_multiple_datasets , "[[" , "features_mean" )
  features_mean5_list = lapply( fr_lists_multiple_datasets , "[[" , "features_mean5" )
  features_matrix_list = lapply( fr_lists_multiple_datasets , "[[" , "features_matrix" )
  features_array_list = lapply( fr_lists_multiple_datasets , "[[" , "features_array" )
  response_list = lapply( fr_lists_multiple_datasets , "[[" , "response" )

  merged_features_mean <- unlist(features_mean_list)
  dim(merged_features_mean)
  length(merged_features_mean)

  merged_features_mean5 <- unlist("cbind", features_mean5_list)
  dim(merged_features_mean5)
  length(merged_features_mean5)

  merged_features_matrix <- do.call("rbind", features_matrix_list)
  options(stringsAsFactors=FALSE)
  merged_features_array <- do.call(abind, args = list(features_array_list, along = 1))
  merged_response <- unlist(response_list)
  dim(merged_response)
  length(merged_response)

  merged_fr_list = list(features_mean = merged_features_mean
                        , features_mean5 = merged_features_mean5
                        , features_matrix = merged_features_matrix
                        , features_array = merged_features_array
                        , response = merged_response
  )
  return(merged_fr_list)
}











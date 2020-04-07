
class_meta_cells <- setClass("class_meta_cells",
                             slots = c(
                               df_met_counts = "data.frame",
                               df_demet_counts = "data.frame",
                               df_rate_meta = "data.frame",
                               df_rate_meta_mean = "data.frame",
                               df_rate_meta_mean_5 = "data.frame",
                               
                               df_exp_meta = "data.frame",
                               num_neighbors = "numeric",
                               pca_coords = "matrix"
                             )  
                          )


class_model_pred <- setClass("class_model_pred",
                             slots = c(
                               list_models = "list",
                               list_abs_error = "list",
                               list_error_rates = "list",
                               list_predictions = "list",
                               response_test = "numeric",
                               global_correlations = "numeric",
                               median_cors = "numeric",
                               
                               list_gene_correlations = "list"
                             )  
                          )





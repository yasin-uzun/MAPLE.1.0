# MethylPredict
An R package for predicting gene activity level for single cell DNA Methylation data

MethylPredict is a supervised learning algorithm that is developed the predict the gene activity levels of individual cells from single cell DNA Methylation. It is implemented in R programming language and uses CpG methylation sites in promoter regions to infer gene activity. 

## Installation

To install MethylPredict, type the following commands in R command windows:
```
library(devtools)
install_github("yasin-uzun/MethylPredict.1.0")
```
## Annotation Data

To run MethylPredict, you need some genome annotation files. You can download those annotation files using these links:

## Example Data

We provide an example dataset for you to test MethylPredict. This data is a subset of CpG files from the study of Luo et al (1).   
The original data has over 3000 cells. In here, we just provide data 100 cells to give a quick run.

## Running
Assuming you have downloaded the test data and you have the CpG cov files (in Bismark format), you can run MethylPredict as follows:

```R
library(MethylPredict)

#Set directory names
annot_dir = 'data/annot/'
model_dir = 'data/models/'
cov_dir = 'data/example/cov_files/'

#Set input files
annot_file = paste0(annot_dir,'/gencode.mm10.vM22.genes.bed')
cpg_content_file =  paste0(annot_dir,'/regions.genes.tss_ud_5K.cpg_ratio.bin_size_500.mm10.rds')

#Compute binned data
binned_list = compute_binned_met_counts(cov_dir = cov_dir, annot_file = annot_file )

#Compute meta cells
meta_object = compute_meta_cells(df_met =  binned_list[["df_binned_met"]],
                                 df_demet =  binned_list[["df_binned_demet"]])
                                 
#Generate features                               
fr_list = get_fr_list(meta_data = meta_object, cpg_content_file = cpg_content_file)

#Load CNN model and predict
cnn_model_file = paste0(model_dir,'/Clark/cnn_model.patience_10.hidden_act_relu.output_act_linear.loss_mse.hd5')
predict_cnn = cnn_predict(fr_list, cnn_model_file)

#Load Elastic model and predict
elastic_model_file = paste0(model_dir,'/Clark/elastic_model.alpha_0.5.rds')
predict_elastic = elastic_predict(fr_list, elastic_model_file)

#Load RF model and predict
rf_model_file = paste0(model_dir,'/Clark/rf_model.sampsize_0.8.rds')
predict_rf = rf_predict(fr_list, rf_model_file)

#Compute Ensemble prediction
prediction_list = list(predict_cnn, predict_elastic, predict_rf)
predict_ensem = ensemble_predict(prediction_list)

#Convert prediction into matrix format (genesxcells)
dna_mat = convert_preds_to_matrix(predict_ensem)


```

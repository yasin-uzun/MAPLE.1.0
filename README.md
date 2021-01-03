# MAPLE: Methylome Association by Predictive Linkage to Expression
An R package for predicting gene activity level for single cell DNA Methylation data

MAPLE is a supervised learning algorithm that is developed tO predict the gene activity levels of individual cells from single cell DNA Methylation. It is implemented in R and uses CpG methylation sites in promoter regions to infer the activity levels of the genes. 

## Installation

To install MAPLE, type the following commands in R command prompt:
```R
library(devtools)
install_github("yasin-uzun/MAPLE.1.0")
```
## Annotation Data

To run MethylPredict, you need some genome annotation files. You can download those annotation from here: 

[Annotation Data](https://drive.google.com/drive/folders/1E2UH8K02n-RSY67JCDN__ghArtwZ3e0x?usp=sharing)

## Built Models
We trained multiple models using different multi-omics training data. You can download them from here:

[Trained Models](https://drive.google.com/drive/folders/1lLHOsT-gaPwwdYWjexDhfXSH9QjKQDAa?usp=sharing)

## Example Data

We provide an example dataset o test MethylPredict. This data is a subset of CpG files from the study of Luo et al \[1\] (GSE97179).   
The original data has over 3000 cells. In here, we just provide data 100 cells to give a quick run:

[Example Dataset](https://drive.google.com/drive/folders/1GMJpbrXHJ91iPmyclGIPVhxoLcfqyjum?usp=sharing)

## Running
Assuming you have downloaded the test data and you have the CpG cov files (in Bismark format), you can run MethylPredict as follows:

```R
library(MAPLE)

#Set directory names
annot_dir = 'data/annot/'
model_dir = 'data/models/Clark/'
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
cnn_model_file = paste0(model_dir, '/cnn_model.hd5')
predict_cnn = cnn_predict(fr_list, cnn_model_file)

#Load Elastic model and predict
elastic_model_file = paste0(model_dir, '/elastic_model.rds')
predict_elastic = elastic_predict(fr_list, elastic_model_file)

#Load RF model and predict
rf_model_file = paste0(model_dir, '/rf_model.rds')
predict_rf = rf_predict(fr_list, rf_model_file)

#Compute Ensemble prediction
prediction_list = list(predict_cnn, predict_elastic, predict_rf)
predict_ensem = ensemble_predict(prediction_list)

#Convert prediction into matrix format (genesxcells)
gene_activity_matrix = convert_preds_to_matrix(predict_ensem)


```
## Citation
If you use our software, please cite [our paper](https://genome.cshlp.org/content/early/2020/11/20/gr.267047.120.abstract) \[2\] .


## Contact
For comments and questions, please contact Yasin Uzun (uzuny at email chop edu)

## References
\[1\] Luo, Chongyuan, Christopher L. Keown, Laurie Kurihara, Jingtian Zhou, Yupeng He, Junhao Li, Rosa Castanon, et al. 2017. “Single-Cell Methylomes Identify Neuronal Subtypes and Regulatory Elements in Mammalian Cortex.” Science 357 (6351): 600–604

\[2\] Yasin Uzun, Hao Wu, Kai Tan. Predictive modeling of single-cell DNA methylome data enhances integration with transcriptome data. Genome Research. Published in Advance. November 20, 2020.


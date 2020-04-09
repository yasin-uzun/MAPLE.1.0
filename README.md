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

```
library(MethylPredict)

```

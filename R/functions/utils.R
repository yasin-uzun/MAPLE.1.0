

convert_dataset_names <- function(ds_vec, sample_name = NULL)
{
  if(is.null(sample_name))
  {
    ds_vec = gsub('scnmt_gastr', 'scNMT-Gastr.', ds_vec)
    
    ds_vec = gsub('scmandt_muscle', 'scM&T - Muscle', ds_vec)
    ds_vec = gsub('scmandt', 'scM&T - ESC', ds_vec)
    ds_vec = gsub('scnmt', 'scNMT - EB', ds_vec)
    print(123)
    
  }else
  {
    ds_vec = gsub('scmandt_muscle', 'scM&T', ds_vec)
    ds_vec = gsub('scmandt', 'scM&T', ds_vec)
    ds_vec = gsub('scnmt', 'scNMT', ds_vec)   
    
    sample_name_nice = sample_name
    if(sample_name == 'Serum') {sample_name_nice = 'ESC'}
    if(sample_name == 'GSM') {sample_name_nice = 'Muscle'}
    
    ds_vec = paste(ds_vec, sample_name_nice)
  }
  return(ds_vec)
}


convert_dataset_names_by_author <- function(ds_vec, sample_name = NULL)
{
  if(is.null(sample_name))
  {
    ds_vec = gsub('scnmt_gastr', 'Argelaguet et al.', ds_vec)
    
    ds_vec = gsub('scmandt_muscle', 'Hernando-Herraez et al.', ds_vec)
    ds_vec = gsub('scmandt', 'Angermueller et al.', ds_vec)
    ds_vec = gsub('scnmt', 'Clark et al.', ds_vec)
    print(123)
    
  }else
  {
    ds_vec = gsub('scmandt_muscle', 'scM&T', ds_vec)
    ds_vec = gsub('scmandt', 'scM&T', ds_vec)
    ds_vec = gsub('scnmt', 'scNMT', ds_vec)   
    
    sample_name_nice = sample_name
    if(sample_name == 'Serum') {sample_name_nice = 'ESC'}
    if(sample_name == 'GSM') {sample_name_nice = 'Muscle'}
    
    ds_vec = paste(ds_vec, sample_name_nice)
  }
  return(ds_vec)
}



get_protocol <- function(dataset)
{
  protocol= dataset
  if(dataset == 'scmandt') { protocol = 'scM&T'}
  if(dataset == 'scmandt_muscle') { protocol = 'scM&T'}
  if(dataset == 'scnmt') { protocol = 'scNMT'}
  if(dataset == 'sctrio') { protocol = 'scTrio'}
  return(protocol)
}

get_tissue <- function(sample_name)
{
  tissue = sample_name
  if(sample_name == 'Serum') { tissue = 'ESC'}
  if(sample_name == 'GSM') { tissue = 'Muscle'}
  if(sample_name == 'EB') { tissue = 'EB'}
  if(sample_name == 'CRC01') { tissue = 'CRC (HS)'}
  return(tissue)
}


get_model_name <- function(model_name)
{
  model = model_name
  if(model_name == 'cnn_model') { model = 'CNN'}
  if(model_name == 'elastic_model') { model = 'Elastic'}
  if(model_name == 'rf_model') { model = 'RF'}
  return(model)
}





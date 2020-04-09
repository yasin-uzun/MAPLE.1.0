library(data.table)
library(GenomicRanges)
library(rdist)
library(doSNOW)
library(reshape2)

#source('/mnt/isilon/tan_lab/uzuny/scripts/bed_ops/set_ops.R')



solidify <- function(df_input)
{
  rowsums = rowSums(df_input[, 3:ncol(df_input)])
  head(rowsums)
  sum(rowsums>5)

  df_totals = data.frame(cell = df_input$cell, gene = df_input$gene_symbol, count = rowsums)
  head(df_totals)
  df_matrix = dcast(df_totals, gene ~ cell)
  df_matrix[1:5, 1:5]
  rownames(df_matrix) = df_matrix$gene
  df_matrix$gene = NULL
  df_matrix[1:5, 1:5]
  df_matrix[is.na(df_matrix)] = 0
  return(as.matrix(df_matrix))
}

calc_met_rate <- function(df_met, df_demet)
{

  mat_met = solidify(df_met)
  mat_demet = solidify(df_demet)

  mat_met[1:5, 1:5]
  mat_demet[1:5, 1:5]

  shared_genes = intersect(rownames(mat_met), rownames(mat_demet))
  length(shared_genes)

  mat_met = mat_met[shared_genes, ]
  mat_demet = mat_demet[shared_genes, ]

  dim(mat_met)
  dim(mat_demet)

  mat_rate = mat_met / (mat_met + mat_demet)
  overall_means = rowSums(mat_met) / (rowSums(mat_met) + rowSums(mat_demet))

  sum(is.na(mat_rate))

  for(i in 1:nrow(mat_rate))
  {
    mat_rate[i, which(is.na(mat_rate[i,]) )]  = overall_means[i]
  }

  dim(mat_rate)
  length(overall_means)

  return(mat_rate)
}


compute_pca <- function(df_met, df_demet)
{

  met_rate <- calc_met_rate(df_met, df_demet)
  met_rate[1:5, 1:5]
  dim(met_rate)


  variance <- apply(met_rate, 1, var, na.rm=TRUE)
  mean <- apply(met_rate, 1, mean, na.rm=TRUE)
  bcv = variance/ mean
  head(bcv)
  bcv_sorted = sort(bcv, decreasing = T)

  veg_count = 5000

  top_bcv = head(bcv_sorted, veg_count)

  vegs = names(top_bcv)

  met_rate_veg = met_rate[vegs, ]

  pca <- prcomp(t(met_rate_veg), scale.=TRUE)#run the PCA analysis

  return(pca)
}


cos.sim <- function(A, B)
{
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}

compute_pca_distance_matrix <- function(pca, distance_dim_count, metric = 'euclidean')
{
  print('***')
  distance_dim_count = min(distance_dim_count, ncol(pca))
  pca_top_dims = pca[, 1:distance_dim_count]
  head(pca_top_dims)

  pca_distance_matrix = pdist(pca_top_dims, metric =  metric,  p = 2)
  rownames(pca_distance_matrix) = rownames(pca_top_dims)
  colnames(pca_distance_matrix) = rownames(pca_top_dims)
  pca_distance_matrix[1:5, 1:5]

  return(pca_distance_matrix)
}

#pca_distance_matrix = pca_distance_matrix_part_keep
#df_binned = df_met_part_keep
compute_meta_cells_methyl <- function(df_binned, pca_distance_matrix, num_neighbors, run_parallel = F)
{

  cells = colnames(pca_distance_matrix)
  head(cells)




  #cell = 'scTrioSeq2Met_CRC01_LN1_162'
  counter = 0
  dt_binned = data.table(df_binned)
  n = ncol(dt_binned)
  quant_cols = colnames(dt_binned)[3:n]
  df_meta = NULL

  if(run_parallel)
  {
    df_meta <- foreach(i=1:length(cells), .combine=rbind) %dopar%
      #df_meta <- foreach(i=1:4, .combine=rbind) %dopar%
    {
      start = Sys.time()

      library(data.table)
      library(rdist)

      #counter = counter + 1
      cell = as.character(cells[i])

      print(paste(i, cell))

      if(  sum(rownames(pca_distance_matrix) == cell) == 0 ) {next}

      closest_neighbors = names(head(sort(pca_distance_matrix[cell, ]), (num_neighbors + 1) ))

      dt_binned_sub = data.table(dt_binned[dt_binned$cell %in% closest_neighbors, ])
      dim(dt_binned_sub)
      print(class(dt_binned_sub))

      #dt_sub_aggr = df_binned_sub[, lapply(.SD, function(x) mean(x,na.rm=T) ) , by = gene_symbol, .SDcols = quant_cols  ]
      dt_sub_aggr = dt_binned_sub[, lapply(.SD, function(x) sum(x,na.rm=T) ) , by = gene_symbol, .SDcols = quant_cols  ]
      Sys.time()


      df_sub_aggr = data.frame(cell = cell, dt_sub_aggr)
      Sys.time()

      end = Sys.time()
      print(end - start)

      head(df_sub_aggr)
      #df_sub_aggr[1:5, 1:5]
      #imp_list[[cell]] = data.frame(cell = cell, df_sub_aggr)

      df_sub_aggr
    }


  }else
  {
    temp_list = list()
    start = Sys.time()
    for(i in 1:length(cells))
      #for(i in 1:30)
    {


      #counter = counter + 1
      cell = as.character(cells[i])

      print(paste(i, cell))

      if(  sum(rownames(pca_distance_matrix) == cell) == 0 ) {next}

      closest_neighbors = names(head(sort(pca_distance_matrix[cell, ]), (num_neighbors + 1) ))

      dt_binned_sub = data.table(dt_binned[dt_binned$cell %in% closest_neighbors, ])
      dim(dt_binned_sub)
      print(class(dt_binned_sub))

      #dt_sub_aggr = df_binned_sub[, lapply(.SD, function(x) mean(x,na.rm=T) ) , by = gene_symbol, .SDcols = quant_cols  ]
      dt_sub_aggr = dt_binned_sub[, lapply(.SD, function(x) sum(x,na.rm=T) ) , by = gene_symbol, .SDcols = quant_cols  ]
      Sys.time()


      df_sub_aggr = data.frame(cell = cell, dt_sub_aggr)
      Sys.time()



      head(df_sub_aggr)
      #df_sub_aggr[1:5, 1:5]
      #imp_list[[cell]] = data.frame(cell = cell, df_sub_aggr)

      temp_list[[i]] = df_sub_aggr

      end = Sys.time()
      print(end - start)

    }#for

    print("Combining all the cells")
    df_meta = do.call(rbind, temp_list)
    end = Sys.time()
    print(end - start)

  }#else

  print('------------------------------------------------')
  print(dim(df_meta))
  print(df_meta[1:5, 1:5])
  print('------------------------------------------------')



  return(df_meta)

}


compute_meta_cells_expr <- function(exp_matrix, pca, num_neighbors, aggr_type = 'sum', run_parallel = F)
{
  pca_distance_matrix <- compute_pca_distance_matrix(pca, distance_dim_count = 10)

  cells = colnames(pca_distance_matrix)
  head(cells)


  if(run_parallel)
  {
    df_meta_expr <- foreach(i=1:length(cells), .combine=cbind) %dopar%
      #df_meta <- foreach(i=1:10, .combine=cbind) %dopar%
      #df_meta <- foreach(i=1:4, .combine=rbind) %dopar%
    {
      library(data.table)
      library(rdist)

      #counter = counter + 1
      cell = as.character(cells[i])

      print(paste(i, cell))

      if(  sum(rownames(pca_distance_matrix) == cell) > 0 )
      {
        closest_neighbors = names(head(sort(pca_distance_matrix[cell, ]), (num_neighbors + 1) ))

        expr_sub = exp_matrix[, closest_neighbors]
        if(aggr_type == 'mean')
        {
          meta_expr = rowMeans(expr_sub)
        }
        else if(aggr_type == 'sum')
        {
          meta_expr = rowSums(expr_sub)
        }else
        {
          stop('aggr_type can be either sum or mean!')
        }

        head(meta_expr)
        length(meta_expr)

        result = meta_expr
      }else
      {
        result = exp_matrix[,cell]
      }#else
      df_result = data.frame(result)
      colnames(df_result) = c(cell)

      df_result

    }#foreach

  }else
  {
    temp_list = list()
    for(i in 1:length(cells))

    {


      #counter = counter + 1
      cell = as.character(cells[i])

      print(paste(i, cell))

      if(  sum(rownames(pca_distance_matrix) == cell) > 0 )
      {
        closest_neighbors = names(head(sort(pca_distance_matrix[cell, ]), (num_neighbors + 1) ))

        expr_sub = exp_matrix[, closest_neighbors]
        if(aggr_type == 'mean')
        {
          meta_expr = rowMeans(expr_sub)
        }
        else if(aggr_type == 'sum')
        {
          meta_expr = rowSums(expr_sub)
        }else
        {
          stop('aggr_type can be either sum or mean!')
        }

        head(meta_expr)
        length(meta_expr)

        result = meta_expr
      }else
      {
        result = exp_matrix[,cell]
      }#else
      df_result = data.frame(result)
      colnames(df_result) = c(cell)

      temp_list[[cell]] = df_result

    }#for

    df_meta_expr = do.call(rbind, temp_list)

  }#else

  return(df_meta_expr)

}#compute_meta_cells_expr


compute_meta_bin_counts <- function(df_binned, pca, distance_dim_count = 10, num_neighbors = 20, run_parallel = F)
{
  print(dim(pca))

  distance_dim_count = min(distance_dim_count, ncol(pca))

  pca_distance_matrix <- compute_pca_distance_matrix(pca, distance_dim_count = 10)
  df_binned_meta <- compute_meta_cells_methyl(df_binned, pca_distance_matrix, num_neighbors, run_parallel)

  print('***************************************')
  print(dim(df_binned_meta))
  print('***************************************')

  met_meta_sums <- rowSums(df_binned_meta[, 3:ncol(df_binned_meta)])
  df_binned_sum = data.frame(df_binned_meta[, 1:2], count = met_meta_sums)
  head(df_binned_meta)
  head(df_binned_sum)
  df_sum_matrix = dcast(df_binned_sum, gene_symbol ~ cell)
  df_sum_matrix[1:5, 1:5]
  rownames(df_sum_matrix) = df_sum_matrix$gene_symbol
  df_sum_matrix$gene_symbol = NULL


  met_meta_sums <- rowSums(df_binned_meta[, (3+5):(ncol(df_binned_meta)- 5 )])
  df_binned_sum = data.frame(df_binned_meta[, 1:2], count = met_meta_sums, run_parallel)
  head(df_binned_meta)
  head(df_binned_sum)
  df_sum_matrix_5 = dcast(df_binned_sum, gene_symbol ~ cell)
  df_sum_matrix_5[1:5, 1:5]
  rownames(df_sum_matrix_5) = df_sum_matrix_5$gene_symbol
  df_sum_matrix_5$gene_symbol = NULL

  list_result = list('df_binned_meta' = df_binned_meta, 'df_sum_matrix' = df_sum_matrix, 'df_sum_matrix_5' = df_sum_matrix_5)
  return(list_result)
}



compute_meta_cells <- function(df_met, df_demet, df_exp = NULL, exp_aggr_type = 'sum', run_parallel = F,  num_neighbors = 20)
{

  all_cells = unique(df_met$cell)


  df_met     = df_met[df_met$cell %in% all_cells, ]
  df_demet   = df_demet[df_demet$cell %in% all_cells, ]


  dim(df_met)
  dim(df_demet)

  print('Computing PCA')

  pca <- compute_pca(df_met, df_demet)
  pca_coords = pca$x


  print('PCA computed')

  met_meta_result = compute_meta_bin_counts(df_binned = df_met, pca = pca_coords,
                                            distance_dim_count = 10,  num_neighbors = num_neighbors,
                                            run_parallel)

  df_met_counts_meta = met_meta_result$df_binned_meta
  df_met_sums_meta = met_meta_result$df_sum_matrix
  df_met_sums_meta_5 = met_meta_result$df_sum_matrix_5

  print('Met. computed')

  demet_meta_result = compute_meta_bin_counts(df_binned = df_demet, pca = pca_coords,
                                              distance_dim_count = 10,  num_neighbors = num_neighbors,
                                              run_parallel)

  df_demet_counts_meta = demet_meta_result$df_binned_meta
  df_demet_sums_meta = demet_meta_result$df_sum_matrix
  df_demet_sums_meta_5 = demet_meta_result$df_sum_matrix_5

  print('De-Met. computed')


  df_rate_meta = df_demet_counts_meta
  n = ncol(df_rate_meta)
  df_rate_meta[,3:n] <- df_met_counts_meta[,3:n] / (df_met_counts_meta[,3:n] + df_demet_counts_meta[,3:n] )
  head(df_rate_meta)

  df_rate_meta_mean = df_met_sums_meta / (df_met_sums_meta + df_demet_sums_meta)
  df_rate_meta_mean[1:5, 1:5]

  df_rate_meta_mean_5 = df_met_sums_meta_5 / (df_met_sums_meta_5 + df_demet_sums_meta_5)
  df_rate_meta_mean_5[1:5, 1:5]

  if(is.null(df_exp))
  {
    df_exp_meta = data.frame()
  }else
  {
    df_exp_meta <- compute_meta_cells_expr(df_exp, pca_coords, num_neighbors = num_neighbors, aggr_type = exp_aggr_type, run_parallel = F)
  }



  meta_object = new("class_meta_cells"
                    , df_met_counts = df_met
                    , df_demet_counts = df_demet
                    , df_rate_meta = df_rate_meta
                    , df_rate_meta_mean = df_rate_meta_mean
                    , df_rate_meta_mean_5 = df_rate_meta_mean_5
                    , df_exp_meta = df_exp_meta
                    , num_neighbors = num_neighbors
                    , pca_coords = pca_coords)


  return(meta_object)


}



compute_binned_met_counts <- function(cov_dir, annot_file, sample_name = 'combined')
{

  cl <- makeCluster(24, outfile="", type = 'SOCK')
  registerDoSNOW(cl)


  if(sample_name == "combined")
  {
    cov_files = list.files(cov_dir, pattern = paste0('',  '.*cov.gz'))
  }else{
    cov_files = list.files(cov_dir, pattern = paste0(sample_name,'.*cov.gz'))

  }
  length(cov_files)

  dt_body = read.table(annot_file)


  colnames(dt_body) = c('chr', 'start', 'end',  'strand',
                        #'transcript_id',
                        'gene_id',
                        'gene_symbol', 'gene_type')

  #head(dt_body)

  dim(dt_body)
  length(unique(dt_body$gene_id))


  #dt_body = dt_body[dt_body$start > 15000, ]
  dt_tss_flank = dt_body[dt_body$start > 5000, ]


  dt_tss_flank$tss = dt_tss_flank$start

  strand_info = as.character(dt_tss_flank$strand)
  neg_strand = which(strand_info == '-')
  dt_tss_flank$tss[neg_strand] <- dt_tss_flank$end[neg_strand]

  dt_tss_flank$start = dt_tss_flank$tss - 5000
  dt_tss_flank$end = dt_tss_flank$tss + 4999


  head(dt_tss_flank)

  dt_tss_flank_coding = dt_tss_flank[dt_tss_flank$gene_type == 'protein_coding', ]
  head(dt_tss_flank_coding)

  bed_gr_body <- with(dt_tss_flank_coding, GRanges(chr, IRanges(start+1, end), strand, tss,
                                                   gene_id, gene_symbol, gene_type))
  bed_gr_body
  #bin_size = 200
  bin_size = 500


  shuf = F

  # i = 1

  counter = 0


  #df_binned <- foreach(i=1:8, .combine=rbind) %dopar%
  df_list_met = list()
  df_list_demet = list()

  result_list = list()

  result_list <- foreach(i=1:length(cov_files)) %dopar%
  #for(i in 1:length(cov_files))
  {
    library(data.table)
    library(GenomicRanges)
    #source('/mnt/isilon/tan_lab/uzuny/projects/MethylPredict/software_package/R_package/v12/MethylPredict/R/set_ops.R')

    print(i)
    cov_file = cov_files[i]
    print(cov_file)
    dt_cov = fread(paste0(cov_dir, cov_file) )
    colnames(dt_cov) = c('chr', 'start', 'end', 'met_rate', 'met', 'demet')
    #dt_cov$chr = paste0('chr', dt_cov$chr)

    head(dt_cov)
    dim(dt_cov)
    unique(dt_cov$chr)
    chromes_with_underscore = grepl('_', dt_cov$chr)
    head(dt_cov[chromes_with_underscore, ])
    regular_chromes = !grepl('_', dt_cov$chr)
    dt_cov_regular = dt_cov[regular_chromes, ]
    head(dt_cov_regular)
    dt_cov_regular = dt_cov_regular[dt_cov_regular$chr != 'chrM', ]
    dt_cov_regular = dt_cov_regular[dt_cov_regular$chr != 'chrMT', ]

    unique(dt_cov_regular$chr)


    bed_gr_cov <- with(dt_cov_regular, GRanges(chr, IRanges(start+1, end), strand = '*', met_rate, met, demet)  )
    bed_gr_cov

    print("intersecting")

    dt_inter = data.table(intersect_bed(bed_gr_body, bed_gr_cov))
    dim(dt_inter)
    head(dt_inter)

    if( nrow(dt_inter) == 0 )
    {
      NULL
    }else
    {

      quant_cols = c('b_met', 'b_demet')
      head(quant_cols)

      dt_aggr_by_gene <- dt_inter[, lapply(.SD, sum), by = .(gene_symbol, gene_type), .SDcols = quant_cols  ]
      head(dt_aggr_by_gene)


      attach(dt_inter)
      tss_distance = b_start - tss
      strand_info = as.character(dt_inter$strand)
      tss_distance[strand_info == '-'] = -1 * tss_distance[strand_info == '-']
      detach(dt_inter)


      head(tss_distance)
      dt_inter$tss_distance = tss_distance
      dt_inter$binned_distance = floor(tss_distance / bin_size)

      max(dt_inter$binned_distance)
      min(dt_inter$binned_distance)

      head(dt_inter)
      tail(dt_inter)

      print('Quantifying')


      dt_aggr <- dt_inter[, lapply(.SD, sum), by = .(gene_symbol, binned_distance, gene_type), .SDcols = quant_cols  ]
      head(dt_aggr)
      print('*****************I computed dt_aggr****************************************')
      dt_aggr$met_rate = dt_aggr$b_met / (dt_aggr$b_met + dt_aggr$b_demet) * 100
      head(dt_aggr)

      colnames(dt_aggr) = c('gene_symbol', 'binned_distance', 'gene_type', 'met_count', 'demet_count',  'met_rate')


      metrics = c('met_count', 'demet_count')

      count_list = list()

      for(metric in metrics)
      {

        metric_temp = metric
        dt_features = dt_aggr[, c('gene_symbol', 'binned_distance', ..metric_temp )]
        head(dt_features)

        dt_features_denorm = dcast(dt_features, formula = gene_symbol~binned_distance)
        head(dt_features_denorm)
        dim(dt_features_denorm)

        cell = cov_file

        cell = gsub('.cov.gz.sorted', '', cell)
        cell = gsub('.cov.gz', '', cell)
        cell = gsub('CpG_context.', '', cell)
        cell = gsub('_', '___', cell)


        df_denorm = data.frame(cell = cell, data.frame(dt_features_denorm))

        head(df_denorm)
        dim(df_denorm)
        count_list[[metric]] = df_denorm

      }#for(metric in metrics)

      cell_name = gsub('cov.gz', '', cov_file)
      #result_list[[cell_name]] = count_list

      count_list

    }#else



  }#result_list

  length(result_list)
  names(result_list)
  names(result_list[[1]])

  result_list[[1]][['met_count']] [1:5, 1:5]
  result_list[[1]][['demet_count']] [1:5, 1:5]

  #length(df_list)

  met_count_list = lapply( result_list , "[[" , 'met_count' )
  demet_count_list = lapply( result_list , "[[" , 'demet_count' )

  df_binned_met = do.call("rbind", met_count_list)
  dim(df_binned_met)
  head(df_binned_met)
  df_binned_met[is.na(df_binned_met)] = 0
  rownames(df_binned_met) = paste(df_binned_met$cell, df_binned_met$gene_symbol, sep = '_x_')
  head(df_binned_met)

  df_binned_demet = do.call("rbind", demet_count_list)
  dim(df_binned_demet)
  head(df_binned_demet)
  df_binned_demet[is.na(df_binned_demet)] = 0
  rownames(df_binned_demet) = paste(df_binned_demet$cell, df_binned_demet$gene_symbol, sep = '_x_')
  head(df_binned_demet)

  stopCluster(cl)

  binned_list = list('df_binned_met' = df_binned_met, 'df_binned_demet' = df_binned_demet)

  return(binned_list)
}








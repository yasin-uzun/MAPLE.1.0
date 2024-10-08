library(data.table)
library(GenomicRanges)


source('/mnt/isilon/tan_lab/uzuny/scripts/bed_ops/set_ops.R')


#dataset = 'scmandt_muscle'
#sample_name = 'GSM'

#dataset = 'scnmt'
#sample_name = 'combined'

#dataset = 'scnmt_gastr'
#sample_name = 'combined'

dataset = 'scnmt_gastr'
sample_name = 'combined'


cov_dir = paste0('public_met/',dataset,'/data/cov_files/')

binned_met_dir = paste0('public_met/',dataset,'/data/binned_met/')

dir.create(binned_met_dir, showWarnings = F)


if(organism == 'mouse')
{
  body_file = 'mouse/mm10/annotation/gencode.vM22.genes.protein_coding.bed'
}

if(organism == 'human')
{
  body_file = 'human/hg38/annotation/gencode.v29.genes.protein_coding.bed'
}


library(doSNOW)
cl <- makeCluster(24, outfile="", type = 'SOCK')
registerDoSNOW(cl)

#for(sample_name in sample_names)
#{

if(sample_name == "combined")
{
  cov_files = list.files(cov_dir, pattern = paste0('',  '.*cov.gz'))
}else{
  cov_files = list.files(cov_dir, pattern = paste0(sample_name,'.*cov.gz'))

}
length(cov_files)


dt_body = read.table(body_file)


colnames(dt_body) = c('chr', 'start', 'end',  'strand',
                      'gene_id',
                      'gene_symbol', 'gene_type')

#head(dt_body)

dim(dt_body)
length(unique(dt_body$gene_id))


dt_tss_flank = dt_body[dt_body$start > 5000, ]

dt_tss_flank$tss = dt_tss_flank$start
dt_tss_flank$tss[dt_tss_flank$strand == '-'] <- dt_tss_flank$end[dt_tss_flank$strand == '-']

dt_tss_flank$start = dt_tss_flank$tss - 5000
dt_tss_flank$end = dt_tss_flank$tss + 4999


head(dt_tss_flank)

dt_tss_flank_coding = dt_tss_flank[dt_tss_flank$gene_type == 'protein_coding', ]
head(dt_tss_flank_coding)

bed_gr_body <- with(dt_tss_flank_coding, GRanges(chr, IRanges(start+1, end), strand, tss,
                                                 gene_id, gene_symbol, gene_type))
bed_gr_body
bin_size = 500


shuf = F

counter = 0

shuf = F

counter = 0

df_list_met = list()
df_list_demet = list()


result_list <- foreach(i=1:length(cov_files)) %dopar%
{
  library(data.table)
  library(GenomicRanges)

  print(i)
  cov_file = cov_files[i]

  dt_cov = fread(paste0(cov_dir, cov_file) )
  colnames(dt_cov) = c('chr', 'start', 'end', 'met_rate', 'met', 'demet')

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
    tss_distance[strand == '-'] = -1 * tss_distance[strand == '-']
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


      counter = counter + 1

      cell = cov_file

      cell = gsub('.cov.gz.sorted', '', cell)
      cell = gsub('.cov.gz', '', cell)
      cell = gsub('CpG_context.', '', cell)

      cell = gsub('_', '___', cell)


      df_denorm = data.frame(cell = cell, data.frame(dt_features_denorm))

    
      head(df_denorm)
      dim(df_denorm)
      count_list[[metric]] = df_denorm


    }

    count_list

  }
}



length(result_list)
names(result_list)
names(result_list[[1]])

result_list[[1]][['met_count']] [1:5, 1:5]
result_list[[1]][['demet_count']] [1:5, 1:5]



met_count_list = lapply( result_list , "[[" , 'met_count' )
demet_count_list = lapply( result_list , "[[" , 'demet_count' )

df_binned_met = do.call("rbind", met_count_list)
dim(df_binned_met)
head(df_binned_met)
df_binned_met[is.na(df_binned_met)] = 0
rownames(df_binned_met) = paste(df_binned_met$cell, df_binned_met$gene_symbol, sep = '_x_')
head(df_binned_met)
binned_met_rds = paste0(binned_met_dir, '/binned_met.bin_size_',bin_size, '.', 'met_count', '.',sample_name,'.rds')
saveRDS(df_binned_met, file = binned_met_rds)

df_binned_demet = do.call("rbind", demet_count_list)
dim(df_binned_demet)
head(df_binned_demet)
df_binned_demet[is.na(df_binned_demet)] = 0
rownames(df_binned_demet) = paste(df_binned_demet$cell, df_binned_demet$gene_symbol, sep = '_x_')
head(df_binned_demet)
binned_demet_rds = paste0(binned_met_dir, '/binned_demet.bin_size_',bin_size, '.', 'demet_count', '.',sample_name,'.rds')
saveRDS(df_binned_demet, file = binned_demet_rds)



stopCluster(cl)


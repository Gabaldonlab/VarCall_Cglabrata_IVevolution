#!/gpfs/projects/bsc40/mschikora/anaconda3/envs/VarCall_CNV_env/bin/Rscript

# define environment
library(argparser, quietly=TRUE)
library(RSVSim, quietly=TRUE) # doc in https://rdrr.io/bioc/RSVSim/man/simulateSV.html

# parse cmd line args
argp = arg_parser("Takes known SV files and benchmarks some SVs of interest. It writes a file with the accuracy of each SV and some statistics like TP, FP, TN ... This only runs for those predicted SVs provided as arguments")

argp = add_argument(argp, "--known_tandemDuplications", help="path to the known SV as output by simulateSV")
argp = add_argument(argp, "--known_deletions", help="path to the known SV as output by simulateSV")
argp = add_argument(argp, "--known_insertions", help="path to the known SV as output by simulateSV")
argp = add_argument(argp, "--known_translocations", help="path to the known SV as output by simulateSV")
argp = add_argument(argp, "--known_inversions", help="path to the known SV as output by simulateSV")

argp = add_argument(argp, "--predicted_tandemDuplications", help="path to a bed file with the predicted tandemDuplications to benchmark")
argp = add_argument(argp, "--predicted_deletions", help="path to a bed file with the predicted deletions to benchmark")
argp = add_argument(argp, "--predicted_insertions", help="path to the known SV as output by simulateSV")
argp = add_argument(argp, "--predicted_translocations", help="path to the known SV as output by simulateSV")
argp = add_argument(argp, "--predicted_inversions", help="path to the known SV as output by simulateSV")

argp = add_argument(argp, "--outfile_benchmark", help="a file where to write the table of the benchmark.")
argp = add_argument(argp, "--tolerance_bp", default=50, help="The number of tolerance bp when comparing SVs")

argv = parse_args(argp)

# edit the argv (ONLY FOR TESTING ON SPECIFIC FILES)
#argv$known_tandemDuplications = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/tandemDuplications.tab"
#argv$predicted_tandemDuplications = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/benchmark_GridssClove_haploid/benchmark_max90x_ignoreRegionsTrue/several_parameter_combinations_filter_small_af1.00/filters_0/maxDELpctCovered10_minDUPcoverage7.tandemDuplications.bed"
#argv$known_deletions = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/deletions.tab"
#argv$predicted_deletions = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/benchmark_GridssClove_haploid/benchmark_max90x_ignoreRegionsTrue/several_parameter_combinations_filter_small_af1.00/filters_0/maxDELpctCovered10_minDUPcoverage7.deletions.bed"
#argv$outfile_benchmark = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/benchmark_GridssClove_haploid/benchmark_max90x_ignoreRegionsTrue/several_parameter_combinations_filter_small_af1.00/filters_0/predictions_benchmark.tab"
#argv$known_inversions = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/inversions.tab"
#argv$predicted_inversions = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/benchmark_GridssClove_haploid/benchmark_max90x_ignoreRegionsTrue/several_parameter_combinations_filter_small_af1.00/filters_0/complexVariants.inversions.bed"
#argv$known_translocations = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/translocations.tab.corrected"
#argv$predicted_translocations = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/benchmark_GridssClove_haploid/benchmark_max90x_ignoreRegionsTrue/several_parameter_combinations_filter_small_af1.00/filters_0/complexVariants.translocations.bedpe.withBlancedINFO"
#argv$known_insertions = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/insertions.tab.corrected"
#argv$predicted_insertions = "/home/mschikora/samba/scripts/VarCall_CNV_ReadProcessing/test_dir/RUN1_CST34_2G_FLZ_VarCallresults/gridss_output/simulation_1/uniform/benchmark_GridssClove_haploid/benchmark_max90x_ignoreRegionsTrue/several_parameter_combinations_filter_small_af1.00/filters_0/complexVariants.insertions.bedpe.withCopiedINFO"

# define functions
get_benchmark_array = function(pred_SV_df, known_SV_df){
  
  # takes dataframes of predicted and known dfs, and it returns a list of precision, recall, TP, FP, FN
  
  # perform comparison
  comp = compareSV(pred_SV_df, known_SV_df, tol=argv$tolerance_bp)
  
  # define the overlap field if not existing (for insertions)
  if (!"Overlap" %in% colnames(comp)){ comp$Overlap = comp$Overlap_5prime }

  # define the events that are found and missed
  found_events = paste(unique(comp[comp$Overlap!="",]$Name),collapse="||")
  missed_events = paste(unique(comp[comp$Overlap=="",]$Name),collapse="||")
  
  # define the true positive events
  overlapping_coordinates = unique(unlist(lapply(comp[comp$Overlap!="",]$Overlap, function(x) strsplit(x,", "))))
  
  # for each type of file, make a df that has the coordinate and the positive calls
  if ("chr2" %in% colnames(pred_SV_df)){
    df_1 = pred_SV_df[,c("chr", "start1", "end1")]; colnames(df_1) = c("chr", "start", "end")
    df_2 = pred_SV_df[,c("chr2", "start2", "end2")]; colnames(df_2) = c("chr", "start", "end")
    
    # add the ID
    df_1$ID = rownames(df_1); df_2$ID = rownames(df_2)
    df_calls = rbind(df_1, df_2)
    
  } else{
    df_calls = pred_SV_df
    df_calls$ID = rownames(df_calls)
  }

  # add the loc in text format
  df_calls$location = apply(df_calls, 1, function(r) sprintf("%s:%d-%d", r["chr"], as.integer(r["start"]), as.integer(r["end"])))
  
  # define the true positive IDs
  true_positive_events = unique(df_calls[df_calls$location%in%overlapping_coordinates, "ID"])
  true_positive_events_str =  paste(true_positive_events, collapse="||")
  
  # define the false positive IDs
  false_positive_events = setdiff(rownames(pred_SV_df), true_positive_events)
  false_positive_events_str = paste(false_positive_events, collapse="||")

  # get statistics
  all_predicted = length(row.names(pred_SV_df))
  all_known = length(row.names(comp))
  
  TP = sum(comp$Overlap!="")
  FP = all_predicted - TP
  FN = all_known - TP
  
  precision = TP/(TP + FP)
  recall = TP/(TP + FN)

  benchmark_array = c(precision, recall, as.integer(TP), as.integer(FP), as.integer(FN), as.integer(length(true_positive_events)), as.integer(length(false_positive_events)), as.integer(all_known), found_events, missed_events, true_positive_events_str, false_positive_events_str)
  
  return(benchmark_array)
  
}

# initialize_df_all, which will be constantly added with each of the types of precision and recall values
df_benchmark = data.frame()

# go through each SVTYPE
for (svtype in c("deletions", "tandemDuplications", "inversions", "insertions", "translocations")) {
  
  # define the filenames
  predicted_filename = argv[paste("predicted", svtype ,sep="_")][[1]]
  known_filename = argv[paste("known", svtype ,sep="_")][[1]]
  
  # continue if it is NA
  if (is.na(predicted_filename)){next} 
  
  # simple events
  if (svtype %in% c("deletions", "tandemDuplications", "inversions")) {
    
    # get the predicted SV into a df
    pred_SV_df = read.table(predicted_filename, header=TRUE)
    row.names(pred_SV_df) = pred_SV_df[,"ID"]
    pred_SV_df = pred_SV_df[, c("chr", "start", "end")]
  
    # get the known SV into a df
    known_SV_df = read.table(known_filename, header=TRUE)
    
    # perform prediction
    benchmark_array = as.data.frame(t(c(get_benchmark_array(pred_SV_df, known_SV_df), svtype)))
    
    # put into df  
    df_benchmark = rbind(df_benchmark, benchmark_array)
    
  
  }
  
  # complex events
  if (svtype %in% c("insertions", "translocations")) {
    
    # define the differential field, which is TRUE/FALSE
    if (svtype=="insertions"){dif_field="Copied"}
    if (svtype=="translocations"){dif_field="Balanced"}
    
    # load the dataframes
    all_pred_SV_df = read.table(predicted_filename, header=TRUE)
    all_known_SV_df = read.table(known_filename, header=TRUE)
    
    # go through each of the values of dif_field (TRUE/FALSE)
    for (value_dif_field in unique(all_known_SV_df[,dif_field])) {
      
      # define the type of event
      svtype_dif_field = paste(svtype, dif_field, value_dif_field, sep="_")
      
      # get the dfs with this value
      pred_SV_df = all_pred_SV_df[all_pred_SV_df[dif_field]==value_dif_field,]
      row.names(pred_SV_df) = pred_SV_df[,"ID"]
      pred_SV_df = pred_SV_df[,c("chr", "start1", "end1", "chr2", "start2", "end2")]
      known_SV_df  = all_known_SV_df[all_known_SV_df[dif_field]==value_dif_field,]

      # perform prediction
      benchmark_array = as.data.frame(t(c(get_benchmark_array(pred_SV_df, known_SV_df), svtype_dif_field)))
      
      # put into df  
      df_benchmark = rbind(df_benchmark, benchmark_array)
      

    } 
  }
}

# change the column names
colnames(df_benchmark) = c("precision", "recall", "TP", "FP", "FN", "TP_number_eventIDs", "FP_number_eventIDs", "all_events", "found_events", "missed_events", "TP_eventIDs", "FP_eventIDs", "typeSV")

# write to outfile
write.table(df_benchmark, file=argv$outfile_benchmark, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



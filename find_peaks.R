# #####################################################################################
# Script by Martha Kandziora, 2020
# v.0.11
# 
# Find second peak of two overlapping distributions.
# Used to find the optimal divergence between sequences that are potential paralogs.
# Can be run with Rscript, set working dir to the folder that 
# contains a subfolder "dist" with all the distance matrices produced by Biopython distmat2
# #####################################################################################

setwd=("./")

library(mixtools)

# functions to use
run_mixEM_multiple <- function(data, repeats){
  all_md = list()
 
  for(i in 1:repeats){
    mD <- normalmixEM(data, maxit = 10000)
    all_md[[i]] <- mD
  }
  
  maxmd=-100000000000000000
  for(i in 1:repeats){
    entry = all_md[[i]]
    if(entry$loglik > maxmd){
      max_md = entry
    }
  }
  return(max_md)
}

get_second_peak <- function(all_da, name, repeats){
  max_md = run_mixEM_multiple(all_da, repeats)
  fn = paste(getwd(), paste(name,'png', sep = '.'), sep = '/')
  
  png(fn)
  plot(max_md, which=2)
  title(sub=name)
  dev.off()
  
  
  if(max_md$mu[1] > max_md$mu[2]){
    second_peak = max_md$mu[1]
  }else{
    second_peak = max_md$mu[2]
  }
 
  all_second_peaks = append(all_second_peaks, second_peak)
  return(all_second_peaks)
}


##########################


folder = "./dist"
files = list.files(folder)


all_second_peaks = c()
for (fn in files)
{
  fn_rel = paste(folder, fn, sep="/")
  # read the data
  tryCatch({
    all_da = read.csv(fn_rel, sep='\n', header=FALSE)
    all_da = unlist(all_da)
    locus = stringr::str_split(fn, '.fasta.mafft.trim.fas.dist.txt')[[1]][1]
    
    # turn values into bimodal and get min value of second peak
    tryCatch(
      {all_second_peaks = get_second_peak(all_da, locus, repeats = 10)},
      error = function(e) {
        print(paste('error:', e))
      }
    )
  },error = function(e) {
    print(paste('empty file?:', e))
  }
  )
}


# hist(all_second_peaks)
# standard_deviation = sd(all_second_peaks)
# standard_deviation
mean_val = round(mean(all_second_peaks), digits = 2)
mean_val
# median_val = median(all_second_peaks)
# median_val



mD <- run_mixEM_multiple(all_second_peaks, repeats = 10)

fn = "histogramm_secondpeak.png"
png(fn)
plot(mD, which=2)
title(sub="all second peaks")
dev.off()

if(mD$mu[1] > mD$mu[2]){
  second_peak = mD$mu[1]
}else{
  second_peak = mD$mu[2]
}

peak2 = round(second_peak)
write(c(mean_val, peak2), stdout())

fileConn <- file('paralog_divergence_range.txt')
writeLines(c(toString(mean_val), toString(peak2)), fileConn)
close(fileConn)


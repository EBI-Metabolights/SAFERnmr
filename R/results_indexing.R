# unzip_studies #########################################################################################################################
#' Unzips studies found in a subdirectory
#'
#' There are several operations that need to be done for SAFER results to be used
#' (i.e.) with browse_evidence().
#' 
#' Results are written as zip files for portability. 
#' 
#' The idea is that you have a local location where you keep your SAFER results, 
#' and re-run this function to unpack any new results. 
#' 
#' This command
#'  - unzips files which haven't previously been unzipped
#'  - excludes any studies you explicitly name
#'  - checks the heatmap html file to ensure it was written correctly, and 
#'    corrects it if needed
#' 
#' Examples:
#'  data.dir <- '/Users/mjudge/Documents/safer/results/' # this is where you put your zip files
#'  unzip_studies(data.dir) 
#'  unzip_studies(data.dir, exclude = c('1697751863','1697759923')) # you can also exclude problematic runs
#'
#' @param data.dir directory containing zip files of your run results Other files and unzipped results can be in there, but any zip file will be assumed to be a run result
#' @param exclude (optional) list of run names that are problematic (i.e. cannot be unzipped correctly for some reason, or which you don't want to unzip)
#' 
#' @return nothing
#'          
#' @importFrom magrittr %>%
#'
#' @export
unzip_studies <- function(data.dir, exclude = NULL){
  
    # If a study fails because of odd zipping, re-zip it in shell with the following:
      # zip -r -j file/path/output.zip file/path/input
    # Unzip any non-unzipped files ####
      # data.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/'
      data <- dir(data.dir)
        excluded <- which(grepl(exclude %>% paste(collapse = "|"), data))
        zipped <- grepl('.zip', data)
        not.zip <- which(!zipped)
        zipped <- which(zipped)
          zip.names <- data[zipped] %>% stringr::str_remove_all('.zip')
        already.unzipped <- (zip.names %in% data[not.zip]) %>% zipped[.]
        dont.unzip <- c(already.unzipped, excluded, not.zip)
        do.unzip <- data[-dont.unzip]
          
      if (length(do.unzip) > 0){
        fail <- lapply(do.unzip, function(d){
          # d <- data[-dont.unzip] %>% .[1]
          message('Unzipping ', d, ' ...')
          
          # New ####
            tryCatch({
              run.id <- d %>% stringr::str_remove('.zip$')
              out <- paste0(data.dir, '/',run.id) %>% stringr::str_replace_all('//','/')
              unzip(paste0(data.dir, '/',d) %>% stringr::str_replace_all('//','/'), junkpaths = FALSE, exdir = out)
              check_heatmap_html_output(out)
              return(0)
              
            }, error = function(cond){
              1
            })
        })
        
      } else {
        
        fail <- list('There are no zipped studies which need to be unzipped.')
      }

  if (is.numeric(fail)){
    return(fail %>% unlist)
  } else {
    return("nothing to unzip")
  }
  
}

# index_studies #########################################################################################################################
#' Indexes unzipped studies found in a subdirectory
#'
#' Sister function to unzip_studies()
#' 
#' After running unzip_studies(), it may be useful to get an overview of all the results in the directory.
#' This function compiles all run summaries into a single data frame, and saves the table in a csv file in the directory. 
#' 
#' Examples:
#'  run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923','1713439191',''))
#'
#' @param data.dir directory containing zip files of your run results
#' @param exclude (optional) list of run names that are problematic (i.e. cannot be unzipped correctly for some reason, or which you don't want to include in the index)
#' 
#' @return data frame containing the compiled run summary data for the studies in the directory. Also written to ./run.summary.csv.
#'          
#' @importFrom magrittr %>%
#'
#' @export
index_studies <- function(data.dir, exclude = NULL){
  
  # Extract out run summary data ####
  
      data <- dir(data.dir)
        dontuse <- grepl(paste('.zip$','.csv$','.pdf$','.xlsx$',exclude %>% paste(collapse = "|"), sep = "|"), data)
        unzipped <- !dontuse
      data.unzipped <- paste0(data.dir,data[unzipped])


    df <- lapply(data.unzipped, function(x){
        dat <- read.csv(paste0(x, '/run.summary.csv'))
        names(dat) <- names(dat) %>% stringr::str_replace_all('\\.', '_')
        dat$local_id <- dat$run_id
        dat
    }) %>% bind_rows
    
    # df$total_time[df$total_time > 10]  <- df$total_time[df$total_time > 10] / 60
    # df$matching_elapsed_time[df$matching_elapsed_time > 10]  <- df$matching_elapsed_time[df$matching_elapsed_time > 10] / 60
    df$local_path <- data.unzipped
    df <- df[,!(names(df) == 'X')]
    df$start <- df$run_id %>% as.POSIXct(origin = "1970-01-01")
    df
}


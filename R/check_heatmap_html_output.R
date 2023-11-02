#' check_heatmap_html_output
#'
#' Check the provided results file for correct standalone html heatmap output.
#' 
#' If file is a zip, unzip first to a temp directory, replace the file, 
#' then re-zip and remove the temp directory. 
#' 
#' If file is a directory, simply re-write the html heatmap file as a standalone.
#' 
#'
#' @param pars a list containing necessary parameters for the function
#' @param copy.to (optional) full filepath location where you want the html file copied to
#' @return RDS file containing match scores between library reference spectra and the best subset spectrum score for each.
#' 
#' @importFrom magrittr %>%
#' @importFrom htmlwidgets saveWidget
#'
#' @export
check_heatmap_html_output <- function(filename, copy.to = NULL){
  # filename <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1698408936'
  # filename <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1698599655.zip'
  
  # Accessory functions ####
    can_pandoc <- function(){
      
      # Check if the system is Linux
        is_linux <- tolower(Sys.info()["sysname"]) == "linux"
        
      if (is_linux){
        
      
      # Check if pandoc is installed using system command
        # pandoc_check <- system("pandoc --version", intern = TRUE)
        
        # Check if pandoc is installed using system2 command
        pandoc_check <- system2("pandoc", "--version", stdout = TRUE, stderr = TRUE)
        
        # Check the output for pandoc version
        
        return(length(pandoc_check) > 0)
        
      } else {
        
        # If not linux, assume pandoc functionality works:
        
          return(TRUE)
        
      }
    }
    
  # Main function ####
  
    looks.like.zip <- filename %>% stringr::str_detect('\\.zip$')
      
    if (looks.like.zip){
      
      # Open up the zipfile
        tmp <- filename %>% dir_pop %>% paste0('/safer_tmp_', Sys.time() %>% as.numeric %>% round, '/')
        
        dir.create(tmp, showWarnings = FALSE)
        
        # unzipped.file <- filename %>% stringr::str_remove('\\.zip$')
        unzipped.file <- tmp
        
        unzip(filename, exdir = unzipped.file)
      
    } else {
      
      unzipped.file <- filename
      
    }
  
    # Open the heatmap.file and check to see if it has the correct header for 
    # self-contained html file
    
      heatmap.file <- paste0(unzipped.file,'/match_scores_sample_x_compound.html') %>% stringr::str_replace_all('//', '/') 
      # - same as printing fcn uses
      
      first.line <- readLines(heatmap.file, n = 1)
      
      file.is.html <- first.line == '<!DOCTYPE html>'
      
      if (file.is.html){
        
        # do nothing
        written <- FALSE
        
      } else {
  
        if (can_pandoc()){
          
          # Remake the heatmap and print
          
            message('Printing self-contained heatmap html file to: ', heatmap.file)
            n.cmpds <- print_html_heatmap_from_scoresfile(unzipped.file)
            written <- TRUE
            
          # If requested, copy the file to alternate location
          
            if (!is.null(copy.to)){
              
              message('Also copying file to: ', copy.to)
              file.copy(from = heatmap.file, 
                        to = paste0(copy.to), 
                        overwrite = TRUE)
            }
          
        } else {
          
          message("\nSelf-contained heatmap html file cannot be drawn because this is a Linux system which does not have pandoc installed, and that is a system dependency for htmlwidgets::saveWidget(), and therefore print_html_heatmap_from_scoresfile().\n")
          written <- FALSE
          
        }

      }
    
      if (looks.like.zip){
        
        if (written){
          
          # zip up the tmp file
          
            # zip(zipfile = '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1698599655_2.zip', #filename, 
            #     files = tmp)
            
        }
        
        # Either way, cleanup the tmp file
         
          unlink(tmp, recursive = TRUE)
          
      }
      
      
}

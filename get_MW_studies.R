# Get mwtab text files for NMR-only studies 
# For each study where 
# NMR is the sole instrument used 
# and where the species is a mammal 
# and where submitted before July 1st 2022.

  library(rvest)
  library(dplyr)
  
  # Get the list of urls to hit
  
  library(xml2)
  library(RCurl)
  library(dplyr)
  library(rvest)

  base.url <- "https://www.metabolomicsworkbench.org/data/"
  exclude.after <- "2022-07-01"
  file.dest <- "/Users/mjudge/Documents/metabolomics_workbench/"
  
  # *** Sort by Analysis, then save the 2000-row html page (Save As html) and provide:
  summary.table.file <- "/Users/mjudge/Desktop/Metabolomics Workbench _ NIH Data Repository _ Overview.html"
  
##########################################################################################      
    # Load the table ####
      html <- read_html(summary.table.file)
      table <- html %>%
          html_nodes(xpath = "//*[@id=\"content\"]/table[2]") %>%
          html_table(header=T) %>% bind_rows
      # table$Analysis %>% unique
    
    # Filter for NMR ####
    is.nmr <- table$Analysis %in% c("NMR", "NMR#")
      # table$Species[is.nmr] %>% tolower %>% unique
    
    # Filter for mammals ####
      genuses <- c("Homo ",
                   "mus ",
                   "sus ",
                   "equus ",
                   "felis ",
                   "rattus ",
                   "ovis ",
                   "hydrochoerus ",
                   "capra ",
                   "tursiops "
                   ) %>% paste(collapse = "|")
      
      is.mammal <- grepl(pattern = genuses, table$Species, ignore.case = TRUE)
    
    # Combine filters ####
      study.rows <- which(is.nmr & is.mammal)
      study.names <- table$`Study ID`[study.rows]
      
    # Get the urls for the studies ####
    table.nodes <-  lapply(study.rows, function(s){
                      # s <- study.rows[1]
                      selector <- paste0("#content > table:nth-child(4) > tbody", # the table should be pretty static.
                                         " > tr:nth-child(",s+1,")",             # table row (including header)
                                         " > td:nth-child(",1,") > a")           # table column (static)
                
                      # For some reason, the href selectors don't work with this page. Manually dig in:
                            return(html %>%
                                        html_nodes("#content > table") %>%
                                        .[[2]] %>%
                                        html_nodes("tr"))
                    }) %>% .[[1]] %>% .[-1]
      
   # for each desired row, get the link ####
        desired.col <- 1
      
        # rows s+1
        df.studies <- lapply(table.nodes[study.rows],
                             function(node)
                               {
                                  # n <- 2
                                  # node <- table.nodes[[n]]
                                  element.study <- node %>% html_nodes("td") %>% .[[desired.col]] %>% as.character
                                    study.url.part2 <- element.study %>% stringr::str_split(., '"') %>% .[[1]] %>% .[2]
                                    
                                  study.url <- paste0(base.url, study.url.part2) %>% stringr::str_replace_all("amp;","")
                                  # print(study.url)
                                  return(study.url)
                                }
                             ) %>% unlist %>% # put in data frame:
          data.frame("Study" = study.names,
                     "study.url" = .,
                     "species" = table$Species[study.rows],
                     "analysis" = table$Analysis[study.rows])
        
############################################################################################        
        
        # For each url, go in and get the date
        # - check to see if it's old enough
        # - if so, get the mwtab file page and the download links for mwtab files
        # - open each mwtab file
        # - if Analysis = NMR | Instrument Type [contains] "NMR", write the mwtab file
        # - if multiple files, indicate with a counter in the filename
        
        tab.file.dest <- paste0(file.dest, "mwtab_files_nmr/")
          dir.create(tab.file.dest)
            pblapply(1:nrow(df.studies),
            # pblapply(144,
                          function(s) 
                            {
                            # s <- 144
                            # Keep track of time between server pings
                              tictoc::tic()
                            # Read the webpage
                              st.html <- read_html(df.studies$study.url[s])
                              
                            # Get the submission date ####
                              st.str <- st.html %>% toString
                              sub.date <- extract.str.between(st.str,
                                                              preceded.by = "<b>Submit Date</b></td>\n<td>",
                                                              followed.by = "</td>") %>% as.Date
                              if (is.na(sub.date)){sub.date <- exclude.after %>% as.Date}
                              release.date <- extract.str.between(st.str,
                                                              preceded.by = "<b>Release Date</b></td>\n<td>",
                                                              followed.by = "</td>") %>% as.Date
                              
                            # If the study is old enough, go get its mwtab files ####
                              if (any(  c(sub.date,release.date) < as.Date(exclude.after) )){
                                # Derive the mwtab page urls using pattern ####
                                 res.type <- extract.str.between(df.studies$study.url[s],
                                                              preceded.by = 'ResultType=',
                                                              followed.by = '')
                                 mwtab.url <-  paste0(base.url,
                                                      "study_textformat_list.php?STUDY_ID=",
                                                      df.studies$Study[s])
                                 mwtab.html <- read_html(mwtab.url)
                                 mwtab.html.str <- mwtab.html %>% toString
                                 
                                 # Get any mwtab file download on the page  ####
                                   mwtab.link.nodes <- mwtab.html %>%
                                      html_nodes("a") %>% as.character
                                   
                                   # Get the second part of the url (specific to each mwtab file)  ####
                                   # - only use the view link so we get a text file output
                                     mwtab.tags <- lapply(mwtab.link.nodes, function(l){
                                     is.mwtab <- stringr::str_detect(toString(l), 'study_textformat_view')
                                     is.view <- stringr::str_detect(toString(l), '\">View</a>')
                                     if(is.mwtab & is.view){
                                       return(l %>% stringr::str_replace(., '<a href=\"','') %>%
                                                     stringr::str_replace(., '&amp','') %>%
                                                     stringr::str_replace(., ';','&') %>%
                                                     stringr::str_replace(., '\">View</a>','')
                                       )
                                     } 
                                   }) 
                                     
                                   # Build full URL  ####
                                     mwtab.links <- mwtab.tags %>% unlist %>% 
                                                     paste0(base.url, 
                                                            .)
                                     
                                   # Download the mwtab files for this study, and write if they are NMR files  ####
                                     mwtab.info <- 
                                     lapply(1:length(mwtab.links), 
                                            function(x) 
                                            {
                                              # x <- 1
                                               # Don't hammer the server ####
                                                 Sys.sleep(runif(1, min = crawl_delay_range[1], 
                                                            max = crawl_delay_range[2])/2)
                                               # Read the mwtab file from the web ####
                                                 str <- readr::read_file(mwtab.links[x])
                                               # Check to make sure this is an NMR file. If not, don't write the file ####
                                                 lns <- str %>% stringr::str_split(., '\\n') %>% unlist
                                                 
                                                 # Get the analysis and inst type lines
                                                   an.line <- grepl(pattern = "AN:ANALYSIS_TYPE", lns) %>% which %>% lns[.]
                                                   ins.line <- grepl(pattern = "INSTRUMENT_TYPE", lns) %>% which %>% lns[.]
                                                   ins.name <- grepl(pattern = "INSTRUMENT_NAME", lns) %>% which %>% lns[.]
                                                   
                                                 # # Trim the lines to get the types ####
                                                 #   label.end <- stringr::str_locate(an.line, "AN:ANALYSIS_TYPE\\s*") %>% as.numeric %>% .[2]
                                                 #   an.type <- stringr::str_sub(an.line, label.end, nchar(an.line)) %>% 
                                                 #                stringr::str_remove(., "\\t")
                                                 #   label.end <- stringr::str_locate(ins.line, "INSTRUMENT_TYPE\\s*") %>% as.numeric %>% .[2]
                                                 #   ins.type <- stringr::str_sub(ins.line, label.end, nchar(ins.line)) %>% 
                                                 #                stringr::str_remove(., "\\t")
                                                 print(c(s,x))
                                                   
                                               # If it is, write the file ####
                                               downloaded <- FALSE
                                               if (any(c(stringr::str_detect(tolower(an.line), tolower("NMR")),F)) |
                                                   any(c(stringr::str_detect(tolower(ins.name), tolower("NMR")),F)) |
                                                   any(c(stringr::str_detect(tolower(ins.line), tolower("NMR")),F)))
                                                 {
                                                    downloaded <- TRUE
                                                    write(str, 
                                                    paste0(tab.file.dest,
                                                            df.studies$Study[s],"_",x,".mwtab"))    
                                               }
                                              
                                               # Return info
                                                return(data.frame(Analysis.Type = an.line,
                                                                  Instrument.Type = ins.line,
                                                                  Instrument.Name = ins.name,
                                                                  downloaded = downloaded))
                                             }) %>% rbind
                              }
                              
                              # write(str, file = "/Users/mjudge/Desktop/myfile.txt")
                                   

                            # Don't hammer the server - add about a second, if needed ####
                              elapsed <- tictoc::toc(quiet = TRUE)
                              if (round(as.numeric(elapsed$toc - elapsed$tic)) < 1){
                                  Sys.sleep(runif(1, min = crawl_delay_range[1], 
                                                max = crawl_delay_range[2]))
                              }
                            })
        
  f.names <- dir(tab.file.dest, pattern = "*.mwtab")
  downloaded.study.names <- f.names %>% stringr::str_remove(., pattern = '_[:digit:].mwtab')
  extra <- which(!(downloaded.study.names %in% df.studies$Study))
  missing <- which(!(df.studies$Study %in% downloaded.study.names))

  # Report
    message(length(downloaded.study.names), " files downloaded.")
    if(any(missing)){
      message("Missing files: ", paste(c(df.studies$Study[missing]),collapse = " | "))
    } else {message("All files downloaded successfully.")}
    if(any(extra)){
      message("Extra files present in directory: ", paste(c(df.studies$Study[extra]),collapse = " | "))
      } else {message("No extra files present.")}
    
  df.studies$study.url[missing]

  
  
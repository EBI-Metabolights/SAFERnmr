devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')#1697793655
# pipeline('/Users/mjudge/Documents/ftp_ebi/params/params_1709742511.yaml')

# Accessory ####

data.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/'

    unzip_studies(data.dir, exclude = c('1697751863','1697759923'))
    run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923','1713439191','1726497509','1726501646'))

    
    
########### Just look at latest one ######
  
  browse_evidence('/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1719502487')
  browse_evidence(run.idx$local_path[nrow(run.idx)], select.sample = select.sample)
  browse_evidence(run.idx$local_path[nrow(run.idx)], select.compounds = 'Citrate')
  
  
########### Runs by end time ######

  run.idx <- run.idx %>% filter(write_time >= '2024-09-23 15:00:00') %>% arrange(write_time)
  i <- 4
  browse_evidence(run.idx$local_path[i], select.compounds = 'Citrate')

  
  


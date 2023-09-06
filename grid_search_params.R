# 
in.files <- c('/Users/mjudge/Desktop/MTBLS1_params.yaml',
              '/Users/mjudge/Desktop/MTBLS395_params.yaml',
              '/Users/mjudge/Desktop/MTBLS430_params.yaml',
              '/Users/mjudge/Desktop/MTBLS424_params.yaml')


adjust_params <- function(change.pars, in.files, extension){
  
  lapply(in.files, function(f){
    
    # Read the parameters in 
      pars <- yaml::yaml.load_file(paste0(results.dir,'/params.yaml'), eval.expr = TRUE)
    
    # Change the parameters
      
      
    # Write params.yaml file
    
      
  })
  yaml::write_yaml(pars, )
  
}
devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')#1697793655
# pipeline('/Users/mjudge/Documents/ftp_ebi/params/params_1709742511.yaml')

# Accessory ####

data.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/'

    unzip_studies(data.dir, exclude = c('1697751863','1697759923'))
    run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923','1713439191'))

par.templ <- '/Users/mjudge/Documents/GitHub/SAFERnmr/params.template.yaml'
pars <- yaml::yaml.load_file(paste0(par.templ), eval.expr = TRUE)

study.data <- readxl::read_xlsx('/Users/mjudge/Edison_Lab@UGA Dropbox/Michael Judge/MJ_UGA_Root/Scheduling/ml.studies.xlsx')

dir.create('output_params')

# Helper function to update a nested list using a path
update_nested_list <- function(params, yaml_keys, value) {
  if (length(yaml_keys) == 1) {
    params[[yaml_keys]] <- value  # Base case: if only one key, update the value
  } else {
    # Recursively go deeper into the structure
    key <- yaml_keys[1]
    params[[key]] <- update_nested_list(params[[key]], yaml_keys[-1], value)
  }
  return(params)
}

# Function to create params.yaml for each row in study.data using _ as a delimiter for YAML list levels
generate_params_yaml <- function(template_file, study_data, output_dir, filename_column) {
  
  # Read the template params.yaml with eval.expr = TRUE
  template <- yaml.load_file(template_file, eval.expr = TRUE)
  
  # Loop through each row in study_data
  for (i in 1:nrow(study_data)) {
    
    # Make a copy of the template for the current row
    params <- template
    
    # Loop through each column in study_data
    for (col_name in colnames(study_data)) {
      # Skip the filename column since it's used for the output filename
      if (col_name == filename_column) next
      
      # Get the YAML path from the column name (split _ into list levels)
      yaml_path <- strsplit(col_name, "_")[[1]]
      
      # Get the value from the current row
      value <- study_data[[col_name]][i]
      
      # Update the params using the yaml_path
      params <- update_nested_list(params, yaml_path, value)
    }
    
    # Use the value in filename_column to name the output YAML file
    output_file <- file.path(output_dir, study_data[[filename_column]][i])
    
    # Write the updated params to a new YAML file
    write_yaml(params, output_file)
    
    # Inform the user
    cat("Generated:", output_file, "\n")
  }
}

# Call the function
generate_params_yaml(par.templ, study.data, "output_params", "filename")

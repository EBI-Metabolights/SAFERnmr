#' browse_evidence : results viewer shiny app
#' 
#' Browse reference spectra (PRCSs) and metabolite-sample summed scores, and drill
#' down to the evidence for their associations at the multiple levels. 
#' 
#' To use, call with the results directory. 
#'
#' @param results.dir A list of input parameters.
#' @return spins up a shiny app from the results directory
#' @import yaml
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace
#' 
#' @import shiny
#' @import plotly
#' 
#' @export
browse_evidence <- function(results.dir = NULL, select.compounds = NULL, select.samples = NULL){
#######################################################################################
  if (is.null(results.dir)) {
    results.dir <- getwd()
  }
  
  # for the selectizer:
  sort.choices <- c("Scores - cluster", "Search for compound...") # "Compound names - cluster", 'Compound names - alphabetical', 
  
  
## Params ###############################################################################
  
  message('Evidence Viewer: reading data...')
  # Locate us ####
  
    # Handle whether or not user adds /
      if (stringr::str_sub(results.dir, start= -1) != "/"){
        results.dir <- paste0(results.dir, '/')
      }

  pars <- yaml::yaml.load_file(paste0(results.dir,'params.yaml'), eval.expr = TRUE)
  
  # Apply same validation function as in the pipeline (so pars are interpreted the same throughout)
    par.list <- pars %>% valid_pars(quiet = TRUE)
    pars <- par.list$pars
  
  study <- pars$study$id
  ppm.tolerance = pars$matching$filtering$ppm.tol
  hshift = 0
  

##########################     Setup/Read Data    ####################################

    # Read in spectral matrix data
        fse.result <- readRDS(paste0(results.dir, "fse.result.RDS"))
          xmat <- fse.result$xmat
          ppm <- fse.result$ppm
          rm(fse.result)

    # Read in the features 
    
        features.c <- readRDS(paste0(results.dir, "feature.final.RDS"))

    # Read in scores matrix 
      scores <- readRDS(paste0(results.dir,"scores.RDS"))
        scores.matrix <- scores$ss.ref.mat
        # colnames(scores.matrix) <- 1:ncol(scores.matrix)
        colnames(scores.matrix) <- colnames(scores.matrix) %>% lapply(function(x) strsplit(x, '\\s\\|\\s') %>% .[[1]] %>% .[c(1,3)] %>% paste(collapse = ' | ')) %>% unlist
          
          if (is.null(select.compounds)) {
            select.compounds <- 1:nrow(scores.matrix)
          }
        
          if (is.null(select.samples)){
            select.samples <- 1:ncol(scores.matrix)
          }
        
    # Read in match data
        backfit.results <- readRDS(paste0(results.dir,"smrf.RDS"))
          backfits <- backfit.results$backfits
          match.info <- backfit.results$match.info
          

    # Read in library data
        lib.data.processed <- readRDS(paste0(results.dir, "lib.data.processed.RDS"))
                
          # ppm isn't needed anymore; using spectral matrix ppm.
            
            lib.data.processed <- lib.data.processed %>% lapply(function(x) {x$data <- NULL; x$ppm <- NULL; return(x)})
            
          # add compound names as column in match.info
            compound.names.refmat <- lapply(lib.data.processed, function(x) x$compound.name) %>% unlist
            match.info$ref.name <- match.info$ref %>% compound.names.refmat[.]
          
      rfs.used <- scores$rfs.used
      
##########################     Filter    ####################################          
          
      # Remove all compounds without any scores > 0.5

        score.cutoff <- 0
        keeprefs <-  intersect(which(apply(scores.matrix, 1, max) > score.cutoff), select.compounds)
        
        keepsamples <-  which(apply(scores.matrix, 2, max) > score.cutoff)
          message('Out of ', ncol(scores.matrix), ' samples, ', length(keepsamples), ' passed max score cutoff of ', score.cutoff, '. ')
          
          # If select.samples is a sample name, then convert to an index
          if (is.character(select.samples)){
            select.sample.inds <- grepl(pattern = select.samples, colnames(scores.matrix)) %>% which
          } else {
            select.sample.inds <- select.samples
          }
          
          keepsamples <- intersect(keepsamples, select.sample.inds)
          if (length(keepsamples) == 0){
            error('The selected sample did not have any scores > ', score.cutoff, '. Quitting.')
          }
          
        refs.used <- keeprefs
        compound.names.refmat <- compound.names.refmat[refs.used]
        samples.used <- keepsamples
        
        scores.matrix <- scores.matrix[keeprefs,,drop=F]
        scores.matrix <- scores.matrix[,keepsamples,drop=F]
        lib.data.processed <- lib.data.processed[keeprefs]
        
        fits.keep <- match.info$ref %in% refs.used
        match.info <- match.info[fits.keep, ]
        backfits <- backfits[fits.keep]

##########################     Cluster    ####################################      

      
      if (nrow(scores.matrix) > 1){
        clust.refs <- T
        clust.samples <- T
      } else {
        clust.refs <- F
        clust.samples <- T
      }

      # hclust: column subsetting doesn't work when names are used. ### #
      ref.order <- 1:nrow(scores.matrix)
      sample.order <- 1:ncol(scores.matrix)
      
      if (clust.refs){
        if (nrow(scores.matrix) > 1){
          ref.order <- hclust(dist(scores.matrix))$order
        } else {
          ref.order <- 1
        }
          
      }
      
      if (clust.samples){
        if (ncol(scores.matrix) > 1){
          sample.order <- hclust(dist(t(scores.matrix)))$order
        } else {
          sample.order <- 1
        }
          
      }
      
      refs <- data.frame(number = seq_along(refs.used), # this is the initial row number (before sort). lib info matches this. 
                         id = refs.used %>% as.numeric, # this is the ref number in the full library (also match.info$ref)
                         name = compound.names.refmat) # name
        refs <- refs[ref.order, ]
        refs$row.mat <- 1:nrow(refs)         # this is the row number in mat
        
      if (is.null(colnames(scores.matrix))){
        colnames(scores.matrix) <- samples.used %>% as.numeric
      }
      
      samples <- data.frame(number = seq_along(samples.used),              # number is column number upon sort
                            id = samples.used %>% as.numeric,              # id is the sample number upon import
                            name = colnames(scores.matrix)) # just use column number for now                            
                            # name = 1:ncol(scores.matrix) %>% as.character) # just use column number for now
        samples <- samples[sample.order, ]
        samples$col.mat <- 1:nrow(samples)                    

      scores.matrix<- scores.matrix[ref.order, ,drop =F]
      scores.matrix<- scores.matrix[, sample.order, drop =F]
      mat <- scores.matrix
      
      
# Layout ##################
ui <-
  fluidPage(
    
    titlePanel(h1(paste0(study))),
    
    # Left side stuff ####
    column(6,
           
           # Reference Spectrum Pane ####
             
            fluidRow(
              
             h3("Reference spectrum: "),
             
             plotlyOutput("ref.plot")
             
            ),
    
           # Stackplot pane ####
           
            fluidRow(
              h3(textOutput("stackplotTitle")),
              
              radioButtons("plot_switch", "Plot type:",
               c("Features on Reference" = "features",
                 "Ref-Features on Samples" = "stackplot"), inline = TRUE),
              
              # verbatimTextOutput("state.description"),
              plotOutput("stack.ref.feats"),
              
              sliderInput(inputId = 'vshift.slide', label = "vshift", min = 0, max = 5, value = 1, step = .05),
              # sliderInput(inputId = 'hshift.slide', "hshift", -.1, .1, 0.5, step = 0.001),
            )
    ),

    # Right side stuff ####
    column(6,
          # Heatmap controls ####
            fluidRow(
              
              selectInput("outputType", "Sort by: ", choices = sort.choices),
              
              conditionalPanel(condition = "input.outputType=='Search for compound...'",
                               selectizeInput(inputId="compound.selection",
                                              label=NULL,multiple=T,choices=NULL))
              
            ),
          
          # Heatmap pane ####
           
            fluidRow(
              
             h3('Evidence Score Matrix'),
             plotlyOutput("heat")
             
            ),

          # Sample scores and selection pane ####
          
            fluidRow(
              
              h3(textOutput("selectedRow_name")),
              
              plotlyOutput("scatterScores")
              
            )
    )
  )


# ------------------------------------------------

########################################################################################

# Server ####
server <- function(input, output, session) {
  
  ####### Initialize values ####
    values <- reactiveValues(plotType = 'features',
                             selectedRow = NULL,
                             selectedCols = NULL,
                             selectedRange = NULL,
                             refplot.xlim = NULL,
                             refplot.xlim.previous = NULL,
                             metab.evidence = NULL,
                             bfs = NULL)

  ####### Right-side stuff ######
  
        # Make heatmap ####
          output$heat <- renderPlotly({
            
            # Plot_ly add_heatmap (but how to select???) ####
              drawHeatmap(mat, source.name = "heatmap", dropRowNames = F) %>% 
                event_register("plotly_click")
            
          })
          
        # Select row from heatmap ####
        
          sel.row.from.heatmap <- reactive({
            
            clickData <- event_data("plotly_click", source = "heatmap")
            
            if (length(clickData) == 0) {
              
              return(NULL)
              
            } else {
              
              
              clickData <- clickData %>% as.list
              message('Row ', clickData$pointNumber[[1]][1] + 1, ' selected from heatmap click')
              clickData$pointNumber[[1]][1] + 1

            }
            
          })
        
        # Update values$selectedRow if sel.row.from.heatmap was updated  ####
          observeEvent(sel.row.from.heatmap(), {
                         values$selectedRow <- sel.row.from.heatmap()
                         
                        # Reset plotly data for our feeder plots
                          values$selectedCols <- NULL
                          values$selectedRange <- NULL

                       })
          
        # Format compound (row) name for printing as title  ####
          output$selectedRow_name <- renderPrint({
            
            req(values$selectedRow)
            if (is.null( values$selectedRow )) {
              
              "Click on a cell in the heatmap to look an individual row (compound)"
              
            } else {
              
              paste0(values$selectedRow %>% refs$name[.],
                     ' match scores in each sample')
              
            }
            
          })
          
        # Make scatterplot of samples x score for the selected reference  ####
          output$scatterScores <- renderPlotly({ ####

            # Check that a row was selected

              req(values$selectedRow)

            # Subset the matrix

              mat.sel <- mat[values$selectedRow, , drop = F]

            # Return the scatterplot of scores x samples

              drawScatterScores(mat.sel, dropRowNames = T, source.name = 'scatter') %>%
                event_register('plotly_selected')

          })
          
        # Gather scatterplot selection if changed  ####
          selectedCols <- reactive({
                
              # Get the column (sample) selection from the scatterplot
                select.box <- event_data("plotly_selected", 
                                         source = "scatter",
                                         priority = "event") # this will reset to null with new ref
                selectedCols <- select.box$x %>% unique

                # Which samples? - is there a sample selection? ####
    
                  if (is.null(selectedCols)) {return(NULL)}
    
              # Report that columns were selected
                message('Selected ', length(selectedCols), ' samples...')

              # Return the column inds
                selectedCols

          })
          
          # Update in reactive variable as well
            observeEvent(selectedCols(),{
              values$selectedCols <- selectedCols()
            })
          
  ####### Left-side stuff  ######

          # Plot type selection ####
        
          update_plotType <- reactive({ input$plot_switch })
          
        # Update values$plotType if needed  ####
          observeEvent(update_plotType(), {
                         values$plotType <- update_plotType()

                       })

        # Plot selected ref, record region selection in event_register ####

          # Make the actual plot
          output$ref.plot = renderPlotly({
            
              req(values$selectedRow)
            
              if (is.null(values$selectedRow)) {message('null row selection'); return(NULL)}

                # The lib.data.processed is filtered, but not ordered. Use inds.
                
                ref <- lib.data.processed[[ values$selectedRow %>% refs$number[.] ]]
                dat <- ref$mapped$data.compressed %>% expand_stacklist
                ref$mapped$data <- dat$data
                ref$mapped$ppm <- ppm
                
              # Plot all the features with their cluster assignments

                plot_spec(ref$mapped$data, ref$mapped$ppm, aucs = NULL, title = ref$compound.name, source.name = 'refspec') %>%
                  layout(dragmode = "zoom") %>%
                  event_register(event = "plotly_brushed")

            })
          
          # Change the title if the selectizer was used
          output$ref.plot.title <- renderPrint({

            req(values$selectized)

            if (is.null(values$selectized)){return(NULL)}

            paste0("Reference Spectrum for ", values$selectized)

          })
          
        # Get the ref range selection from the ref plot, if updated
          selectedRange <- reactive({

            # Identify selected ppm range within ref - is there a selected range in the ref plot?? ####

              sel.range <- event_data(source = 'refspec', 
                                      event = "plotly_brushed", 
                                      priority = "event") # this will reset to null with new ref

                if (is.null(sel.range)){return(NULL)} # and null will be returned

                selectedRange <- sel.range$x[1:2] # ppm units

                  message('selected ref region : ', round(selectedRange[1], 3), '-', round(selectedRange[2], 3), ' ppm')

                selectedRange
          })
          
          # Update in reactive variable as well
            observeEvent(selectedRange(),{
              values$selectedRange <- selectedRange()
            })

        # Update slider.vshift for stackplot

          slider.vshift <- reactive({
                vshift <- input$vshift.slide
                
          })

          # Update in reactive variable as well
          
            observeEvent(slider.vshift(),{
              values$slider.vshift <- slider.vshift()
            })
            
        # Update xlim for stackplot if ref plot zoom changes

           refplot.xlim <- reactive({

            # Identify selected ppm range within ref - is there a selected range in the ref plot?? ####
              # if plotly ref spec zoom is updated, do this:
             
              zoom <- event_data("plotly_relayout", "refspec")
                                       # priority = "event") # this will reset to null with new ref

                #   Unfortunately, the plotly relayout output changes depending on the situation, 
                #   so go through some checks. It will only give a useful region if the user zooms 
                #   in. It will reset to NULL when selection tool is activated. 
                # 
                #   see: https://stackoverflow.com/questions/49268902/get-axis-ranges-after-plotly-resize-in-shiny
                #   (borrowed some code)
              
                if(is.null(zoom) || names(zoom[1]) %in% c("xaxis.autorange", "width")) {
                  # For cases when a new ref plot is made
                  
                  xlim <- NULL

                } else {

                    xlim <- c(zoom$`xaxis.range[0]`, zoom$`xaxis.range[1]`)
                    
                    if (length(xlim) == 0){
                      
                      # For cases when plotly errantly reports null zoom range because selection was activated
                      
                      message('refplot zoom data length(xlim) == 0, but plotly_relayout != NULL:')
                      message('\tselection box activation detected. Using previous zoom region.')
                      
                      # The observeEvent(refplot.xlim()) will parse this:
                        
                        xlim <- "use.previous"
                      
                    } else {
                      # When it works as expected and zoom range is reported. 
                      # the value was set above; just report the new range to message.
                      # NULL values or closures should never make it to this code. 
                      
                      message('zoomed to ref region : ')
                      message(round(xlim[1], 3), '-', round(xlim[2], 3), ' ppm') #
                      
                    }
                }

                xlim
          })

          # Update in reactive variable as well
            # Make decision about whether or not to actually update the range, depending on whether
            # it looks like plotly errantly reset the range. If so, use the previous value (stored
            # as a temp value). 
            
            observeEvent(refplot.xlim(),{
              
              if (is.character(refplot.xlim())){
                
                # If it's a character, that's because it's 'use.previous'. 
                
                  values$refplot.xlim <- values$refplot.xlim.previous
                
              } else {
                
                # Set the value if it looks like a real range
                
                  values$refplot.xlim <- refplot.xlim()
                  
                # Also record as 'previous' value. In the event that values$refplot.xlim
                # is set to NULL because selection tool was used, we can still access it 
                # there. 
                  
                  values$refplot.xlim.previous <- values$refplot.xlim
                
              }
            })
            
        # Print a description of the selections ####
            output$state.description <- renderText({
                  if (is.null(values$selectedRange))
                    {
                      sel.rng <- values$selectedRange

                  } else {
                      sel.rng <- round(values$selectedRange, 3)
                      sel.rng <- paste0(sel.rng[1], '-', sel.rng[2], ' ppm')
                    }

                  paste0('\n\tselected row (compound): ', values$selectedRow %>% refs$name[.],
                         '\n\tselected cols: ', length(values$selectedCols),
                         '\n\tselectedRange: ', sel.rng)#,
                         # '\n\tref plot zoom: ', paste(values$refplot.xlim))
            })
            
          
          
        # Feature Plot: plot a stack plot for each group of overlapping ref feats mapped to bounds ####
          #  Make title ####
          #  Use the stackplot title to communicate instructions on selecting data. 
          
            output$stackplotTitle <- renderText({
              
              # "PCRS features fit to individual spectra"
              
              # Can't do anything without a row selection:

                if (is.null( values$selectedRow )) {

                  return("Click on a cell in the heatmap to look an individual row (compound)")

                }

              # Next, select the range:

                if (is.null( values$selectedRange )) {

                  return("Use the drag box selection tool to select a feature in the reference spectrum above")

                }

              # Finally, if no samples are selected, indicate that's necessary

                if (is.null( values$selectedCols ) ) {

                  return("Select samples (not too many!) in the scores plot pane (lower right)")

                }

                  return(paste0("PCRS features fit to individual spectra for selected region: ",
                                paste0(round(values$selectedRange[1], 3), 
                                       '-', 
                                       round(values$selectedRange[2], 3)), ' ppm'))
    
                
            })
            
          # Feat plot stuff ####
            
            output$stack.ref.feats <- renderPlot({
              
              # Check we have the necessary filters ####
                  
                  req(values$selectedRow)
                  req(values$selectedCols)
                  req(values$selectedRange)
                  req(values$slider.vshift)
                  req(values$refplot.xlim)
                    # No need to filter on this - plot is triggered by selection box (which is more specific)
                    if(is.null(values$refplot.xlim)){return(NULL)}
                  req(ppm.tolerance)
                  # req(vshift)
                  # req(hshift)
                  
                  message('\n\nTrying stackplot for')
                  
              # Which compound? ####
              
                  if (is.null(values$selectedRow)) {return(NULL)}
  
                  message('\t', values$selectedRow %>% refs$name[.])
  
              # Which samples? - is there a sample selection? ####
  
                  if (is.null(values$selectedCols)) {return(NULL)}
  
                  message('\tin ', length(values$selectedCols), ' samples...')
                  
              # If all those pieces are in place, go on to selecting evidence: ####
              
                      # selectedRow <- which(refs$name %in% 'R-Lactate')
                      # selectedCols <- 1:4
                      # if (all(values$selectedCols %in% 1:4) & values$selectedRow == 96){browser()}
                      metab.evidence <- select_evidence_refmet(ref = values$selectedRow %>% refs[.,], 
                                                               sample = values$selectedCols %>% samples[.,],
                                                               # Big objects to subset using ref.ind:
                                                                 features.c,
                                                                 match.info,
                                                                 backfits,
                                                                 rfs.used, # new object with inds of rfs which contributed to scores
                                                                 lib.data.processed,
                                                                 score.name = 'fsaxrval',
                                                               # # Spectral data:
                                                               #   xmat, 
                                                                 ppm)
  
                      rf.fits <- metab.evidence$rf.fits
                      
                      # Does any evidence exist for that ref in those samples?? ####
                        
                        if(is.null(rf.fits)){message('\tno evidence for this region'); return(NULL)
                        } #else {message('\tselected ', length(rf.fits$fit.xrow), ' spec-features.')}

                        # Which ref feats ranges intersect selected ppm rng within tolerance? ####
  
                            # message('\t', c(rf.fits$fit.rngs) %>% range %>% round(3) %>% paste(., collapse = " - "), ' ppm')
                            out.of.range <- apply(rf.fits$fit.rngs, 2, function(rf.rng)
                              {range_intersect(rf.rng, values$selectedRange) %>% is.na %>% all}) %>% unlist
  
                            in.range <- !out.of.range
  
                        # Plot the backfit stackplot ref features that matched to those ppm points ####
  
                            if (any(in.range)){
  
                                message('Plotting ', sum(in.range), ' ref-feature x spectrum associations\n',
                                        '\tfor region ', 
                                        round(values$selectedRange[1], 3), '-', round(values$selectedRange[2], 3), ' ppm\n',
                                        '\tin ', length(values$selectedCols), ' samples...\n',
                                        '\tPlease wait (plotting more spectra and/or wider regions takes more time)...\n')
                                 
                                
                                bfs <- list(fit.feats = rf.fits$fit.feats[in.range, , drop = F],
                                            fit.positions = rf.fits$fit.positions[in.range, , drop = F],
                                            fit.xrow = rf.fits$fit.xrow[in.range],
                                            pass.fit = T)
  
                                # Compute the feature stackplots ####
  
                                  plt.pars <- list(vshift = values$slider.vshift, 
                                                   pixels = c(512, 512), # inc.res
                                                   pointsize = 0, 
                                                   interpolate = T, 
                                                   # exp.by = 0.05,
                                                   xlim = values$refplot.xlim # this should never arrive here as NULL - screen out beforehand
                                                   )
                                  
                                # Decide which plot to make
                                
                                  if (values$plotType == "features"){
                                    
                                    # EST-style plot showing where features annotate the ref
                                      
                                      feature_est_plot(reg = values$refplot.xlim, 
                                                       metab.evidence, 
                                                       features.c,
                                                       ppm,
                                                       plt.pars)
                                    
                                  } else {
                                    
                                    if (values$plotType == "stackplot"){
                                      # browser()
                                      return(
                                        # png(filename= '/Users/mjudge/Desktop/out.png')
                                          fastStack.withFeatures(xmat, ppm, raster = T, bfs = bfs, plt.pars, 
                                                                 res.ratio = .1)
                                        # dev.off()
                                      )

                                    } else {
                                      
                                      return(NULL)
                                      
                                    }
                                  
                                  }
                            }
                        message("\tnothing in range (+\\- ppm match tolerance; ", ppm.tolerance,')...'); return(NULL)
  
            })


  ####### Server-side compound selectizing ####
        
          # Update selectizer ####
          
            observeEvent(input$outputType, # look in outputType (from selectInput)
              {
                if(req(input$outputType == "Search for compound...")){

                  updateSelectizeInput(session,"compound.selection",
                                       "Type/Select metabolite name (first one will be used)",
                                       choices=unique(refs$name), server=T,
                                       options = list(maxOptions = 15,
                                                      maxItems = 1,
                                                      placeholder = 'Start typing...'))

                  # Set working values based on this selection 
                  
                    output$selectized.compound <- renderText({paste0('selectized: ', input$compound.selection)})

                }
              }
            )
          
          # Get value from selectizer choice, map to row index ####
            observeEvent(input$compound.selection, {
              message("selectizer: ", input$compound.selection)

                selectedRows <- which(refs$name == input$compound.selection) # in terms of matrix row inds

                if (length(selectedRows) == 0){
                  message('compound name mismatch?')
                  values$selectedRow <- NULL
                }

                if (length(selectedRows) > 1){

                  message('Multiple compounds have name: ', refs$name[selectedRows[1]])
                  message('\tchoosing one with best overall score...')
                  row.bestscore <- mat[selectedRows, , drop = F] %>% rowSums %>% which.max %>% selectedRows[.]
                  message('Row ', row.bestscore, ' selected')
                  values$selectedRow <- row.bestscore

                } else {

                    message('Row ', selectedRows[1], ' selected')
                    values$selectedRow <- selectedRows[1]
                  
                }
                
              # Reset plotly data for our feeder plots
                values$selectedCols <- NULL
                values$selectedRange <- NULL

            })

}

# Run app ####



shinyApp(ui, server)
  
  
}

  
# Troubleshooting evidence selection ####
# selectedRow <- which(refs$name %in% 'Citrate')
# selectedCols <- 1:4
# metab.evidence <- select_evidence_refmet(ref = selectedRow %>% refs[.,], 
#                                          sample = selectedCols %>% samples[.,],
#                                          # Big objects to subset using ref.ind:
#                                            features.c,
#                                            match.info,
#                                            backfits,
#                                            rfs.used, # new object with inds of rfs which contributed to scores
#                                            lib.data.processed,
#                                          # # Spectral data:
#                                          #   xmat, ppm,
#                                          # Filtering thresholds:
#                                            ppm.tolerance = ppm.tolerance,
#                                            cutoff.residuals.feat = cutoff.residuals.feat,
#                                            cutoff.residuals.spec = cutoff.residuals.spec)
      
      
      
      
      
      



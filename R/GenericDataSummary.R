#' Generate an html R Markdown document briefly summarizing a dataset
#' 
#' \code{GenericDataSummary} creates a Markdown document for any dataset that 
#' separately summarizes continuous variables in that dataset with a summary 
#' table and boxplots and provides frequency counts for all categorical data. 
#' There is also an option to include a data dictionary table in this simple 
#' summary document. Uses for this function include quickly inspecting a 
#' dataset for aberrant values and getting a sense for which variables are 
#' included.
#' 
#' @param filename A file name and path for the .html R Markdown document to be 
#' generated. Defaults to "DataSummary_X.html" where X is the system date and 
#' time when the document is saved. You can specify the ".html" suffix here or
#' omit it, in which case it will be added for you.
#' @param directory The folder in which you'd like to save the html document. 
#' This defaults to your current working directory.
#' @param Data A data frame to be summarized.
#' @param Dict An optional data frame for the data dictionary to be included for 
#' reference as a table in the Markdown document. This defaults to NULL.
#' @param idvar An optional single character value indicating the column name of 
#' a column of IDs in Data. If specified, a short statement on the number of 
#' unique IDs and counts of rows per ID will be added to the introduction of 
#' the document. This defaults to NULL.
#' 
#' @details
#' This Markdown document makes use of \link[DT]{datatable} from the DT package 
#' for the data dictionary and continuous summary statistics table, which can 
#' be sorted and searched. Categorical variable counts are summarized in a 
#' datatable with a column showing a bar plot of percentages.
#' 
#' @return \code{GenericDataSummary} saves an html document to the specified 
#' filename and directory.
#' 
#' @import DT ggplot2 cowplot htmltools knitr rmarkdown dplyr
#' @export GenericDataSummary
#' 
#' @examples
#' # Make a data dictionary for the built-in iris dataset
#' irisdd <- data.frame(Variable = names(iris), 
#' Description = gsub("\\.", " ", names(iris)), Levels = "")
#' irisdd[irisdd$Variable == "Species", "Levels"] <- 
#' paste(unique(iris$Species), collapse = ", ")
#' 
#' # create the summary document
#' GenericDataSummary("NewIrisSummary", getwd(), iris, irisdd)
#' 

GenericDataSummary <- function(filename = NULL, directory = getwd(), 
                               Data, Dict = NULL, idvar = NULL){
  #QC filename and Data/Dict
  if(is.null(filename)){
    filename <- paste0("DataSummary_", format(Sys.time(), "%Y%m%d_%H%M"), ".html")
  }
  if(!grepl("\\.html$", filename)) filename <- paste0(filename, ".html")
  if(!is.data.frame(Data)) stop("Data must be an object of class data.frame")
  if(!is.null(Dict)){
    if(!is.data.frame(Dict)) stop("Dict must be an object of class data.frame")
  }
  if(!is.null(idvar)) if(!is.character(idvar) | length(idvar) > 1){
    stop("idvar must be a single character value indicating an ID column name.")
  }
  
  #grab file path for Rmd
  rmdpath <- system.file("rmd", "GenericDataSummary.Rmd", 
                         package = "DrewDayRFunctions")
  
  rmarkdown::render(
    input = rmdpath,
    output_dir = directory,
    output_file = filename)
}

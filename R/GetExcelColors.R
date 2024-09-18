#' \code{GetExcelColors} is a function for reading an xlsx file and generating
#' a data frame with all the column names preserved but the contents of the 
#' columns replaced with hexidecimal color codes for the background colors in
#' the Excel files. The function utilizes \code{\link[readxl]{read_xlsx}} to 
#' read in the Excel file and passes arguments to that function. Then it uses
#' \code{\link[tidyxl]{xlsx_formats}} to get the formatting for the cells.
#' 
#' @param filepath This is a "X.xlsx" file path pointing to an Excel .xlsx file.
#' @param sheet This is either the integer index or name of an Excel sheet to 
#' process. This defaults to 1, meaning the function will process the first 
#' sheet/tab in the .xlsx file.
#' @param skip This defines the number of rows to skip when reading the .xlsx 
#' sheet. This is useful if there are extra rows before the column header row. 
#' This defaults to 0, meaning that no rows are skipped.
#' @param skiprows This defines additional rows to be omitted from the data frame
#' of background colors. For example, there may be subtitle rows below the column
#' header row that one would want to omit while keeping that row of column names
#' above it. This defaults to NULL, meaning no rows are skipped beyond what is
#' defined in the \code{skip} argument.
#' 
#' @return \code{GetExcelColors} returns a list of the following:
#' \item{ColorData}{A data frame containing column names
#' and background Excel colors. NA values indicate that no background was present
#' in those cells. This object is of class "data.frame".}
#' \item{UniqueColors}{This is a vector of the unique hexadecimal colors in 
#' the dataset.}
#' 
#' @import tidyxl readxl
#' @export GetExcelColors
#' 
#' @examples 
#' 
#' # The following is just an example, but won't actually run since the filepath
#' #  doesn't lead to any actual .xlsx file.
#' 
#' \dontrun{
#' dir <- "X:/"
#' xlfile <- "ExampleData.xlsx"
#' ExampleData_Colors <- GetExcelColors(paste(dir, xlfile, sep = "/"), sheet = 2,
#' skip = 4, skiprows = 6:7)
#' # this grabs color data from the 2nd sheet, skipping the first 4 rows as well
#' #  as rows 6 to 7.
#' }


GetExcelColors <- function(filepath, sheet = 1, skip = 0, skiprows = NULL){
  Data <- readxl::read_xlsx(filepath, sheet = sheet, skip = skip)
  fill_colors <- 
    tidyxl::xlsx_formats(filepath)$local$fill$patternFill$fgColor$rgb
  fills <- 
    xlsx_cells(filepath, sheets = sheet) %>%
    mutate(fill_color = fill_colors[local_format_id]) %>%
    dplyr::select(row, col, fill_color) %>%
    spread(col, fill_color) 
  if(skip != 0) fills <- fills[which(fills$row > skip), ]
  if(!is.null(skiprows)) fills <- fills[which(!fills$row %in% (skiprows - 1)), ]
  fills <- fills %>% select(-row)
  names(fills) <- names(Data)
  return(list(ColorData = as.data.frame(fills), UniqueColors = unique(fill_colors)))
}
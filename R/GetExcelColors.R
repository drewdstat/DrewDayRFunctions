#' \code{GetExcelColors} is a function for reading an xlsx file and generating
#' a data frame with all the column names preserved but the contents of the 
#' columns replaced with Excel color codes for the background colors in
#' the Excel files. Additional output includes a data frame and heatmap plot 
#' showing the counts of cells with each color. The function utilizes 
#' \code{\link[readxl]{read_xlsx}} to read in the Excel file, and then it uses
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
#' and background Excel colors. NA values indicate that no background was 
#' present in those cells. This object is of class "data.frame".}
#' \item{CellCounts}{This is a data frame listing out the unique colors detected
#' in the spreadsheet and how many cells have those colors, organized in 
#' descending order by count.}
#' \item{CountPlot}{This plots the counts of cells with each color in a heatmap 
#' that is colored by the Excel colors present in the spreadsheet. This figure 
#' is a useful way of seeing what the colors look like and how frequent they 
#' are.}
#' 
#' @details
#' Note that the color codes exported by Excel do not directly correspond to the 
#' hexademical color codes utilized by html, R, and other languages. To convert 
#' the 8-character Excel color codes to hexadecimal color codes, simply replace 
#' the leading 2 "FF" characters with "#". 
#' 
#' @import ggplot2 tidyxl readxl
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
  clrsumm <- summary(as.factor(unlist(fills, use.names = F)))
  if("NA's" %in% names(clrsumm)) clrsumm <- 
    clrsumm[-grep("NA's", names(clrsumm), fixed = T)]
  clrsumm <- data.frame(Color = names(clrsumm), Count = clrsumm)
  clrsumm <- clrsumm[with(clrsumm, order(-Count)), ]
  row.names(clrsumm) <- NULL
  plotdat <- clrsumm
  hexclr <- paste0("#", substr(plotdat$Color, 3, 8))
  names(hexclr) <- plotdat$Color
  plotdat$Color <- factor(plotdat$Color, levels = rev(plotdat$Color))
  gg1 <- ggplot(plotdat, aes(x = 1, y = Color)) + theme_minimal() + 
    geom_tile(aes(fill = Color)) + 
    geom_text(aes(label = Count)) + 
    scale_fill_manual(values = hexclr) + 
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(), legend.position = "none")
  return(list(ColorData = as.data.frame(fills), CellCounts = clrsumm, 
              CountPlot = gg1))
}
# The following function is from code here: 
# https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames

read_excel_allsheets <- function(filename, tibble = FALSE) {
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(sheet_i) readxl::read_excel(filename,
                                                             sheet = sheet_i))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    return(x)
}

find_depth <- function(x, threshold = 10000) {
    n <- length(x)
    for (depth in 1:(n - 1)) {
        for (i in 1:(n - depth)) {
            if (abs(x[i + depth] - x[i]) <= threshold) {
                break
            }
            if (i == n - depth) {
                return(depth)
            }
        }
    }
    return(n)  # fallback if all pairs are within threshold
}

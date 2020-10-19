sourceDir <- function(path = "C:/RChipVersion3.0/RChip/R", trace = TRUE, ...)
{
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if (trace) cat(nm,":")
       source(file.path(path, nm), ...)
       if (trace) cat("\n")
    }
}

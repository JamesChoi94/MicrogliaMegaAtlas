
# Generate summary statistics

library(Seurat)
library(ggplot2)
library(dplyr)

seurat_data <- 'data/mg.rds'
stats_dir <- 'SummaryStats/'
summarise_these <- c('study','tissue','injury_model','sex','age')

mg <- readRDS(file = seurat_data)
dir.create(path = stats_dir)

summary_stats <- lapply(
  X = summarise_these,
  FUN = function(x) {
    return(sort(table(mg@meta.data[x]), decreasing = TRUE))
  }
)
names(summary_stats) <- summarise_these

parse_csv_to_md_table <- function(x) {
  # input is a named vector where names are first row of table and elements are
  # second row of table.
  header <- paste(names(x), collapse = ' | ')
  spacer <- paste(rep('---', length(x)), collapse = ' | ')
  contents <- paste(x, collapse = ' | ')
  return(c(header, spacer, contents, ''))
}


add_summary_stats_to_readme <- function(x, readme_path) {
  outs <- lapply(
    X = x,
    FUN = parse_csv_to_md_table
  )
  names(outs) <- names(x)
  readme_text <- readLines(con = readme_path)
  print(readme_text)
  stats_start <- grep(pattern = 'Summary Statistics', x = readme_text)
  if (length(stats_start) == 0) {
    if (nchar(readme_text[length(readme_text)]) != 0) {
      write(file = readme_path, x = '', append = TRUE, sep = '\n')
    }
    write(file = readme_path, 
          x = c('## Summary Statistics', ''),
          append = TRUE)
    for (i in 1:length(outs)) {
      write(x = c(tools::toTitleCase(names(outs)[i]), unlist(outs[[i]])),
            file = readme_path,
            sep = '\n',
            append = TRUE)
    }
  } else {
    print(c(readme_text[1:stats_start[1]], ''))
    write(file = readme_path,
          x = c(readme_text[1:stats_start[1]], ''),
          sep = '\n',
          append = FALSE)
    for (i in 1:length(outs)) {
      write(x = c(tools::toTitleCase(names(outs)[i]), unlist(outs[[i]])),
            file = readme_path,
            sep = '\n',
            append = TRUE)
    }
  }
}


add_summary_stats_to_readme(x = summary_stats, readme_path = './README.md')


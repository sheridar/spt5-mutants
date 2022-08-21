library(rmarkdown)
library(docopt)

doc <- "Usage: knit_Rmd.R [--help] [--input INPUT] [--output OUT]

-i --input INPUT    path to rmarkdown
-o --output OUTPUT  name of output html file to write
-h --help           display this help message"

opts <- docopt(doc)

print(opts)


# Render Rmd
render(
  input       = opts$input,
  output_file = opts$output
)



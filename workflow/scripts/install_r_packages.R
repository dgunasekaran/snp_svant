options(repos = c(CRAN = "https://cloud.r-project.org"))

# Source: https://github.com/trevorld/r-argparse
if (!require("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

# Source: https://ggplot2.tidyverse.org/
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Source: https://bioconductor.org/packages/release/bioc/html/apeglm.html
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Source: https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html
if (!require("VariantAnnotation", quietly = TRUE)) {
  BiocManager::install("VariantAnnotation")
}

# Source: https://bioconductor.org/packages/release/bioc/html/StructuralVariantAnnotation.html
if (!require("StructuralVariantAnnotation", quietly = TRUE)) {
  BiocManager::install("StructuralVariantAnnotation")
}

# Source: https://www.tidyverse.org/
if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

# Source: https://cran.r-project.org/web/packages/stringr/readme/README.html
if (!require("stringr", quietly = TRUE)) {
  install.packages("stringr")
}




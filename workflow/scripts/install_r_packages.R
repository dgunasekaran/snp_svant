BiocManager::install(version = '3.17', ask=FALSE)

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Source: https://github.com/trevorld/r-argparse
if (!require("argparse", quietly = TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/argparse/argparse_2.2.2.tar.gz", repos=NULL, type="source")
}

# Source: https://ggplot2.tidyverse.org/
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!require("lattice", quietly = TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.21-9.tar.gz", repos=NULL, type="source")
}

if (!require("Matrix", quietly = TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.1.tar.gz", repos=NULL, type="source")
}

# Installing older version of BiocManager
install.packages("https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.22.tar.gz", repos=NULL, type="source")
BiocManager::install(version = '3.17', ask=FALSE)

# Source: https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html
if (!require("VariantAnnotation", quietly = TRUE)) {
  BiocManager::install("VariantAnnotation", version="3.17")
}

# Source: https://bioconductor.org/packages/release/bioc/html/StructuralVariantAnnotation.html
if (!require("StructuralVariantAnnotation", quietly = TRUE)) {
  BiocManager::install("StructuralVariantAnnotation", version="3.17")
}

# Source: https://www.tidyverse.org/
if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

# Source: https://cran.r-project.org/web/packages/stringr/readme/README.html
if (!require("stringr", quietly = TRUE)) {
  install.packages("stringr")
}




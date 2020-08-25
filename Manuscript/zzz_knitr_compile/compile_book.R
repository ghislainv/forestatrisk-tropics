#!/usr/bin/Rscript

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :CC-BY-SA 4.0
# ==============================================================================

# Libraries
library(bookdown)
library(knitr)
library(here)

# ==============================================================================
# Article
# ==============================================================================

# Working directory
setwd(here::here("Manuscript/Article"))

# Bookdown
# pdf
options(knitr.table.format="latex")
pdf_format <- bookdown::pdf_document2(citation_package="natbib", fig_caption=TRUE, keep_tex=FALSE,
																			latex_engine="pdflatex", number_sections=TRUE, toc=FALSE,
																			includes=list(in_header="header.tex", before_body="doc_prefix.tex"))
params <- list(title="",author="",date="")
bookdown::render_book("index.Rmd", output_format=pdf_format)

# params
title_html <- 'Forecasting tropical moist forest cover change in the 21st century under a "business-as-usual" scenario'
author_html <- "Ghislain VIEILLEDENT, Christelle VANCUTSEM, and Frédéric ACHARD"
date_html <- format(Sys.time(), "%d %B, %Y")
params <- list(title=title_html,author=author_html, date=date_html)

# docx
# Don't indicate output_format to take into account YAML options
options(knitr.table.format="html")
# Dynamic YAML options
docx_format <- bookdown::word_document2(fig_caption=TRUE, toc=FALSE)
bookdown::render_book("index.Rmd", output_format=docx_format)

# # html
# # Don't indicate output_format to take into account YAML options
# options(knitr.table.format="html")
# # Dynamic YAML options
# bookdown::render_book("index.Rmd")

# ==============================================================================
# Supplementary materials
# ==============================================================================

# Libraries
require(bookdown)
require(knitr)
require(here)

# Working directory
setwd(here::here("Manuscript/Supplementary_Materials"))

# Bookdown
# pdf
options(knitr.table.format="latex")
pdf_format <- bookdown::pdf_document2(citation_package="natbib", fig_caption=TRUE, keep_tex=FALSE,
                                      latex_engine="pdflatex", number_sections=TRUE, toc=FALSE,
                                      includes=list(in_header="header.tex", before_body="doc_prefix.tex"))
params <- list(title="",author="",date="")
bookdown::render_book("index.Rmd", output_format=pdf_format)

# params
title_html <- 'SUPPLEMENTARY MATERIALS\nForecasting tropical moist forest cover change in the 21st century under a "business-as-usual" scenario'
author_html <- "Ghislain VIEILLEDENT, Christelle VANCUTSEM, and Frédéric ACHARD"
date_html <- format(Sys.time(), "%d %B, %Y")
params <- list(title=title_html,author=author_html, date=date_html)

# docx
# Don't indicate output_format to take into account YAML options
options(knitr.table.format="html")
# Dynamic YAML options
docx_format <- bookdown::word_document2(fig_caption=TRUE, toc=FALSE)
bookdown::render_book("index.Rmd", output_format=docx_format)

# # html
# # Don't indicate output_format to take into account YAML options
# options(knitr.table.format="html")
# # Dynamic YAML options
# bookdown::render_book("index.Rmd")

# EOF

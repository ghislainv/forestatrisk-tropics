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
pdf_format <- bookdown::pdf_document2(citation_package="natbib", fig_caption=TRUE, keep_tex=TRUE, keep_md=FALSE,
																			latex_engine="pdflatex", number_sections=FALSE, toc=FALSE,
																			includes=list(in_header="header.tex", before_body="doc_prefix.tex"))
params <- list(title="",author="",date="")
full_width_type=FALSE
font_size_type=10
bookdown::render_book("index.Rmd", output_format=pdf_format)

# params
title_html <- "Spatial forecasting of forest cover change in the humid tropics over the 21^st^ century"
author_html <- "Ghislain VIEILLEDENT, Christelle VANCUTSEM, and Frédéric ACHARD"
date_html <- format(Sys.time(), "%d %B, %Y")
params <- list(title=title_html,author=author_html, date=date_html)

# html
options(knitr.table.format="html")
full_width_type=TRUE
font_size_type=NULL
# Dynamic YAML options
html_format <- bookdown::html_document2(number_sections=FALSE, fig_caption=TRUE, toc=FALSE)
bookdown::render_book("index.Rmd", output_format=html_format)

# Then the html can be imported in Google Doc

# # book
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
pdf_format <- bookdown::pdf_document2(citation_package="natbib", fig_caption=TRUE, keep_tex=TRUE, keep_md=FALSE,
                                      latex_engine="pdflatex", number_sections=FALSE, toc=FALSE,
                                      includes=list(in_header="header.tex", before_body="doc_prefix.tex"))
params <- list(title="",author="",date="")
bookdown::render_book("index.Rmd", output_format=pdf_format)

# # params
# title_html <- "SUPPLEMENTARY MATERIALS\nSpatial forecasting of forest cover change in the humid tropics over the 21st century"
# author_html <- "Ghislain VIEILLEDENT, Christelle VANCUTSEM, and Frédéric ACHARD"
# date_html <- format(Sys.time(), "%d %B, %Y")
# params <- list(title=title_html,author=author_html, date=date_html)

# # html
# # Don't indicate output_format to take into account YAML options
# options(knitr.table.format="html")
# # Dynamic YAML options
# bookdown::render_book("index.Rmd")

# # ==============================================================================
# # Combined pdf
# # ==============================================================================
# 
# require(glue)
# require(here)
# 
# Manuscript_pdf <- here("Manuscript", "Article", "doc", "article.pdf")
# SM_pdf <- here("Manuscript", "Supplementary_Materials", "doc", "sm.pdf")
# Combined_pdf <- here("Manuscript", "Combined", "combined.pdf")
# cmd <- glue("pdftk {Manuscript_pdf} {SM_pdf} cat output {Combined_pdf}")
# system(cmd)

# ==============================================================================
# Data S
# ==============================================================================

# # Libraries
# require(bookdown)
# require(knitr)
# require(here)
# 
# # Working directory
# setwd(here::here("Manuscript/Supplementary_Data"))
# 
# # params
# title_html <- "Spatial forecasting of forest cover change in the humid tropics over the 21^st^ century"
# subtitle_html <- "Supplementary Data"
# author_html <- "Ghislain VIEILLEDENT, Christelle VANCUTSEM, and Frédéric ACHARD"
# date_html <- ""
# params <- list(title=title_html, subtitle=subtitle_html, author=author_html, date=date_html)

# # html
# options(knitr.table.format="html")
# full_width_type=FALSE
# font_size_type=NULL
# # Dynamic YAML options
# html_format <- bookdown::html_document2(number_sections=FALSE, fig_caption=TRUE, toc=TRUE, toc_float=TRUE)
# bookdown::render_book("index.Rmd", output_format=html_format)

# Move files to server
# system("scp ~/Code/forestatrisk-tropics/Manuscript/Supplementary_Data/supplementary-data.html fdb:/home/www/forestatrisk/tropics/supplementary-data/index.html")
system("scp ~/Code/forestatrisk-tropics/Manuscript/Supplementary_Data/tables_website/* fdb:/home/www/forestatrisk/tropics/supplementary-data/")

# ==============================================================================
# Combined with LaTex
# ==============================================================================

require(here)

# Working directory
setwd(here::here("Manuscript/LaTeX"))

# Copy article
f_in <- here("Manuscript", "Article", "doc", "article.tex")
f_out <- here("Manuscript", "LaTeX", "article.tex")
file.copy(f_in, f_out, overwrite=TRUE)

# Copy sm
f_in <- here("Manuscript", "Supplementary_Materials", "doc", "sm.tex")
f_out <- here("Manuscript", "LaTeX", "sm.tex")
file.copy(f_in, f_out, overwrite=TRUE)

# Combine article and sm
# article
article <- readLines(f_in <- here("Manuscript", "LaTeX", "article.tex"))
article <- gsub(pattern="\\{figures/", replace="\\{figs_article/", x=article)
l_a <- length(article)
article <- article[-l_a]
# sm
sm <- readLines(f_in <- here("Manuscript", "LaTeX", "sm.tex"))
sm <- gsub(pattern="\\{figures/", replace="\\{figs_sm/", x=sm)
l_sm <- length(sm)
sm <- sm[-c(1:110, (l_sm-8):(l_sm-1))]
# sm header
sm_head <- c("\\renewcommand{\\thetable}{S\\arabic{table}}",
             "\\renewcommand{\\thefigure}{S\\arabic{figure}}",
             "\\renewcommand{\\theequation}{S\\arabic{equation}}",
             "\\setcounter{figure}{0}",
             "\\setcounter{table}{0}",
             "\\newpage")
# Combine
combined <- c(article, sm_head, sm)
writeLines(combined, here("Manuscript", "LaTeX", "combined.tex"))

# Compile with pdflatex
tools::texi2dvi(here("Manuscript", "LaTeX", "combined.tex"), pdf=TRUE, clean=TRUE)

# EOF
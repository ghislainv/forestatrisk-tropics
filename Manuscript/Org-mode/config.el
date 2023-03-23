;;; config.el --- compiling the article in pdf -*- lexical-binding: t -*-
;;; Copyright (C) 2021 Ghislain Vieilledent

;;; Commentary:
;;; Use this configuration file through build.sh and Makefile

;;; Code:

;; --------------------------------------
;; INSTALL PACKAGES
;; --------------------------------------

(require 'package)
(add-to-list 'package-archives '("org" . "https://orgmode.org/elpa/") t)
(add-to-list 'package-archives '("nongnu" . "https://elpa.nongnu.org/nongnu/") t)
(package-initialize)

;; Fetch the list of packages available
(unless package-archive-contents
  (package-refresh-contents))

;; Install packages
(defvar my-package-list)
(setq my-package-list '(org
			org-contrib))
(dolist (i-package my-package-list)
  (unless (package-installed-p i-package)
    (package-install i-package)))

;; Load packages
(require 'org)
(require 'ox-latex)
(require 'ox-extra)
(require 'oc-csl)
(require 'oc-natbib)

;; --------------------------------------
;; ORG --> PDF
;; --------------------------------------

;; Export to LaTeX
(unless (boundp 'org-latex-classes)
  (setq org-latex-classes nil))
;; LaTeX classes
(add-to-list 'org-latex-classes
             '("koma-article" "\\documentclass{scrartcl}"
               ("\\section{%s}" . "\\section*{%s}")
               ("\\subsection{%s}" . "\\subsection*{%s}")
               ("\\subsubsection{%s}" . "\\subsubsection*{%s}")
               ("\\paragraph{%s}" . "\\paragraph*{%s}")
               ("\\subparagraph{%s}" . "\\subparagraph*{%s}")))

;; Export process from orgmode to LaTeX to PDF
(setq org-latex-pdf-process '("texi2dvi --pdf --clean --verbose --batch %f"))
;; (setq org-latex-pdf-process
;;       '("pdflatex -shell-escape -interaction nonstopmode -output-directory %o %f"
;;         "pdflatex -shell-escape -interaction nonstopmode -output-directory %o %f"
;;         "pdflatex -shell-escape -interaction nonstopmode -output-directory %o %f"))

;; In org-mode 9 you need to have #+PROPERTY: header-args :eval never-export
;; in the beginning or your document to tell org-mode not to evaluate every
;; code block every time you export.
(setq org-confirm-babel-evaluate nil) ;; Do not ask for confirmation all the time!!
(setq org-src-preserve-indentation t)
(setq org-edit-src-content nil)
(setq org-export-with-smart-quotes t)

;; Languages evaluated
(org-babel-do-load-languages
 'org-babel-load-languages
 '((emacs-lisp . t)
  (shell . t)
  (python . t)
  (R . t)
  (screen . t)
  (org . t)
  (makefile . t)))

;; Add :ignore: tag to ignore headline
(ox-extras-activate '(ignore-headlines))

;; Processor for org-cite
(defvar bibdir (expand-file-name "biblio" default-directory)
  "Biblio directory.")
(setq org-cite-csl-styles-dir bibdir)

;; Function to export to pdf
(defun export-to-pdf nil
  "Export to pdf.
This function exports an org file to a pdf"
  (find-file "paper.org")
  (org-latex-export-to-pdf nil nil nil nil nil))

;; Call to the function
(export-to-pdf)
(message "----OooooOOooooO----")

;;; config.el ends here

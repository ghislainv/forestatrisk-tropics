;;; org2pdf.el --- compiling the article pdf -*- lexical-binding: t -*-
;; Copyright (C) 2021 Ghislain Vieilledent

;;; Commentary:
;; None

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

;; --------------------------------------
;; ORG MODE
;; --------------------------------------

;; Open in org-mode
(require 'org)

;; Export to LaTeX
(require 'ox-latex)
(unless (boundp 'org-latex-classes)
  (setq org-latex-classes nil))
;; LaTeX classes
;; Pinp Is Not PNAS
(add-to-list 'org-latex-classes
	     '("pinp-article" "\\documentclass{pinp}"
               ("\\section{%s}" . "\\section*{%s}")
               ("\\subsection{%s}" . "\\subsection*{%s}")
               ("\\subsubsection{%s}" . "\\subsubsection*{%s}")
               ("\\paragraph{%s}" . "\\paragraph*{%s}")
               ("\\subparagraph{%s}" . "\\subparagraph*{%s}")))
;; Export process from orgmode to LaTeX to PDF
(setq org-latex-pdf-process '("texi2dvi --pdf --clean --verbose --batch %f"))

;; In org-mode 9 you need to have #+PROPERTY: header-args :eval never-export
;; in the beginning or your document to tell org-mode not to evaluate every
;; code block every time you export.
(setq org-confirm-babel-evaluate nil) ;; Do not ask for confirmation all the time!!

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
(require 'ox-extra)
(ox-extras-activate '(ignore-headlines))

;; Processor for org-cite
(require 'oc-csl)
(require 'oc-natbib)

;;; org2pdf.el ends here

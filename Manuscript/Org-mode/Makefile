# Makefile example at https://github.com/schnorr/ieeeorg/blob/master/Makefile

paper:
	emacs -Q --batch -l org2pdf.el paper.org -f org-latex-export-to-pdf
	rm paper.tex~

website: paper.pdf
	cp paper.pdf ../../Website/article/Vieilledent2022-preprint.pdf
	Rscript -e 'rmarkdown::render_site("~/Code/forestatrisk-tropics/Website")'
	rsync -azP ~/Code/forestatrisk-tropics/Website/_site/. fdb:/home/www/forestatrisk/

clean:
	rm paper.tex paper.pdf

# Makefile example at https://github.com/schnorr/ieeeorg/blob/master/Makefile

paper:
	./build.sh

website: paper.pdf
	cp paper.pdf ../../Website/article/Vieilledent2022-preprint.pdf
	Rscript -e 'rmarkdown::render_site("~/Code/forestatrisk-tropics/Website")'
	rsync -azP ~/Code/forestatrisk-tropics/Website/_site/. fdb:/home/www/forestatrisk/

clean:
	rm -f paper.tex paper.tex~

TEX = latexme.tex xgap.tex xgaplatt.tex ../doc.x11/install.tex
AUX = latexme.aux xgap.aux xgaplatt.aux ../doc.x11/install.aux

all:
	@echo
	@echo "The following targets are supported:"
	@echo
	@echo " 'doc'                   latex the file 'latexme.tex'"
	@echo " 'clean'                 remove all created files"
	@echo

doc: latexme.dvi

latexme.dvi: $(TEX)
	latex latexme.tex

clean:
	rm -f *.aux *.log *.dvi *.lj *~

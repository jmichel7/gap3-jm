# Making the .pdf and .html manuals has become something of a black art,
# so I automated it by this makefile. JM January 2015.

man:
	pdflatex manual
	cp ~/gap3-jm/doc/manual.bib .
	bibtex manual
	makeindex manual.idx
	pdflatex manual
	pdflatex manual
	ruby ~/gap3-jm/tools/convert.rb -c .
	rm manual.bib
	mv index.htm index.html
	mv *.htm* ../htm
	mv manual.pdf ..

# erase auxiliary files
clean:
	-rm *.aux *.log *.ind *.ilg *.blg *.bbl *.idx *.toc *.htm err outgap outman script index.html
test:
	ruby ~/gap3-jm/tools/maketest.rb curv

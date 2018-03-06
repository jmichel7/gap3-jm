# Making the .pdf and .html manuals has become something of a black art,
# so I automated it by this makefile. JM January 2015.

man:
	pdflatex manual
	cp ~/gap3-jm/doc/manual.bib .
	bibtex manual
	makeindex manual.idx
	pdflatex manual
	pdflatex manual
	mv manual.pdf algebra.pdf
	ruby ~/gap3-jm/tools/convert.rb -c .
	rm manual.bib
	cp index.htm index.html

# erase auxiliary files except manual.toc manual.idx necessary for GAP help
clean:
	-rm *.aux *.log *.ind *.idx *.ilg *.blg *.bbl *.toc *.htm* err outgap outman script

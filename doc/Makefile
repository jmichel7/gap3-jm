# Making the .pdf and .html manuals has become something of a black art,
# so I automated it by this makefile. JM January 2015.

man:
	pdflatex manual
	bibtex manual
	makeindex manual.idx
	pdflatex manual
	pdflatex manual
	ruby ../tools/convert.rb -c .
	cp index.htm index.html
	mv *.htm* ../htm
	mv manual.pdf ..

# erase auxiliary files except manual.toc manual.idx necessary for GAP help
clean:
	-rm *.aux *.log *.ind *.ilg *.blg *.bbl err outgap outman script

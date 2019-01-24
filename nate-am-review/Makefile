TEX_FILES  =  $(wildcard *.tex) $(wildcard ./Section*/*.tex)
TARGET= Main
BIB=APR

$(TARGET).pdf: $(TEX_FILES) $(TARGET).bbl $(TARGET).aux
	pdflatex -shell-escape -draftmode $(TARGET)
	pdflatex -shell-escape $(TARGET)

$(TARGET).bbl $(TARGET).aux: $(BIB).bib
	# ensure we make a .aux file
	pdflatex -shell-escape -draftmode $(TARGET)
	pdflatex -shell-escape -draftmode $(TARGET)
	# use the .aux file
	bibtex $(TARGET)
	# rerun to fix up the .aux file
	pdflatex -shell-escape -draftmode $(TARGET)
	bibtex $(TARGET)

clean:
	rm -f *~ *.ilg *bak *.idx *.ind *.aux *.toc *.ps *.log *.lof *.loa
	rm -f *.bbl *.blg *.dvi *.out $(TARGET).pdf *.ps  *.los *.lot *.tdo
	rm -f *Notes.bib

.PHONY: clean


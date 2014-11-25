all: he3en.pdf

%.pdf: %.ps
	ps2pdf $<
%.ps: %.dvi
	dvips $<
%.dvi: %.tex
	latex $<

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.dvi *.eps *.bak *.bbl *.blg

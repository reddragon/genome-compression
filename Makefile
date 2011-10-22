proposal:
	pdflatex proposal.tex
	bibtex proposal.tex
	pdflatex proposal.tex

clean:
	rm -f proposal.pdf proposal.aux proposal.log proposal.bbl proposal.toc proposal.blg

SRCP = hdnum-tutorial-beamer
SRCA = hdnum-tutorial-article

beamer: 
	latex $(SRCP)
	bibtex $(SRCP)
	bibtex $(SRCP)
	bibtex $(SRCP)
	dvipdf $(SRCP)

art:
	latex $(SRCA)
	bibtex $(SRCA)
	bibtex $(SRCA)
	bibtex $(SRCA)
	dvipdf $(SRCA)

full:
	latex $(SRCP)
	bibtex $(SRCP)
	bibtex $(SRCP)
	bibtex $(SRCP)
	latex $(SRCP)
	latex $(SRCP)
	dvipdf $(SRCP)
	latex $(SRCA)
	bibtex $(SRCA)
	bibtex $(SRCA)
	bibtex $(SRCA)
	latex $(SRCA)
	latex $(SRCA)
	dvipdf $(SRCA)


picture:
	pdflatex picture
	xpdf picture.pdf

split: 
	pdftk A=numstoch-article.pdf cat A1-15   output pdf/numstoch01.pdf
	pdftk A=numstoch-article.pdf cat A16-25  output pdf/numstoch02.pdf
	pdftk A=numstoch-article.pdf cat A26-33  output pdf/numstoch03.pdf
	pdftk A=numstoch-article.pdf cat A34-48  output pdf/numstoch04.pdf
	pdftk A=numstoch-article.pdf cat A49-58  output pdf/numstoch05.pdf
	pdftk A=numstoch-article.pdf cat A59-73  output pdf/numstoch06.pdf
	pdftk A=numstoch-article.pdf cat A74-86  output pdf/numstoch07.pdf
	pdftk A=numstoch-article.pdf cat A87-99  output pdf/numstoch08.pdf
	pdftk A=numstoch-article.pdf cat A100-101  output pdf/inhalt.pdf

clean:
	rm *.bbl *.blg *.dvi *.aux *.idx *.log *.lof *.toc *.lot *.ind *.ilg *.snm *.nav *.out *.vrb

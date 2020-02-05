figures = figures/Fig1.pdf figures/Fig2.pdf figures/Fig3.pdf \
	figures/Fig4.pdf figures/Fig5.pdf

all : figures.pdf

figures.pdf : $(figures) figures.tex
	pdflatex figures.tex

%.pdf : %.R
	Rscript $< $@

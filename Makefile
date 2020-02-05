figures = figures/Fig1.pdf figures/Fig2.pdf figures/Fig3.pdf \
	figures/Fig4.pdf figures/Fig5.pdf
common = figures/results.R figures/plot.R

all : figures.pdf

figures.pdf : $(figures) figures.tex
	pdflatex figures.tex
	rm *.aux *.log

figures/Fig1.pdf : figures/Fig1.R schematics/Fig1A.pdf schematics/Fig1B.pdf $(common)
	Rscript $< $@

figures/Fig2.pdf : figures/Fig2.R schematics/Fig2A.pdf $(common)
	Rscript $< $@

figures/Fig3.pdf : figures/Fig3.R $(common)
	Rscript $< $@

figures/Fig4.pdf : figures/Fig4.R schematics/Fig4A.pdf $(common)
	Rscript $< $@

figures/Fig5.pdf : figures/Fig5.R $(common)
	Rscript $< $@

clean :
	rm *.pdf figures/*.pdf

figures = figures/Fig1.pdf figures/Fig2.pdf figures/Fig3.pdf \
	figures/Fig4.pdf figures/Fig5.pdf

supplement = figures/Suppl1.pdf figures/Suppl2.pdf figures/Suppl3.pdf \
	figures/Suppl4.pdf figures/Suppl5.pdf figures/Suppl6.pdf

common = figures/plot.R

all : output/figures.pdf output/supplement.pdf
	rm -f *.aux *.log Rplots.pdf

output/figures.pdf : figures.tex $(figures)
	pdflatex $<
	cp figures.pdf output/

output/supplement.pdf : supplement.tex $(supplement)
	pdflatex $<
	cp supplement.pdf output/

figures/%.pdf : figures/%.R $(common) schematics/%A.pdf schematics/%B.pdf
	Rscript $< $@

figures/%.pdf : figures/%.R $(common) schematics/%A.pdf
	Rscript $< $@

figures/%.pdf : figures/%.R $(common)
	Rscript $< $@

clean :
	rm -f *.pdf figures/*.pdf

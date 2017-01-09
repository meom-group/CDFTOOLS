# $Rev$
# $Id$
# $Date$
#-----------------------------------------------------------
TEX=latex
CDFTOOLS=cdftools
CDFTOOLSDIR=CDFTOOLS_2.1

all:  $(CDFTOOLS)_prog.dvi $(CDFTOOLS)_user.dvi index

$(CDFTOOLS).dvi: $(CDFTOOLS).tex

index: $(CDFTOOLS)_prog.dvi $(CDFTOOLS)_user.dvi
	makeindex $(CDFTOOLS)_user.idx 
	makeindex $(CDFTOOLS)_prog.idx 

#pdf: $(CDFTOOLS)_user.pdf $(CDFTOOLS)_prog.pdf

pdf: $(CDFTOOLS)_user.tex $(CDFTOOLS)_prog.tex
	pdflatex $(CDFTOOLS)_user.tex
	pdflatex $(CDFTOOLS)_prog.tex

# to force recompilation of tex file when index is updated
touch:
	touch *.tex

clean:
	\rm -f *.dvi *.log *.aux  *~ *.idx *.ind *.ilg *.toc

commit:
	svn ci 

web:  pdf
	hevea $(CDFTOOLS)_user.tex 
	scp $(CDFTOOLS)_user.html molines@meolipc:/var/www/web/CDFTOOLS/cdftools-2.1.html
	scp $(CDFTOOLS)_user.hind molines@meolipc:/var/www/web/CDFTOOLS/cdftools-2.1.hind
	scp $(CDFTOOLS)_user.pdf molines@meolipc:/var/www/web/CDFTOOLS/cdftools-2.1.pdf

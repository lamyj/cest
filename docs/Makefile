SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

all: html

MEDIA = media/z_spectra.png

media:
	mkdir $@

media/z_spectra.png: ../examples/woessner_2005/figure_2.py
	python3 $<
	mv figure_2.png media/z_spectra.png

html: Makefile media $(MEDIA)
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

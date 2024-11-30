SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

all: html

MEDIA = media/z_spectra.png

media:
	mkdir $@

media/z_spectra.png: ../examples/z_spectra.py
	PYTHONPATH=..:$$PYTHONPATH python3 $<
	mv z_spectra.png media/
	rm mtr.png

html: Makefile media $(MEDIA)
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

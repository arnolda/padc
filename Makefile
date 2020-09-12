# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugänglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

.PHONY: plots check all

all: check data plots padc.pdf

data:
	$(MAKE) -C gendata

plots: data
	$(MAKE) -C plots

check:
	$(MAKE) -C tests check
	$(MAKE) -C examples check

lint:
	$(MAKE) -C examples lint
	for src in *.py; do flake8 $$src; done; exit 0

padc.pdf: padc.tex plots *.py
	latexmk -interaction=nonstopmode

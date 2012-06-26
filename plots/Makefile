# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugänglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

default: \
pendel_energie.pdf  pendel_loesung.pdf \
sinus.pdf runge_lagrange.pdf \
splines.pdf leastsq.pdf \
fftsin.pdf fftalias.pdf \
fourier.pdf wavelet.pdf \
fouriertrafo.pdf \
v_ac.pdf error.pdf \
banach.pdf newton.pdf sekanten.pdf bisektion.pdf \
bessel_diff.pdf num_diff.pdf \
trapezregel.pdf mittelpunktsregel.pdf simpsonregel.pdf \
rechteckregel.pdf trapezzusammen.pdf \
bessel_int.pdf romberg.pdf inf_int.pdf pi.pdf \
quasirandom.pdf \
rejection.pdf inversion.pdf \
statistics_test.pdf optical_test.pdf fft_tests.pdf \
rw.pdf dla.pdf percolation.pdf \
2d-laplace.pdf iterative.pdf

trapezregel.pdf mittelpunktsregel.pdf simpsonregel.pdf rechteckregel.pdf trapezzusammen.pdf: num_int.py
	python $<

pendel_energie.pdf pendel_loesung.pdf: pendel_plot.py
	python $<

%.pdf: %.py
	python $<

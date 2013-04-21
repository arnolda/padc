# check some included code
check:
	python check_horner.py
	python check_random.py
	python check_ortho.py 
	python check_qr.py
	python check_sor.py
	python check_steepest.py
	python check_cg.py
	python check_simplex.py

# output the raw tables that are inserted into the lecture
tables:
	python gen_division.py
	python gen_heron.py
	python gen_random_intro.py

# simulations
simanneal:
	~/progs/Espresso/obj-janeway/Espresso lj_liquid_2d.tcl 100 11.5 "0.1 0.05 0.01 0.001"
vac:
	for dens in 0.5 1; do \
		~/progs/Espresso/obj-janeway/Espresso lj_liquid.tcl 100 $$dens 100000; \
	done

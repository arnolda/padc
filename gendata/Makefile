# check some included code
check:
	python3 check_horner.py
	python3 check_random.py
	python3 check_ortho.py 
	python3 check_qr.py
	python3 check_sor.py
	python3 check_steepest.py
	python3 check_cg.py
	python3 check_simplex.py

# output the raw tables that are inserted into the lecture
tables:
	python3 gen_division.py
	python3 gen_heron.py
	python3 gen_random_intro.py

# simulations
simanneal:
	~/progs/Espresso/obj-janeway/Espresso lj_liquid_2d.tcl 100 11.5 "0.1 0.05 0.01 0.001"
vac:
	for dens in 0.5 1; do \
		~/progs/Espresso/obj-janeway/Espresso lj_liquid.tcl 100 $$dens 100000; \
	done

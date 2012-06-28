# check some included code
check:  check_random
	./check_random
	python check_ortho.py 
	python check_qr.py
	python check_sor.py

# output the raw tables that are inserted into the lecture
tables:
	python gen_division.py
	python gen_heron.py
	python gen_random_intro.py

CXXFLAGS=-O3 -Wall

check_random: check_random.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ `gsl-config --libs`
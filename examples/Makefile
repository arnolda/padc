GRAPH_SOURCES = dla.py pendel.py rng_plot_test.py runge_barycentric.py
SOURCES = fibonacci_qr.py root.py simplex_with_output.py

MARKER = "**********************************************************"

check:
	for src in $(GRAPH_SOURCES); do echo $(MARKER); echo $$src; python3 $$src --limit=10; done
	for src in $(SOURCES); do echo $(MARKER); echo $$src; python3 $$src; done

lint:
	for src in $(GRAPH_SOURCES) $(SOURCES); do flake8 $$src || exit 0; done

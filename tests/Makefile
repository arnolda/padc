MARKER="**********************************************************"
check:
	for check in check_*.py; do echo $(MARKER); echo $$check; python3 $$check || exit 1; done

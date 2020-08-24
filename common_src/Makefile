##tephra2_calc: tephra2_calc.c
##	icc -Wall -Ofast -ffast-math -qopt-report=5 -qopt-report-phase=vec,loop -std=c99 -Wall -o tephra2_calc
tephra2_calc: tephra2_calc.c
	icc -Ofast -xbroadwell -restrict -qopt-report=5 -qopt-report-phase=vec -std=c99 -Wall $^ -o $@

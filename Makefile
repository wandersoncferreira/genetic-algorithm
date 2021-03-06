PROJECT=GA_sem_constrains
#COMPILER=gfortran -g -fbounds-check
#COMPILER=gfortran -ffpe-summary
COMPILER=gfortran -O2 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace
BUILD=$(COMPILER) $(PROJECT).f
RUN=./a.out
all:
		$(BUILD) ga_methods.f convol.f crossover.f ray_extract.f ray_nmo.f ray_rc.f ray_tab.f reverse.f ricker.f shuey.f shuffle.f sort.f sswr.f xcross2d.f zero.f 
		$(RUN)

clean-all:
		rm -f *.txt *.out

clean:
		rm -f *.txt *.out

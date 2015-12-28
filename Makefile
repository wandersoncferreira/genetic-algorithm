PROJECT=GA
COMPILER=gfortran
BUILD=$(COMPILER) $(PROJECT).f
RUN=./a.out
all:
		$(BUILD) ga_methods.f convol.f crossover.f ray_extract.f ray_nmo.f ray_rc.f ray_tab.f reverse.f ricker.f shuey.f shuffle.f sort.f sswr.f xcross2d.f zero.f 
		$(RUN)

clean-all:
		rm -f *.txt *.out

clean:
		rm -f *.txt *.out

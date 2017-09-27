#CC=g++
#CC=clang++
CC=~/llvm2/bin/clang++
OPT=-Wall -std=c++1z -fno-omit-frame-pointer -ffast-math -march=native -ggdb
PROFILE=$(CC) $(OPT) -O3 -DNDEBUG
DEBUG=$(CC) $(OPT)
LIB=-I$(HOME)/include -L$(HOME)/lib 

FILE=input/uniform.1000.20
ifdef NSORT
	DEFINES += -DNSORT
endif
ifdef N_RUNS
	DEFINES += -DN_RUNS=$(N_RUNS)
endif
ifdef N_SAMPLES
	DEFINES += -DN_SAMPLES=$(N_SAMPLES)
endif
ifdef SUBSET_SIZE
	DEFINES += -DSUBSET_SIZE=$(SUBSET_SIZE)
endif

.PHONY: run search debug d_lin lin
run: search $(FILE)
	./search bsEq bsLin_1 bsLin_32 isRecurse isLin_1 isLin_2 oracle < $(FILE)

search: search.cc util.h div.h
	$(PROFILE) $(DEFINES) $< -o $@

debug: search.cc util.h
	$(DEBUG) -DN_RUNS=50 search.cc -o $@
	gdb ./$@

div: div.cc
	$(DEBUG) $< -o $@

asm: search.cc
	$(PROFILE) -g2 -S $< -o $@.s
	$(PROFILE) $@.s -o $@

lin: lin.cc
	#$(PROFILE) -o $@ $(LIB) -fprofile-instr-generate $< -lbenchmark -lpthread 
	#./lin
	#llvm-profdata merge -output=default.profdata default.profraw
	#$(PROFILE) -o $@ $(LIB) -fprofile-instr-use $< -lbenchmark -lpthread
	$(PROFILE) -o $@ $(LIB) $< -lbenchmark -lpthread
	./lin --benchmark_out_format=csv --benchmark_out=lin.csv

d_lin: lin.cc
	$(DEBUG) -o $@ $(LIB) $< -lbenchmark -lpthread && gdb ./$@

clean:
	rm -f ./profile ./debug

# https://www.gnu.org/software/make/manual/make.html#Overriding
DISTRIBUTION=uniform
ARRAY_SIZE=1000
.PHONY: input
input: $(foreach distr,$(DISTRIBUTION),$(foreach sz,$(ARRAY_SIZE),input/$(distr).$(sz)))

# get the suffix and drop the pre-pended dot.
input/%:
	python gendata.py $(subst .,,$(suffix $(basename $*))) $(basename $(basename $*)) > $@

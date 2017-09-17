#CC=g++
#CC=clang++
CC=~/llvm2/bin/clang++
OPT=-Wall -std=c++1z -fno-omit-frame-pointer -ffast-math -march=native -ggdb
PROFILE=$(CC) $(OPT) -O3 -DNDEBUG
DEBUG=$(CC) $(OPT)
BENCHMARKS=-DIS2=1 -DIS=2 -DIS_T=2
LIB=-I$(HOME)/include -L$(HOME)/lib 
IS_T=2
BS_T=2
FILE=input/uniform.1000.20

#	./p2iis 10 -1 < input/um_uniform_1m.txt
#	./r2iis 10 -1 < input/um_uniform_1m.txt
.PHONY: profile debug p2iis r2iis d_lin lin p2ibs r2ibs p2ios p1is p1bs q1is ui uil up ub ubl ril
profile: p2iis r2iis input/uniform.1000.20
	./p2iis 10 -1 < input/uniform.1000.20
	./r2iis 10 -1 < input/uniform.1000.20

p1iu: search.cc util.h
	$(PROFILE) -DNSORT -DIS=1 $< -o $@

p1bu: search.cc util.h
	$(PROFILE) -DNSORT -DBS=1 $< -o $@

p2ios: search.cc util.h
	$(PROFILE) -DN_RUNS=50 -DIS=1 -DOS=2 $< -o $@
	./$@ < input/um_uniform_1m.txt

p2ibu: search.cc util.h
	$(PROFILE) -DNSORT -DIS=1 -DBS=2 $< -o $@

p2ibs: search.cc util.h
	$(PROFILE) -DIS_T=$(IS_T) -DBS_T=$(BS_T) -DN_RUNS=50 -DIS=1 -DBS=2 $< -o $@

r2ibs: search.cc util.h
	$(PROFILE) -DIS=2 -DBS=1 -DIS_T=$(IS_T) -DBS_T=3 $< -o $@
	./$@ < input/uniform.1000.20

p1is: search.cc util.h
	$(PROFILE) -DIS_T=$(IS_T) -DIS=1 $< -o $@

ui: search.cc util.h
	$(PROFILE) -DNSORT -DIS_T=1 -DIS=1 $< -o $@

uil: search.cc util.h
	$(PROFILE) -DN_RUNS=1000000 -DNSORT -DIS_T=2 -DIS=1 $< -o $@

ril: search.cc util.h
	$(PROFILE) -DIS_T=2 -DIS=1 $< -o $@
	./$@ < input/uniform.1000.11

up: search.cc util.h
	$(PROFILE) -DNSORT -DIS_T=3 -DIS=1 $< -o $@

ub: search.cc util.h
	$(PROFILE) -DNSORT -DBS_T=1 -DBS=1 $< -o $@

ubl: search.cc util.h
	$(PROFILE) -DNSORT -DBS_T=2 -DBS=1 $< -o $@

uo: search.cc util.h
	$(PROFILE) -DNSORT -DOS=1 $< -o $@

p1bs: search.cc util.h
	$(PROFILE) -DBS_T=$(BS_T) -DBS=1 $< -o $@

p2iis: search.cc util.h
	$(PROFILE) -DN_RUNS=5000 -DIS_T=2 -DIS=1 -DIS2=2 $< -o $@

p2iiu: search.cc util.h
	$(PROFILE) -DNSORT -DIS=1 -DIS2=2 $< -o $@

r2iis: search.cc util.h
	$(PROFILE) -DN_RUNS=5000 -DIS=2 -DIS2=1 $< -o $@

puk.%: input/uniform.1000.1
	$* < $<

debug: search.cc util.h
	$(DEBUG) $(BENCHMARKS) search.cc -o $@
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

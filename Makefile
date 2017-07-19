#CC=g++
#CC=clang++
CC=~/llvm/bin/clang++
OPT=-Wall -std=c++1z -fno-omit-frame-pointer -ffast-math -march=native -ggdb
PROFILE=$(CC) $(OPT) -O3 -DNDEBUG
DEBUG=$(CC) $(OPT)
BENCHMARKS=-DIS2=1 -DIS=2
LIB=-I$(HOME)/include -L$(HOME)/lib 

.PHONY: profile debug p2iis r2iis d_lin lin p2ibs r2ibs
profile: p2ibs r2ibs input/uniform.1000.20
	./p2ibs < input/uniform.1000.20
	./r2ibs < input/uniform.1000.20

p1iu: search.cc util.h
	$(PROFILE) -DNSORT -DIS=1 $< -o $@

p1bu: search.cc util.h
	$(PROFILE) -DNSORT -DBS=1 $< -o $@

p2ios: search.cc util.h
	$(PROFILE) -DIS=1 -DOS=2 $< -o $@

p2ibu: search.cc util.h
	$(PROFILE) -DNSORT -DIS=1 -DBS=2 $< -o $@

p2ibs: search.cc util.h
	$(PROFILE) -DIS=1 -DBS=2 $< -o $@

r2ibs: search.cc util.h
	$(PROFILE) -DIS=2 -DBS=1 $< -o $@

p1is: search.cc util.h
	$(PROFILE) -DIS=1 $< -o $@

p1bs: search.cc util.h
	$(PROFILE) -DBS=1 $< -o $@

p2iis: search.cc util.h
	$(PROFILE) -DIS=1 -DIS2=2 $< -o $@

p2iiu: search.cc util.h
	$(PROFILE) -DNSORT -DIS=1 -DIS2=2 $< -o $@

r2iis: search.cc util.h
	$(PROFILE) -DIS=2 -DIS2=1 $< -o $@

puk.%: input/uniform.1000.1
	$* < $<


debug: search.cc util.h
	$(DEBUG) $(BENCHMARKS) search.cc -o $@

div: div.cc
	$(DEBUG) $< -o $@

asm: search.cc
	$(PROFILE) -g2 -S $< -o $@.s
	$(PROFILE) $@.s -o $@

lin: lin.cc
	$(PROFILE) -o $@ $(LIB) -fprofile-instr-generate $< -lbenchmark -lpthread 
	./lin
	llvm-profdata merge -output=default.profdata default.profraw
	$(PROFILE) -o $@ $(LIB) -fprofile-instr-use $< -lbenchmark -lpthread
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

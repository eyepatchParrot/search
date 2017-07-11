#GC=g++
#GC=clang++
GC=~/llvm/bin/clang++
OPT=-Wall -std=c++11 -fno-omit-frame-pointer -ffast-math -march=native -ggdb
PROFILE=$(GC) $(OPT) -O3 -DNDEBUG
DEBUG=$(GC) $(OPT)
BENCHMARKS=-DBS=1 -DIS=2

.PHONY: profile
profile: p2iis r2iis input/input.uniform.1000.20
	./p2iis < input/input.uniform.1000.20
	./r2iis < input/input.uniform.1000.20

p1iu: search.cc
	$(PROFILE) -DNSORT -DIS=1 $< -o $@

p1bu: search.cc
	$(PROFILE) -DNSORT -DBS=1 $< -o $@

p2ibu: search.cc
	$(PROFILE) -DNSORT -DIS=1 -DBS=2 $< -o $@

p2ibs: search.cc
	$(PROFILE) -DIS=1 -DBS=2 $< -o $@

p1is: search.cc
	$(PROFILE) -DIS=1 $< -o $@

p1bs: search.cc
	$(PROFILE) -DBS=1 $< -o $@

p2iis: search.cc
	$(PROFILE) -DIS=1 -DIS2=2 $< -o $@

p2iiu: search.cc
	$(PROFILE) -DNSORT -DIS=1 -DIS2=2 $< -o $@

r2iis: search.cc
	$(PROFILE) -DIS=2 -DIS2=1 $< -o $@

puk.%: input/input.uniform.1000.1
	$* < $<


debug: search.cc
	$(DEBUG) $(BENCHMARKS) search.cc -o $@

div: div.cc
	$(DEBUG) $< -o $@

asm: search.cc
	$(PROFILE) -g2 -S $< -o $@.s
	$(PROFILE) $@.s -o $@

clean:
	rm -f ./profile ./debug

# https://www.gnu.org/software/make/manual/make.html#Overriding
DISTRIBUTION=uniform
ARRAY_SIZE=1000
.PHONY: input
input: $(foreach distr,$(DISTRIBUTION),$(foreach sz,$(ARRAY_SIZE),input/input.$(distr).$(sz)))

# get the suffix and drop the pre-pended dot.
input/input.%:
	python gendata.py $(subst .,,$(suffix $(basename $*))) $(basename $(basename $*)) > $@

#GC=g++
#GC=clang++
GC=~/llvm/bin/clang++
OPT=-Wall -std=c++11 -fno-omit-frame-pointer -ffast-math -march=native -ggdb
PROFILE=$(GC) $(OPT) -O3 -DNDEBUG
DEBUG=$(GC) $(OPT)
BENCHMARKS=-DIS=1 -DBS=2 -DOS=3
R_BENCHMARKS=-DIS=2 -DBS=1

.PHONY: profile r_profile
profile: search.cc
	$(PROFILE) $(BENCHMARKS) $< -o $@

r_profile: search.cc
	$(PROFILE) $(R_BENCHMARKS) $< -o $@

p1iu: search.cc
	$(PROFILE) -DIS=1 $< -o $@

p1bu: search.cc
	$(PROFILE) -DBS=1 $< -o $@

p2ibu: search.cc
	$(PROFILE) -DIS=1 -DBS=2 $< -o $@

p2ibs: search.cc
	$(PROFILE) -DSORT -DIS=1 -DBS=2 $< -o $@

p1is: search.cc
	$(PROFILE) -DSORT -DIS=1 $< -o $@

p1bs: search.cc
	$(PROFILE) -DSORT -DBS=1 $< -o $@

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
input: $(foreach distr,$(DISTRIBUTION),$(foreach sz,$(ARRAY_SIZE),input.$(distr).$(sz)))

# get the suffix and drop the pre-pended dot.
input.%:
	python gendata.py $(subst .,,$(suffix $(basename $*))) $(basename $(basename $*)) > $@

#GC=g++
GC=clang++
#GC=~/llvm/bin/clang++
OPT=-Wall -std=c++11 -fno-omit-frame-pointer -ffast-math -march=native
PROFILE=$(GC) $(OPT) -O3 -DNDEBUG
DEBUG=$(GC) $(OPT)

profile: search.cc
	$(PROFILE) $< -o $@

r_profile: search.cc
	$(PROFILE) -DREVERSE $< -o $@

debug: search.cc
	$(DEBUG) search.cc -o $@

div: div.cc
	$(DEBUG) $< -o $@

asm: search.cc
	$(PROFILE) -g2 -S $< -o $@.s
	$(PROFILE) $@.s -o $@

clean:
	rm -f ./profile ./debug

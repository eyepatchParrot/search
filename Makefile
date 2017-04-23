#GC=g++
GC=~/llvm/bin/clang++
OPT=-g2 -Wall -std=c++11 -fno-omit-frame-pointer -ffast-math -march=native
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

clean:
	rm -f ./profile ./debug

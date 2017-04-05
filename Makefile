GC=g++
OPT=-g2 -Wall -std=c++11 -fno-omit-frame-pointer -flto -ffast-math
PROFILE=$(GC) $(OPT) -O3 -DNDEBUG
DEBUG=$(GC) $(OPT)

profile: search.cc
	$(PROFILE) search.cc -o $@

debug: search.cc
	$(DEBUG) search.cc -o $@

clean:
	rm -f ./profile ./debug

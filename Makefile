GC=g++ -g2 -Wall -std=c++11 -fno-omit-frame-pointer -flto -ffast-math

profile: search.cc
	#$(GC) -S -DNDEBUG -Og search.cc
	$(GC) -DNDEBUG -Og search.cc -o $@

debug: search.cc
	$(GC) search.cc -o $@

clean:
	rm -f ./profile ./debug

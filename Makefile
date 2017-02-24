GC=g++ -ggdb -Wall

profile: search.cc
	$(GC) -DNDEBUG -Og search.cc -o $@

debug: search.cc
	$(GC) search.cc -o $@

clean:
	rm -f ./profile ./debug

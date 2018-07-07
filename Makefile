#CXX=clang++
CXXFLAGS=-W -Wall -pedantic -std=c++14 -MD -MP

run_main: main.pass

%.pass: %
	./$*
	touch $@
	
main: main.o random.o
	$(CXX) -o $@ $^

clean:
	rm -f *.o *.d *.pass

-include *.d

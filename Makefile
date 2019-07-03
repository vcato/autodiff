CXX=clang++
CXXFLAGS=-W -Wall -pedantic -std=c++14 -MD -MP

all: main.pass kinematic2d.pass

%.pass: %
	./$*
	touch $@
	
main: main.o random.o
	$(CXX) -o $@ $^

kinematic2d: kinematic2d.o random.o
	$(CXX) -o $@ $^

clean:
	rm -f *.o *.d *.pass

-include *.d

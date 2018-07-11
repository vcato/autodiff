CXXFLAGS=-O2 -g -W -Wall -pedantic -std=c++14 -MD -MP

perf_test_disassembly.txt: perf_test.o
	#objdump -S --disassemble $^ >$@
	#$(CXX) $(CXXFLAGS) -S -fverbose-asm -Wa,-adhln perf_test.cpp
	$(CXX) $(CXXFLAGS) -c -g -Wa,-a,-ad perf_test.cpp > $@

perf_test.s: perf_test.o
	#objdump -S --disassemble $^ >$@
	$(CXX) $(CXXFLAGS) -S -fverbose-asm -Wa,-adhln perf_test.cpp

perf_test: perf_test.o
	$(CXX) -o $@ $^

clean:
	rm -f *.o *.d perf_test perf_test_disassembly.txt

-include *.d

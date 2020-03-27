CXX = g++ -std=c++17 -stdlib=libc++ 
CXXFLAGS = -g -lgtest -lgtest_main -lpthread
INCS = -I./ -I../../src -I/opt/gtest/include /usr/local/lib/libgtest_main.a
OBJS = src/blackscholes.o src/newton.o

all: $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) -o  run src/main.cpp $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCS)

clean:
	rm run 
	#rm -rf test/testAll *.o test/testAll.dSYM
	#rm test/*.o 
	rm src/*.o 
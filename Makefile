CC = g++

CXXFLAGS = -O3 -pthread -std=c++11 -Wall -Wextra -pedantic 

default:
	$(CC) $(CXXFLAGS) -o jones src/jones.cc
	$(CC) $(CXXFLAGS) -o motzkin src/motzkin.cc
	$(CC) $(CXXFLAGS) -o kauffman src/kauffman.cc

jones:
	$(CC) $(CXXFLAGS) -o jones src/jones.cc

kauffman:
	$(CC) $(CXXFLAGS) -o kauffman src/kauffman.cc

motzkin:
	$(CC) $(CXXFLAGS) -o motzkin src/motzkin.cc

test: 
	tst/jones.sh
	tst/motzkin.sh
	tst/kauffman.sh

clean:
	rm -f jones
	rm -f motzkin
	rm -f kauffman

.PHONY: default jones kauffman motzkin

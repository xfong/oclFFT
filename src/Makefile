CC=gcc
CXX=g++
CFLAGS=-Wall -I./ -std=c++11
LDFLAGS=

header.o: header.cxx header.h
	$(CXX) $(CFLAGS) -c -o $@ $<

plan.o: plan.cxx plan.h
	$(CXX) $(CFLAGS) -c -o $@ $<

test_print.o: test_print.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

main.o: main.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

main: main.o plan.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	rm -f *.o

test_print: test_print.o header.o
	$(CXX) $(CFLAGS) -o $@ $^
	rm -f *.o

all: main test_print


.PHONY: clean

clean:
	rm -f *.o *~ *.exe


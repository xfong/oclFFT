CC=gcc
CXX=g++
CFLAGS=-Wall -I./ -std=c++11
LDFLAGS=

kernel_header.o: kernel_header.cxx kernel_header.h
	$(CXX) $(CFLAGS) -c -o $@ $<

plan.o: plan.cxx plan.h kernel_header.h
	$(CXX) $(CFLAGS) -c -o $@ $<

test_print.o: test_print.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

main.o: main.cxx
	$(CXX) $(CFLAGS) -c -o $@ $<

main: main.o plan.o kernel_header.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_print: test_print.o kernel_header.o
	$(CXX) $(CFLAGS) -o $@ $^

all: main test_print


.PHONY: clean

clean:
	rm -f *.o *~ *.exe

objclean:
	rm -f *.o *~


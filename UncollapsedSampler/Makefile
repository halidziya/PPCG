CC=cc
CXX=g++

CFLAGS=-c -Wall

cpu : UncollapsedSampler

gpu : UncollapsedSampler
	qsub compile.pbs

UncollapsedSampler : main.cpp
	$(CXX) $(CFLAGS) -o UncollapsedSampler main.cpp

clean:
	rm UncollapsedSampler

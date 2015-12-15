CC=cc
CXX=g++

CFLAGS=-c -Wall

cpu : UncollapsedSampler

gpu : 
	qsub compile.pbs

gpu_compile :
	cc main.cpp  -O3 -o UncollapsedSampler.gpu -hlist=a

UncollapsedSampler : main.cpp
	$(CXX) $(CFLAGS) -o UncollapsedSampler main.cpp

clean:
	rm GPU_compile*
	rm UncollapsedSampler
	rm main.lst*
	rm main_1*

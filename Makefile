# Makefile

all: singlet tryplet coupled

singlet:
	nvcc examples/singlet/singlet.cu -o singlet -O3

tryplet:
	nvcc examples/tryplet/tryplet.cu -o tryplet -O3 

coupled:
	nvcc examples/coupled/coupled.cu -o coupled -O3 


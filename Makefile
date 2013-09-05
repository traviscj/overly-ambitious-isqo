

all:
	cd src && make

doxy:
	cd include && doxygen default-doxy

wc:
	wc -l include/*.hh src/*.cc

clean:
	cd src && make clean

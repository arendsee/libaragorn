all:
	g++ -O -o test -L${PWD}/.. -Wl,-rpath=${PWD}/.. -laragorn -I${PWD}/../src test.c    

test:
	bash test.sh

clean:
	rm -f test-* test

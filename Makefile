VERSION=1.0.0

all:
	g++ -fPIC -shared -Wall -O3 -std=c++11 -o libaragorn.so.${VERSION} -I. *.cpp
	cp -s libaragorn.so.1.0.0 libaragorn.so

test:
	g++ -Wall -O3 -o test *.[ch]
	./test.sh

clean:
	rm -f *.gch test test-trna* a.out

VERSION=1.0.0

all:
	g++ -fPIC -shared -Wall -O -std=c++11 -o libaragorn.so.${VERSION} -I. *.cpp
	cp -fs libaragorn.so.1.0.0 libaragorn.so

clean:
	rm -f *.gch test test-trna* a.out

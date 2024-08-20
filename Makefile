VERSION=1.0.0

all:
	g++ -fPIC -shared -O -std=c++11 -o libaragorn.so.${VERSION} -Isrc/ src/*.cpp
	cp -fs libaragorn.so.1.0.0 libaragorn.so

clean:
	rm -rf src/*.gch libaragorn.so.1.0.0 libaragorn.so

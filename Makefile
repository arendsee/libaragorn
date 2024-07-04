all:
	g++ -Wall -O -o test *.[ch]
	./test
	rm -f *.gch test

SHELL=bash
COMP=g++
CCOMP=gcc
PREFIX ?=/usr/local
FLAGS=-std=c++11 --std=gnu++11 -fPIC
IFLAGS=-I$(PREFIX)/include
LFLAGS=-L$(PREFIX)/lib

all: lib/libmixturedist.so lib/libmixturedist.a

lib/libmixturedist.so: build/mixtureModel.o build/mixtureDist.o build/functions.o build/cdflib.o build/incbeta.o
	$(CCOMP) $(IFLAGS) $(LFLAGS) -shared -o lib/libmixturedist.so build/cdflib.o build/functions.o build/incbeta.o build/mixtureDist.o build/mixtureModel.o -lstdc++

lib/libmixturedist.a: build/mixtureModel.o build/mixtureDist.o build/functions.o build/cdflib.o build/incbeta.o
	ar rcs lib/libmixturedist.a build/mixtureModel.o build/mixtureDist.o build/functions.o build/cdflib.o build/incbeta.o
 	
build/mixtureModel.o: src/mixtureModel.cpp src/mixtureDist.h build/mixtureDist.o
	$(COMP) $(FLAGS) $(IFLAGS) -c src/mixtureModel.cpp -o build/mixtureModel.o

build/mixtureDist.o: src/mixtureDist.cpp src/functions.h build/functions.o
	$(COMP) $(FLAGS) $(IFLAGS) -c src/mixtureDist.cpp -o build/mixtureDist.o

build/functions.o: src/functions.cpp src/functions.h src/incbeta.h src/cdflib.h build/cdflib.o build/incbeta.o
	$(COMP) $(FLAGS) $(IFLAGS) -c src/functions.cpp -o build/functions.o

build/cdflib.o: src/cdflib.c src/cdflib.h
	$(CCOMP) -fPIC -c src/cdflib.c -o build/cdflib.o

build/incbeta.o: src/incbeta.c src/incbeta.h
	$(CCOMP) -Wall -Wshadow -O2 -fPIC -g -c src/incbeta.c -o build/incbeta.o

clean:
	rm build/*.o
	rm lib/*.so

install: | $(PREFIX)/lib $(PREFIX)/include/mixtureDist
	cp lib/*.so $(PREFIX)/lib
	cp lib/*.a $(PREFIX)/lib
	cp src/mixtureDist.h $(PREFIX)/include/mixtureDist
	cp src/mixtureModel.h $(PREFIX)/include/mixtureDist
	cp src/functions.h $(PREFIX)/include/mixtureDist

$(PREFIX)/lib:
	mkdir -p $(PREFIX)/lib

$(PREFIX)/include/mixtureDist:
	mkdir -p $(PREFIX)/include/mixtureDist

SHELL=/bin/bash
COMP=g++
CCOMP=gcc
FLAGS=-std=c++11 --std=gnu++11 -fPIC
PREFIX ?=/usr/local

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CELLAR=$(shell brew info argp-standalone | grep Cellar | cut -d' ' -f1)
	FLAGS += -I$(CELLAR)/include/
	LDFLAGS += -L$(CELLAR)/lib/ -largp
endif

all: lib/libmixturedist.so

lib/libmixturedist.so: build/mixtureModel.o build/mixtureDist.o build/functions.o build/cdflib.o build/incbeta.o
	$(CCOMP) -shared -o lib/libmixturedist.so build/cdflib.o build/functions.o build/incbeta.o build/mixtureDist.o build/mixtureModel.o -lstdc++

build/mixtureModel.o: src/mixtureModel.cpp src/mixtureDist.h build/mixtureDist.o
	$(COMP) $(FLAGS) -c src/mixtureModel.cpp -o build/mixtureModel.o

build/mixtureDist.o: src/mixtureDist.cpp src/functions.h build/functions.o
	$(COMP) $(FLAGS) -c src/mixtureDist.cpp -o build/mixtureDist.o

build/functions.o: src/functions.cpp src/functions.h src/incbeta.h src/cdflib.h build/cdflib.o build/incbeta.o
	$(COMP) $(FLAGS) -c src/functions.cpp -o build/functions.o

build/cdflib.o: src/cdflib.c src/cdflib.h
	$(CCOMP) -fPIC -c src/cdflib.c -o build/cdflib.o

build/incbeta.o: src/incbeta.c src/incbeta.h
	$(CCOMP) -Wall -Wshadow -O2 -fPIC -g -c src/incbeta.c -o build/incbeta.o

clean:
	rm build/*.o
	rm lib/libmixturedist.so

install: | $(PREFIX)/lib $(PREFIX)/include/mixtureDist
	cp lib/libmixturedist.so $(PREFIX)/lib
	cp src/mixtureDist.h $(PREFIX)/include/mixtureDist
	cp src/mixtureModel.h $(PREFIX)/include/mixtureDist
	cp src/functions.h $(PREFIX)/include/mixtureDist

$(PREFIX)/lib:
	mkdir -p $(PREFIX)/lib

$(PREFIX)/include/mixtureDist:
	mkdir -p $(PREFIX)/include/mixtureDist

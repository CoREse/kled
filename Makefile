CC=g++
AR=ar
#CPPFLAGS= -Wall -g -Lhtslib -lhts -Lcrelib -lcre optutils/opthelper.a --std=c++17
CPPFLAGS= --std=c++17 -O3
LDFLAGS= -fopenmp 
HTSLIBDIR=htslib
SUBMODULES=crelib $(HTSLIBDIR) optutils

PROJECT_OBJS=kled.o input.o signature.o contig.o StatsManager.o StatsTracker.o clustering.o report.o
PROJECT_HEADERS=kled.h input.h signature.h contig.h clustering.h report.h
EXAMPLE_OBJS=
HTSLIB=$(HTSLIBDIR)/libhts.a
LAUNCHER=./launcher

ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

HTSPREFIX=htslib/

.PHONY: all clean $(SUBMODULES) htsconf

all: build

build: $(SUBMODULES) kled

$(SUBMODULES):
	make -C $@

htsconf:
	cd $(HTSLIBDIR) && autoreconf -i && ./configure

$(PROJECT_OBJS): $(PROJECT_HEADERS)

HTSLIB_LIBS = -lz -lm -lbz2 -llzma -lcurl -lpthread -lcrypto -ldeflate
kled: $(PROJECT_OBJS)
	$(CC) $^ -o $@ -fopenmp $(HTSPREFIX)/libhts.a $(HTSLIB_LIBS) -Lcrelib -lcre optutils/opthelper.a --std=c++17

clean:
	rm *.o kled

install: kled
	install -d $(PREFIX)/bin/
	install -m 755 kled $(PREFIX)/bin/
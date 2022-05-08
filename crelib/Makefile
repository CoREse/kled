CC=g++
AR=ar
CPPFLAGS= -Wall -O3
LDFLAGS=
LIBS=

LIBCRE_OBJS=crelib.o
LIBCRE_HEADERS=crelib.h defines.h
SAMPLE_OBJS=samples/sample.o
LIBCRE=libcre.a

libs:$(LIBCRE)

all: $(LIBCRE) sample

$(LIBCRE):$(LIBCRE_OBJS)
	$(AR) -rc $@ $(LIBCRE_OBJS)

sample: $(SAMPLE_OBJS) $(LIBCRE)
	$(LINK.cpp) $^ -o $@

$(LIBCRE_OBJS):$(LIBCRE_HEADERS)

clean:
	rm *.o *.a samples/*.o sample

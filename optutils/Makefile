CC=g++
AR=ar
CPPFLAGS= -Wall -O3
LDFLAGS=
LIBS=

OPTHELPER_OBJS=OptHelper.o
OPTHELPER_HEADERS=OptHelper.h
EXAMPLE_OBJS=examples/example.o
OPTHELPER=opthelper.a

libs:$(OPTHELPER)

all: $(OPTHELPER) example

$(OPTHELPER):$(OPTHELPER_OBJS)
	$(AR) -rc $@ $(OPTHELPER_OBJS)

example: $(EXAMPLE_OBJS) $(OPTHELPER)
	$(LINK.c) $^ -o $@

$(OPTHELPER_OBJS):$(OPTHELPER_HEADERS)

clean:
	rm *.o *.a examples/*.o example

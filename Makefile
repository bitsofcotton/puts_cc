CXX=	clang++

# compiler flags.
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-L/usr/local/lib -lomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
LDFLAGS+=	-static

#CXXFLAGS+=	-D_FLOAT_BITS_=32
#CXXFLAGS+=	-D_FLOAT_BITS_=64
#CXXFLAGS+=	-D_FLOAT_BITS_=128
#CXXFLAGS+=	-D_FLOAT_BITS_=256
#CXXFLAGS+=	-D_FLOAT_BITS_=512

CLEANFILES= *.o tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}

tools.o:        tools.cc corpus.hh lienn.hh


CXX=	clang++

# compiler flags.
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-L/usr/local/lib -lomp -fopenmp
#CXXFLAGS+=	-pg
CXXFLAGS+=	-std=c++11
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lgomp -fopenmp
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib
#LDFLAGS+=	-static

#CXXFLAGS+=	-D_FLOAT_BITS_=32
#CXXFLAGS+=	-D_FLOAT_BITS_=64
#CXXFLAGS+=	-D_FLOAT_BITS_=128
#CXXFLAGS+=	-D_FLOAT_BITS_=256
#CXXFLAGS+=	-D_FLOAT_BITS_=512

CLEANFILES= *.o puts putsmp putsc putscmp

all:	puts putsmp putsc putscmp

clean:
	@rm -rf ${CLEANFILES}

puts:
	${CXX} ${CXXFLAGS} -static -o puts tools.cc
putsmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o putsmp tools.cc
putsc:
	${CXX} ${CXXFLAGS} -D_CONTINUOUS_ -static -o putsc tools.cc
putscmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONTINUOUS_ -o putscmp tools.cc


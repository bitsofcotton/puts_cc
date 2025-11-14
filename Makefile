CXX=	clang++
#CXX=	c++

# compiler flags.
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -g3
#CXXFLAGS+=	-O2 -g3
#CXXFLAGS+=	-pg
CXXFLAGS+=	-std=c++11
#CXXFLAGS+=	-std=gnu++98
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

# lieonn.hh flags
# N.B. sed -e s/static\ inline//g | sed -e s/inline//g
#CXXFLAGS+=     -D_OLDCPP_ -ftemplate-depth-99

CLEANFILES= *.o puts putsmp puts32 puts32mp

all:	puts putsmp puts32 puts32mp

clean:
	@rm -rf ${CLEANFILES}

puts:
	${CXX} ${CXXFLAGS} -static -o puts tools.cc
putsmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o putsmp tools.cc
puts32:
	${CXX} ${CXXFLAGS} -D_FLOAT_BITS_=32 -static -o puts32 tools.cc
puts32mp:
	${CXX} ${CXXFLAGS} -D_FLOAT_BITS_=32 ${MPFLAGS} -o puts32mp tools.cc


CXX=	clang++

# compiler flags.
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -g3
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

CLEANFILES= *.o puts puts3 putsmp puts3mp puts32 puts3-32 puts32mp puts3-32mp

all:	puts puts3 putsmp puts3mp puts32 puts3-32 puts32mp puts3-32mp

clean:
	@rm -rf ${CLEANFILES}

puts:
	${CXX} ${CXXFLAGS} -static -o puts tools.cc
puts3:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=3 -o puts3 tools.cc
putsmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o putsmp tools.cc
puts3mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=3 -o puts3mp tools.cc
puts32:
	${CXX} ${CXXFLAGS} -D_FLOAT_BITS_=32 -static -o puts32 tools.cc
puts3-32:
	${CXX} ${CXXFLAGS} -D_PREDV_=3 -D_FLOAT_BITS_=32 -static -o puts3-32 tools.cc
puts32mp:
	${CXX} ${CXXFLAGS} -D_FLOAT_BITS_=32 ${MPFLAGS} -o puts32mp tools.cc
puts3-32mp:
	${CXX} ${CXXFLAGS} -D_PREDV_=3 -D_FLOAT_BITS_=32 ${MPFLAGS} -o puts3-32mp tools.cc


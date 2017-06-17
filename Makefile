# gcc.
CXX=    g++

# compiler flags.
#CXXFLAGS=     -Ofast -mavx
CXXFLAGS=      -O2 -g2 -mavx
LDFLAGS=        -lstdc++

CLEANFILES= *.o *.dSYM tools

all:    tools
tools:  corpus.o tools.o

clean:
        @rm -rf ${CLEANFILES}


/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/*
BSD 3-Clause License

Copyright (c) 2013-2021, bitsofcotton (kazunobu watatsu)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#if !defined(_SIMPLELIN_)

using std::move;
using std::max;
using std::min;
using std::swap;
using std::sort;
using std::binary_search;
using std::pair;
using std::make_pair;
using std::map;
using std::vector;
using std::stringstream;
using std::ifstream;
using std::ofstream;

using std::string;
using std::cerr;
using std::endl;
using std::flush;

#if !defined(_FLOAT_BITS_)
  #include <complex>
  #include <cmath>
  using namespace std;
  typedef uint64_t myuint;
  typedef int64_t  myint;
  typedef long double myfloat;
#else

// Double int to new int class.
template <typename T, int bits> class DUInt {
public:
  inline DUInt();
  inline DUInt(const int& src);
  inline DUInt(const T& src);
  inline DUInt(const DUInt<T,bits>& src);
  inline DUInt(const DUInt<DUInt<T,bits>,bits*2>& src);
  inline DUInt(DUInt<T,bits>&& src);
  inline ~DUInt();
  
  inline DUInt<T,bits>& operator ++ ();
  inline DUInt<T,bits>  operator ++ (int32_t);
  inline DUInt<T,bits>& operator -- ();
  inline DUInt<T,bits>  operator -- (int32_t);
  inline DUInt<T,bits>  operator -  () const;
  inline DUInt<T,bits>  operator +  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator += (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator -  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator -= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator *  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator *= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator /  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator /= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator %  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator %= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator << ( const int& b)            const;
  inline DUInt<T,bits>& operator <<= (const int& b);
  inline DUInt<T,bits>  operator >> ( const int& b)            const;
  inline DUInt<T,bits>& operator >>= (const int& b);
  inline DUInt<T,bits>  operator &  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator &= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator |  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator |= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator ^  (const DUInt<T,bits>& src) const;
  inline DUInt<T,bits>& operator ^= (const DUInt<T,bits>& src);
  inline DUInt<T,bits>  operator ~  ()                         const;
  inline DUInt<T,bits>& operator =  (const DUInt<T,bits>& src);
  inline DUInt<T,bits>& operator =  (const DUInt<DUInt<T,bits>,bits*2>& src);
  inline DUInt<T,bits>& operator =  (DUInt<T,bits>&& src);
  inline bool           operator <  (const DUInt<T,bits>& src) const;
  inline bool           operator <= (const DUInt<T,bits>& src) const;
  inline bool           operator >  (const DUInt<T,bits>& src) const;
  inline bool           operator >= (const DUInt<T,bits>& src) const;
  inline bool           operator == (const DUInt<T,bits>& src) const;
  inline bool           operator != (const DUInt<T,bits>& src) const;
  inline bool           operator && (const DUInt<T,bits>& src) const;
  inline bool           operator || (const DUInt<T,bits>& src) const;
  inline bool           operator !    () const;
  inline                operator bool () const;
  inline                operator int  () const;
  inline                operator T    () const;
  inline                operator DUInt<T,bits> () const;

  // friend std::ostream&  operator << (std::ostream& os, DUInt<T,bits>  v);
  // friend std::istream&  operator >> (std::istream& is, DUInt<T,bits>& v);

  T e[2];
};

template <typename T, int bits> inline DUInt<T,bits>::DUInt() {
  assert(0 < bits && ! (bits & 3));
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(const int& src) {
  const auto abssrc(src < 0 ? - src : src);
  e[0]   = T(abssrc);
  e[1]  ^= e[1];
  if(abssrc != src)
    *this = - *this;
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(const T& src) {
  const auto abssrc(src < T(0) ? - src : src);
  e[0]   = abssrc;
  e[1]  ^= e[1];
  if(abssrc != src)
    *this = - *this;
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(const DUInt<T,bits>& src) {
  *this = src;
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(const DUInt<DUInt<T,bits>,bits*2>& src) {
  *this = src;
}

template <typename T, int bits> inline DUInt<T,bits>::DUInt(DUInt<T,bits>&& src) {
  *this = src;
}

template <typename T, int bits> inline DUInt<T,bits>::~DUInt() {
  ;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator ++ () {
  ++ e[0];
  if(!e[0])
    ++ e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator ++ (int32_t) {
  const auto work(*this);
  ++ *this;
  return work;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator -- () {
  if(!e[0])
    -- e[1];
  -- e[0];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator -- (int32_t) {
  const auto work(*this);
  -- *this;
  return work;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator -  () const {
  auto work(~ *this);
  return ++ work;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator +  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work += src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator += (const DUInt<T,bits>& src) {
  // N.B. assembler can boost dramatically this code. but not here.
  const auto e0(max(e[0], src.e[0]));
  e[0] += src.e[0];
  if(e[0] < e0)
    e[1] ++;
  e[1] += src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator -  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work -= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator -= (const DUInt<T,bits>& src) {
  return *this += - src;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator *  (const DUInt<T,bits>& src) const {
  DUInt<T,bits> result;
  result ^= result;
  for(int i = 0; i < 2 * bits; i ++)
    if(int(src >> i) & 1)
      result += *this << i;
  // N.B.
  //   If we work with multiply with table and summing up with simple window,
  //   and the parallel condition, we can reduce better:
  //     with bit pair [x1, x2, ..., xn], [y1, y2, ..., yn],
  //     make table [[x1*y1, ..., x1*yn],...,[xn*y1, ... xn*yn]],
  //     then counter orthogonal sum-up with parallel.
  //     we get [z1, ... zn] == [x1*y1, ..., sum_i+j=k(x_i*y_j), ..., xn*yn],
  //     then sum-up with certain bit adder and fixing one by one:
  //     r1 := x1*y1, s1 := ((x1*y1) >> 1) + z2, r2 := s2 & 1, ... and so on.
  return result;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator *= (const DUInt<T,bits>& src) {
  return *this = *this * src;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator /  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work /= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator /= (const DUInt<T,bits>& src) {
  const static DUInt<T,bits> one(1);
  if(! src)
    throw "Zero division";
  if(! *this)
    return *this;
  auto cache(*this);
  *this ^= *this;
  for(int i = 2 * bits - 1; 0 <= i; i --)
    if((cache >> i) >= src) {
      *this |= one << i;
      cache -= src << i;
    }
  assert(cache < src);
  return *this;
  // N.B. if we works with newton's method, better speed will be gained.
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator %  (const DUInt<T,bits>& src) const {
  return *this - ((*this / src) * src);
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator %= (const DUInt<T,bits>& src) {
  return *this = *this % src;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator << (const int& b)             const {
  auto work(*this);
  return work <<= b;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator <<= (const int& b) {
  if(! b)
    return *this;
  else if(b < 0)
    return *this >>= (- b);
  else if(b >= bits * 2)
    return *this ^= *this;
  else if(b > bits) {
    e[1]  = e[0] << (b - bits);
    e[0] ^= e[0];
  } else if(b == bits) {
    e[1]  = e[0];
    e[0] ^= e[0];
  } else {
    e[1] <<= b;
    e[1]  |= e[0] >> (bits - b);
    e[0] <<= b;
  }
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator >> (const int& b)             const {
  auto work(*this);
  return work >>= b;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator >>= (const int& b) {
  if(! b)
    return *this;
  else if(b < 0)
    return *this <<= (- b);
  else if(b >= bits * 2)
    return *this ^= *this;
  else if(b > bits) {
    e[0]  = e[1] >> (b - bits);
    e[1] ^= e[1];
  } else if(b == bits) {
    e[0]  = e[1];
    e[1] ^= e[1];
  } else {
    e[0] >>= b;
    e[0]  |= e[1] << (bits - b);
    e[1] >>= b;
  }
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator &  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work &= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator &= (const DUInt<T,bits>& src) {
  e[0] &= src.e[0];
  e[1] &= src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator |  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work |= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator |= (const DUInt<T,bits>& src) {
  e[0] |= src.e[0];
  e[1] |= src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator ^  (const DUInt<T,bits>& src) const {
  auto work(*this);
  return work ^= src;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator ^= (const DUInt<T,bits>& src) {
  e[0] ^= src.e[0];
  e[1] ^= src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator ~  () const {
  DUInt<T,bits> work;
  work.e[0] = ~ e[0];
  work.e[1] = ~ e[1];
  return work;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator =  (const DUInt<T,bits>& src) {
  e[0] = src.e[0];
  e[1] = src.e[1];
  return *this;
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator =  (const DUInt<DUInt<T,bits>,bits*2>& src) {
  return *this = src.e[0];
}

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator =  (DUInt<T,bits>&& src) {
  e[0] = move(src.e[0]);
  e[1] = move(src.e[1]);
  return *this;
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator <  (const DUInt<T,bits>& src) const {
  if(e[1])
    return e[1] != src.e[1] ? e[1] < src.e[1] : e[0] < src.e[0];
  return bool(src.e[1]) || e[0] < src.e[0];
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator <= (const DUInt<T,bits>& src) const {
  return *this < src || *this == src;
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator >  (const DUInt<T,bits>& src) const {
  return ! (*this <= src);
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator >= (const DUInt<T,bits>& src) const {
  return ! (*this < src);
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator == (const DUInt<T,bits>& src) const {
  return ! (*this != src);
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator != (const DUInt<T,bits>& src) const {
  return (*this ^ src).operator bool();
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator && (const DUInt<T,bits>& src) const {
  return this->operator bool() && src.operator bool();
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator || (const DUInt<T,bits>& src) const {
  return this->operator bool() || src.operator bool();
}

template <typename T, int bits> inline bool      DUInt<T,bits>::operator !    () const {
  return ! this->operator bool();
}

template <typename T, int bits> inline           DUInt<T,bits>::operator bool () const {
  return e[0] || e[1];
}

template <typename T, int bits> inline           DUInt<T,bits>::operator int () const {
  return int(e[0]);
}

template <typename T, int bits> inline           DUInt<T,bits>::operator T   () const {
  return e[0];
}

template <typename T, int bits> inline           DUInt<T,bits>::operator DUInt<T,bits> () const {
  return *this;
}

template <typename T, int bits> std::ostream&  operator << (std::ostream& os, DUInt<T,bits> v) {
  const static DUInt<T,bits> ten(10);
  vector<char> buf;
  while(v) {
    const auto div(v / ten);
    buf.emplace_back(int(v - div * ten));
    v = div;
  }
  if(buf.size()) {
    for(int i = 0; 0 <= i && i < buf.size(); i ++)
      os << int(buf[buf.size() - 1 - i]);
    return os;
  }
  return os << '0';
}

template <typename T, int bits> std::istream&  operator >> (std::istream& is, DUInt<T,bits>& v) {
  const static DUInt<T,bits> ten(10);
  v = DUInt<T,bits>(0);
  // skip white spaces.
  while(! is.eof()) {
    const auto buf(is.get());
    if(buf != ' ' && buf != '\t') {
      is.unget();
      break;
    }
  }
  while(! is.eof() ) {
    const auto buf(is.get());
    if('0' <= buf && buf <= '9') {
      v *= ten;
      v += DUInt<T,bits>(int(buf - '0'));
    } else {
      is.unget();
      break;
    }
  }
  return is;
}


// add sign.
template <typename T, int bits> class Signed : public T {
public:
  inline Signed();
  inline Signed(const int& src);
  inline Signed(const T& src);
  inline Signed(const Signed<T,bits>& src);
  inline bool operator <  (const Signed<T,bits>& src) const;
  inline bool operator <= (const Signed<T,bits>& src) const;
  inline bool operator >  (const Signed<T,bits>& src) const;
  inline bool operator >= (const Signed<T,bits>& src) const;
  // friend std::ostream&  operator << (std::ostream& os, Signed<T,bits> v);
};

template <typename T, int bits> inline Signed<T,bits>::Signed() {
  ;
}

template <typename T, int bits> inline Signed<T,bits>::Signed(const int& src) {
  T tsrc(src);
  *this = reinterpret_cast<const Signed<T,bits>&>(tsrc);
}

template <typename T, int bits> inline Signed<T,bits>::Signed(const T& src) {
  *this = reinterpret_cast<const Signed<T,bits>&>(src);
}

template <typename T, int bits> inline Signed<T,bits>::Signed(const Signed<T,bits>& src) {
  *this = src;
}

template <typename T, int bits> inline bool Signed<T,bits>::operator <  (const Signed<T,bits>& src) const {
  const auto mthis(int(*this >> (bits - 1)));
  const auto msrc( int(src   >> (bits - 1)));
  if(mthis ^ msrc)
    return mthis;
  if(mthis)
    return - dynamic_cast<const T&>(src) < - dynamic_cast<const T&>(*this);
  return dynamic_cast<const T&>(*this) < dynamic_cast<const T&>(src);
}

template <typename T, int bits> inline bool Signed<T,bits>::operator <=  (const Signed<T,bits>& src) const {
  return ! (*this > src);
}

template <typename T, int bits> inline bool Signed<T,bits>::operator >  (const Signed<T,bits>& src) const {
  return ! (*this < src) && *this != src;
}

template <typename T, int bits> inline bool Signed<T,bits>::operator >= (const Signed<T,bits>& src) const {
  return ! (*this < src);
}

template <typename T, int bits> std::ostream& operator << (std::ostream& os, Signed<T,bits> v) {
  const static Signed<T,bits> zero(0);
  if(v < zero) {
    os << '-';
    v = - v;
  }
  return os << dynamic_cast<const T&>(v);
}


// integer to integer float part.
template <typename T, typename W, int bits, typename U> class SimpleFloat {
public:
  inline SimpleFloat();
  template <typename V> inline SimpleFloat(const V& src);
  inline SimpleFloat(const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat(SimpleFloat<T,W,bits,U>&& src);
  inline ~SimpleFloat();
  
  inline SimpleFloat<T,W,bits,U>  operator -  () const;
  inline SimpleFloat<T,W,bits,U>  operator +  (const SimpleFloat<T,W,bits,U>& src) const;
         SimpleFloat<T,W,bits,U>& operator += (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator -  (const SimpleFloat<T,W,bits,U>& src) const;
  inline SimpleFloat<T,W,bits,U>& operator -= (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator *  (const SimpleFloat<T,W,bits,U>& src) const;
         SimpleFloat<T,W,bits,U>& operator *= (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator /  (const SimpleFloat<T,W,bits,U>& src) const;
         SimpleFloat<T,W,bits,U>& operator /= (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator %  (const SimpleFloat<T,W,bits,U>& src) const;
  inline SimpleFloat<T,W,bits,U>& operator %= (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>  operator <<  (const U& b) const;
  inline SimpleFloat<T,W,bits,U>& operator <<= (const U& b);
  inline SimpleFloat<T,W,bits,U>  operator >>  (const U& b) const;
  inline SimpleFloat<T,W,bits,U>& operator >>= (const U& b);
  inline SimpleFloat<T,W,bits,U>& operator =  (const SimpleFloat<T,W,bits,U>& src);
  inline SimpleFloat<T,W,bits,U>& operator =  (SimpleFloat<T,W,bits,U>&& src);
  inline bool             operator == (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator != (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator <  (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator <= (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator >  (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator >= (const SimpleFloat<T,W,bits,U>& src) const;
  inline bool             operator !  () const;
  inline                  operator bool () const;
  inline                  operator int  () const;
  inline                  operator T    () const;
  inline                  operator SimpleFloat<T,W,bits,U> () const;
  // XXX: absfloor, absceil implementation.
  inline SimpleFloat<T,W,bits,U>  floor() const;
  inline SimpleFloat<T,W,bits,U>  ceil() const;
  inline SimpleFloat<T,W,bits,U>  abs()  const;
         SimpleFloat<T,W,bits,U>  log()  const;
         SimpleFloat<T,W,bits,U>  exp()  const;
         SimpleFloat<T,W,bits,U>  sin()  const;
         SimpleFloat<T,W,bits,U>  cos()  const;
         SimpleFloat<T,W,bits,U>  atan() const;
  inline SimpleFloat<T,W,bits,U>  sqrt() const;
  
  // friend std::ostream&    operator << (std::ostream& os, const SimpleFloat<T,W,bits,U>& v);
  // friend std::istream&    operator >> (std::istream& is, SimpleFloat<T,W,bits,U>& v);
  
  unsigned char s;
  typedef enum {
    INF = 0,
    NaN = 1,
    SIGN = 2,
    DWRK = 3
  } state_t;
  T m;
  U e;
  const U& uzero() const;
  const SimpleFloat<T,W,bits,U>& zero()   const;
  const SimpleFloat<T,W,bits,U>& one()    const;
  const SimpleFloat<T,W,bits,U>& two()    const;
  const SimpleFloat<T,W,bits,U>& pi()     const;
  const SimpleFloat<T,W,bits,U>& halfpi() const;
  const SimpleFloat<T,W,bits,U>& quatpi() const;
  const SimpleFloat<T,W,bits,U>& twopi()  const;
  const SimpleFloat<T,W,bits,U>& sqrt2()  const;
private:
  template <typename V> inline U normalize(V& src) const;
  inline SimpleFloat<T,W,bits,U>& ensureFlag();
  inline unsigned char safeAdd(U& dst, const U& src);
  inline char residue2() const;

  // XXX: these are NOT threadsafe on first call.
  const vector<SimpleFloat<T,W,bits,U> >& exparray()    const;
  const vector<SimpleFloat<T,W,bits,U> >& invexparray() const;
};

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::SimpleFloat() {
  assert(0 < bits && ! (bits & 1));
  s |= (1 << NaN) | (1 << INF);
}

template <typename T, typename W, int bits, typename U> template <typename V> inline SimpleFloat<T,W,bits,U>::SimpleFloat(const V& src) {
  const static V vzero(0);
  s ^= s;
  m  = T(src < vzero ? - src : src);
  e ^= e;
  s |= safeAdd(e, normalize(m));
  if(src < vzero)
    s |= 1 << SIGN;
  ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::SimpleFloat(const SimpleFloat<T,W,bits,U>& src) {
  *this = src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::SimpleFloat(SimpleFloat<T,W,bits,U>&& src) {
  *this = src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::~SimpleFloat() {
  ;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator -  () const {
  auto work(*this);
  work.s ^= 1 << SIGN;
  return work;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator +  (const SimpleFloat<T,W,bits,U>& src) const {
  auto work(*this);
  return work += src;
}

template <typename T, typename W, int bits, typename U>        SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator += (const SimpleFloat<T,W,bits,U>& src) {
  if((s |= src.s & (1 << NaN)) & (1 << NaN))
    return *this;
  if(s & (1 << INF)) {
    if((src.s & (1 << INF)) && (s ^ src.s) & (1 << SIGN))
      s |= 1 << NaN;
    return *this;
  }
  if(src.s & (1 << INF))
    return *this = src;
  if(! m)
    return *this = src;
  if(! src.m)
    return *this;
  if(! ((s ^ src.s) & (1 << SIGN))) {
    if(e >= src.e) {
      m >>= 1;
      s |= safeAdd(e, 1);
      U se(e);
      if(! safeAdd(se, - src.e) && se < U(bits))
        m += src.m >> int(se);
    } else
      return *this = src + *this;
  } else {
    if(e > src.e) {
      U se(e);
      if(! safeAdd(se, - src.e) && se < U(bits))
        m -= src.m >> int(se);
    } else if(e == src.e) {
      if(m >= src.m)
        m -= src.m;
      else
        return *this = src + *this;
    } else
      return *this = src + *this;
  }
  s |= safeAdd(e, normalize(m));
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator -  (const SimpleFloat<T,W,bits,U>& src) const {
  auto work(*this);
  return work -= src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator -= (const SimpleFloat<T,W,bits,U>& src) {
  s ^= 1 << SIGN;
  *this += src;
  s ^= 1 << SIGN;
  return *this;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator *  (const SimpleFloat<T,W,bits,U>& src) const {
  auto work(*this);
  return work *= src;
}

template <typename T, typename W, int bits, typename U>        SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator *= (const SimpleFloat<T,W,bits,U>& src) {
  s ^= src.s & (1 << SIGN);
  if((s |= src.s & (1 << NaN)) & (1 << NaN))
    return *this;
  if((! m) || (! src.m)) {
    s |= 1 << DWRK;
    return ensureFlag();
  }
  if((s |= src.s & (1 << INF)) & (1 << INF))
    return *this;
  auto mm(W(m) * W(src.m));
  s |= safeAdd(e, src.e);
  s |= safeAdd(e, normalize(mm));
  s |= safeAdd(e, U(bits));
  m  = T(mm >> bits);
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline char SimpleFloat<T,W,bits,U>::residue2() const {
  if(uzero() < e || U(bits) <= - e)
    return 0;
  if(! e)
    return char(int(m) & 1);
  return char(int(m >> - int(e)) & 1);
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::operator /  (const SimpleFloat<T,W,bits,U>& src) const {
  auto work(*this);
  return work /= src;
}

template <typename T, typename W, int bits, typename U>        SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator /= (const SimpleFloat<T,W,bits,U>& src) {
  s ^= src.s & (1 << SIGN);
  if((s |= src.s & (1 << NaN)) & (1 << NaN))
    return *this;
  if(! (s & (1 << INF)) && (src.s & (1 << INF))) {
    s |= 1 << DWRK;
    return ensureFlag();
  }
  if(s & (1 << INF)) {
    if(src.s & (1 << INF))
      s |= 1 << NaN;
    return *this;
  }
  if(! src.m) {
    throw "Zero division";
    s |= 1 << NaN;
    return *this;
  }
  if(! m)
    return *this;
  auto mm((W(m) << bits) / W(src.m));
  s |= safeAdd(e, - src.e);
  s |= safeAdd(e, normalize(mm));
  m  = T(mm >> bits);
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>  SimpleFloat<T,W,bits,U>::operator %  (const SimpleFloat<T,W,bits,U>& src) const {
  return *this - (*this / src).floor() * src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator %= (const SimpleFloat<T,W,bits,U>& src) {
  return *this = *this % src;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>  SimpleFloat<T,W,bits,U>::operator <<  (const U& b) const {
  auto work(*this);
  return work <<= b;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator <<= (const U& b) {
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  s |= safeAdd(e, b);
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>  SimpleFloat<T,W,bits,U>::operator >>  (const U& b) const {
  auto work(*this);
  return work >>= b;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator >>= (const U& b) {
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  s |= safeAdd(e, - b);
  return ensureFlag();
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator =  (const SimpleFloat<T,W,bits,U>& src) {
  s = src.s;
  e = src.e;
  m = src.m;
  return *this;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::operator =  (SimpleFloat<T,W,bits,U>&& src) {
  s = move(src.s);
  e = move(src.e);
  m = move(src.m);
  return *this;
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator == (const SimpleFloat<T,W,bits,U>& src) const {
  return ! (*this != src);
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator != (const SimpleFloat<T,W,bits,U>& src) const {
  return (((s | src.s) & ((1 << INF) | (1 << NaN))) ||
           (s != src.s || e != src.e || m != src.m)) &&
         ! (! m && ! src.m);
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator <  (const SimpleFloat<T,W,bits,U>& src) const {
  if((s | src.s) & (1 << NaN))
    throw "compair NaN";
  const auto s_is_minus(s & (1 << SIGN));
  if(s_is_minus ^ (src.s & (1 << SIGN)))
    return s_is_minus;
  if(s & (1 << INF)) {
    if(src.s & (1 << INF))
      throw "compair INF";
    return s_is_minus;
  }
  if(src.s & (1 << INF))
    return ! s_is_minus;
  if(m && src.m) {
    if(e < src.e)
      return ! s_is_minus;
    if(e == src.e)
      return s_is_minus ? src.m < m : m < src.m;
    return s_is_minus;
  }
  return !m ? (bool(src.m) && ! s_is_minus) : s_is_minus;
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator <= (const SimpleFloat<T,W,bits,U>& src) const {
  return *this < src || *this == src;
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator >  (const SimpleFloat<T,W,bits,U>& src) const {
  return ! (*this <= src);
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator >= (const SimpleFloat<T,W,bits,U>& src) const {
  return ! (*this < src);
}

template <typename T, typename W, int bits, typename U> inline bool             SimpleFloat<T,W,bits,U>::operator !  () const {
  return ! m && isfinite(*this);
}

template <typename T, typename W, int bits, typename U> inline                  SimpleFloat<T,W,bits,U>::operator bool () const {
  return ! (!*this);
}

template <typename T, typename W, int bits, typename U> inline                  SimpleFloat<T,W,bits,U>::operator int  () const {
  return int(this->operator T());
}

template <typename T, typename W, int bits, typename U> inline                  SimpleFloat<T,W,bits,U>::operator T    () const {
  auto deci(*this);
  if(deci.s & (1 << INF))
    throw "Inf to convert int";
  if(deci.s & (1 << NaN))
    throw "NaN to convert int";
  if(! deci.m)
    return T(0);
  if(U(bits) <= deci.e || (uzero() < deci.e && (deci.m << int(deci.e)) >> int(deci.e) != deci.m))
    throw "Overflow to convert int.";
  if(deci.e <= - U(bits))
    return T(0);
  if(deci.e <  uzero())
    deci.m >>= - int(deci.e);
  else if(uzero() < deci.e)
    deci.m <<=   int(deci.e);
  return s & (1 << SIGN) ? - deci.m : deci.m;
}

template <typename T, typename W, int bits, typename U> inline                  SimpleFloat<T,W,bits,U>::operator SimpleFloat<T,W,bits,U> () const {
  return *this;
}

template <typename T, typename W, int bits, typename U> template <typename V> inline U SimpleFloat<T,W,bits,U>::normalize(V& src) const {
  V   bt(1);
  int b(0);
  int tb(0);
  for( ; bt; tb ++) {
    if(src & bt)
      b = tb;
    bt <<= 1;
  }
  const auto shift(tb - b - 1);
  assert(0 <= shift);
  src <<= shift;
  return - U(shift);
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::ensureFlag() {
  if(s & (1 << INF))
    s &= ~ (1 << DWRK);
  else if(! m || (s & (1 << DWRK))) {
    e ^= e;
    m ^= m;
    s &= ~ ((1 << DWRK) | (1 << INF));
  }
  return * this;
}

template <typename T, typename W, int bits, typename U> inline unsigned char SimpleFloat<T,W,bits,U>::safeAdd(U& dst, const U& src) {
  const auto dst0(dst);
  dst += src;
  if((dst0 > uzero() && src > uzero() && dst <= uzero()) ||
     (dst0 < uzero() && src < uzero() && dst >= uzero()))
    return 1 << (dst0 < uzero() ? DWRK : INF);
  return 0;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::floor() const {
  if(uzero() <= e)
    return *this;
  if(e <= - U(bits))
    return zero();
  auto deci(*this);
  deci.m >>= - int(deci.e);
  deci.m <<= - int(deci.e);
  return deci;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::ceil() const {
  const auto fl(this->floor());
  if(*this - fl) {
    auto pmone(one());
    pmone.s |= s & (1 << SIGN);
    return fl + pmone;
  }
  return fl;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::abs() const {
  auto work(*this);
  work.s &= ~ (1 << SIGN);
  return work;
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::log() const {
  const static auto einv(one() / one().exp());
  const static auto one_einv(one() + einv);
  if((s & (1 << SIGN)) && m)
    throw "Negative log";
  if(s & ((1 << INF) | (1 << NaN)))
    return *this;
  if(! m) {
    auto work(*this);
    work.s |= (1 << INF) | (1 << SIGN);
    return work;
  }
  if(einv <= *this && *this <= one_einv) {
    // ln(x) = (x - 1) - (x - 1)^2/2 + (x-1)^3/3- ...
    const auto dx(*this - one());
          auto x(dx);
          auto before(one());
          auto res(zero());
    for(int t = 1; (res - before).m; t ++, x *= dx) {
      const auto abst(x / SimpleFloat<T,W,bits,U>(t));
      before = res;
      res   += (t % 2 ? abst : - abst);
    }
    return res;
  }
  static const auto& ea(exparray());
  static const auto& iea(invexparray());
        auto  result(zero());
        auto  work(*this);
  if(one_einv < work) {
    for(int i = min(ea.size(), iea.size()) - 1; 0 < i; i --)
      if(ea[i] <= work) {
        result += one() << U(i - 1);
        work   *= iea[i];
      }
    if(! (work <= one_einv)) {
      result += one();
      work   *= iea[1];
    }
  } else {
    for(int i = min(ea.size(), iea.size()) - 1; 0 < i; i --)
      if(work <= iea[i]) {
        result -= one() << U(i - 1);
        work   *= ea[i];
      }
    if(! (einv <= work)) {
      result -= one();
      work   *= ea[1];
    }
  }
  assert(einv <= work && work <= one_einv);
  return result += work.log();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::exp() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    if(! (s & (1 << NaN)) && (s & (1 << SIGN)))
      return zero();
    return *this;
  }
  if(this->abs() <= one()) {
    // exp(x) = 1 + x/1! + x^2/2! + ...
    auto denom(one());
    auto x(*this);
    auto before(zero());
    auto res(one());
    for(int t = 1; (res - before).m; t ++, x *= *this) {
      before = res;
      denom *= SimpleFloat<T,W,bits,U>(t);
      res   += x / denom;
    }
    return res;
  }
  static const auto& en(exparray());
  static const auto& ien(invexparray());
        auto  work(this->abs());
        auto  result(one());
  for(int i = 1; 0 <= i && i < min(en.size(), ien.size()) && work.floor(); i ++, work >>= U(1))
    if(work.residue2())
      result *= s & (1 << SIGN) ? ien[i] : en[i];
  if(work.floor()) {
    work.s |= 1 << INF;
    return work;
  }
  const auto residue(*this - this->floor());
  assert(residue.abs() <= one());
  return result *= residue.exp();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::sin() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    auto res(*this);
    res.s |= 1 << NaN;
    return res;
  }
  if(- one() <= *this && *this <= one()) {
    // sin(x) = x - x^3/3! + x^5/5! - ...
    const auto sqx(*this * *this);
          auto denom(one());
          auto x(sqx * *this);
          auto before(zero());
          auto res(*this);
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      SimpleFloat<T,W,bits,U> tt(t);
      tt   <<= U(1);
      before = res;
      denom *= - tt * (tt + one());
      res   += x / denom;
    }
    return res;
  }
  if(- halfpi() <= *this && *this <= halfpi())
    return ((*this - quatpi()).cos() + (*this - quatpi()).sin()) / sqrt2();
  if(this->abs() == pi())
    return zero();
  if(- pi() <= *this && *this <= pi())
    return (halfpi() - *this).cos();
  if(- twopi() <= *this && *this <= twopi())
    return - (*this + pi()).sin();
  return (*this % twopi()).sin();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::cos() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    auto res(*this);
    res.s |= 1 << NaN;
    return res;
  }
  if(- one() <= *this && *this <= one()) {
    // cos(x) = 1 - x^2/2! + x^4/4! - ...
    const auto sqx(*this * *this);
          auto denom(one());
          auto x(sqx);
          auto before(zero());
          auto res(one());
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      SimpleFloat<T,W,bits,U> tt(t);
      tt   <<= U(1);
      before = res;
      denom *= - tt * (tt - one());
      res   += x / denom;
    }
    return res;
  }
  if(this->abs() == halfpi())
    return zero();
  if(- halfpi() <= *this && *this <= halfpi())
    return ((*this - quatpi()).cos() - (*this - quatpi()).sin()) / sqrt2();
  if(this->abs() == pi())
    return - one();
  if(- pi() <= *this && *this <= pi())
    return (halfpi() - *this).sin();
  if(- twopi() <= *this && *this <= twopi())
    return - (*this + pi()).cos();
  return (*this % twopi()).cos();
}

template <typename T, typename W, int bits, typename U> SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::atan() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    if(! (s & (1 << NaN)))
      return s & (1 << SIGN) ? - halfpi() : halfpi();
    return *this;
  }
  static const auto half(one() >> U(1));
  static const auto four(one() << U(2));
  static const auto five((one() << U(2)) + one());
  if(- half <= *this && *this <= half) {
    // arctan(x) = x - x^3/3 + x^5/5 - ...
    const auto sqx(*this * *this);
          auto x(sqx * *this);
          auto before(zero());
          auto res(*this);
    for(int t = 1; (res - before).m; t ++, x *= sqx) {
      const auto abst(x / ((SimpleFloat<T,W,bits,U>(t) << U(1)) + one()));
      before = res;
      res   += (t % 2 ? - abst : abst);
    }
    return res;
  }
  // N.B.
  //  atan(u) + atan(v) = atan((u + v) / (1 - uv)) mod pi, uv != 1.
  //    in u = 0.5, v = x - 0.5 case,
  //  atan(x / (1 - x / 2 + 1 / 4)) = atan(.5) + atan(x - .5) =
  //  atan(x / (1.25 - .5 * x)) 
  //  y := x / (1.25 - .5 * x) then,
  //  (1.25 - .5 * x) * y = x,
  //  (5 - 2x) * y = 4 x
  //  x = 5y / (4 + 2y),
  //     y - x = ((4 + 2y) * y - 5y) / (4 + 2y)
  //           = y * (2y - 1) / (4 + 2y)
  //     so 0 <= y and 0 < y case, this makes decreasing function.
  //       (v = x - .5 and 0 <= 2y - 1)
  if(- two() <= *this && *this <= two()) {
    static const auto atanhalf(half.atan());
    if(s & (1 << SIGN))
      return - (- *this).atan();
    const auto v(five * *this / (four + (*this << U(1))) - half);
    assert(v < *this);
    return atanhalf + v.atan();
  }
  // N.B.
  //    in u = v case,
  //  2 atan(u) = atan(2 * u / (1 - u * u))
  //    in u := x + 1 case,
  //  2 atan(1 + x) = atan(2 * (1 + x) / (x + x * x))
  //                = atan(2 / x)
  //    in Y := 2 / x case,
  //  atan(Y) = 2 atan(1 + 2 / Y)
  const auto y(one() + two() / (*this));
  assert(- two() <= y && y <= two());
  return y.atan() << U(1);
}

template <typename T, typename W, int bits, typename U> const vector<SimpleFloat<T,W,bits,U> >& SimpleFloat<T,W,bits,U>::exparray() const {
  static vector<SimpleFloat<T,W,bits,U> > ebuf;
  if(ebuf.size())
    return ebuf;
  ebuf.emplace_back(one());
  ebuf.emplace_back(ebuf[0].exp());
  for(int i = 1; 0 <= i; i ++) {
    const auto en(ebuf[i] * ebuf[i]);
    if(en && isfinite(en))
      ebuf.emplace_back(en);
    else
      break;
  }
  return ebuf;
}

template <typename T, typename W, int bits, typename U> const vector<SimpleFloat<T,W,bits,U> >& SimpleFloat<T,W,bits,U>::invexparray() const {
  static vector<SimpleFloat<T,W,bits,U> > iebuf;
  if(iebuf.size())
    return iebuf;
  const auto& ea(exparray());
  for(int i = 0; 0 <= i && i < ea.size(); i ++) {
    const auto ien(one() / ea[i]);
    if(ien && isfinite(ien))
      iebuf.emplace_back(ien);
    else
      break;
  }
  return iebuf;
}

template <typename T, typename W, int bits, typename U> const U& SimpleFloat<T,W,bits,U>::uzero() const {
  const static U vuzero(0);
  return vuzero;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::zero() const {
  const static SimpleFloat<T,W,bits,U> vzero(0);
  return vzero;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::one() const {
  const static SimpleFloat<T,W,bits,U> vone(1);
  return vone;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::two() const {
  const static auto vtwo(one() << U(1));
  return vtwo;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::pi() const {
  const static auto vpi(quatpi() << U(2));
  return vpi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::halfpi() const {
  const static auto vhalfpi(quatpi() << U(1));
  return vhalfpi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::quatpi() const {
  const static auto vquatpi(one().atan());
  return vquatpi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::twopi() const {
  const static auto vtwopi(quatpi() << U(3));
  return vtwopi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>& SimpleFloat<T,W,bits,U>::sqrt2() const {
  const static auto vsqrt2((one() << U(1)).sqrt());
  return vsqrt2;
}

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U> SimpleFloat<T,W,bits,U>::sqrt() const {
  if(s & ((1 << INF) | (1 << NaN))) {
    auto res(*this);
    if(s & (1 << SIGN)) res.s |= 1 << NaN;
    return res;
  }
  auto res((this->log() >> U(1)).exp());
  // get better accuracy (is this enough?, double accuracy on one loop.)
  // newton's method: 0 == f'(x_n) (x_{n+1} - x_n) + f(x_n)
  //            x_{n+1} := x_n - f(x_n)/f'(x_n).
  //         where f(x) := x_n * x_n - *this
  if(! res)
    return res;
  return (res + *this / res) >> U(1);
}

template <typename T, typename W, int bits, typename U> std::ostream& operator << (std::ostream& os, const SimpleFloat<T,W,bits,U>& v) {
  if(isnan(v))
    return os << "NaN ";
  if(isinf(v))
    return os << (const char*)(v.s & (1 << v.SIGN) ? "-" : "") << "Inf ";
  return os << (const char*)(v.s & (1 << v.SIGN) ? "-" : "") << T(v.m) << "*2^" << v.e << " ";
}

template <typename T, typename W, int bits, typename U> std::istream& operator >> (std::istream& is, SimpleFloat<T,W,bits,U>& v) {
  const static SimpleFloat<T,W,bits,U> two(2);
  const static SimpleFloat<T,W,bits,U> ten(10);
               SimpleFloat<T,W,bits,U> e(0);
  bool mode(false);
  bool sign(false);
  bool fsign(false);
  v = SimpleFloat<T,W,bits,U>(0);
  // skip white spaces.
  while(! is.eof()) {
    const auto buf(is.get());
    if(buf != ' ' && buf != '\t' && buf != '\n') {
      is.unget();
      break;
    }
  }
  while(! is.eof() && ! is.bad() ) {
    const auto buf(is.get());
    switch(buf) {
    case '-':
      sign  = true;
    case '+':
      if(fsign)
        throw "Wrong input";
      fsign = true;
      break;
    case '*':
      if(mode)
        goto ensure;
      if(sign)
        v   = - v;
      mode  = true;
      sign  = false;
      fsign = false;
      if(is.get() != '2') {
        is.unget();
        goto ensure;
      }
      if(is.get() != '^') {
        is.unget();
        is.unget();
        goto ensure;
      }
      break;
    case '.':
      throw "not implemented now";
      break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      if(mode) {
        e *= ten;
        e += SimpleFloat<T,W,bits,U>(int(buf - '0'));
      } else {
        v *= ten;
        v += SimpleFloat<T,W,bits,U>(int(buf - '0'));
      }
      fsign = true;
      break;
    default:
      goto ensure;
    }
  }
 ensure:
  if(sign) {
    if(mode)
      e = - e;
    else
      v = - v;
  }
  v *= pow(two, e);
  return is;
}

template <typename T, typename W, int bits, typename U> static inline bool isinf(const SimpleFloat<T,W,bits,U>& src) {
  return src.s & (1 << src.INF);
}

template <typename T, typename W, int bits, typename U> static inline bool isnan(const SimpleFloat<T,W,bits,U>& src) {
  return src.s & (1 << src.NaN);
}

template <typename T, typename W, int bits, typename U> static inline bool isfinite(const SimpleFloat<T,W,bits,U>& src) {
  return ! (src.s & ((1 << src.INF) | (1 << src.NaN)));
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> floor(const SimpleFloat<T,W,bits,U>& src) {
  return src.floor();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> ceil(const SimpleFloat<T,W,bits,U>& src) {
  return src.ceil();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> abs(const SimpleFloat<T,W,bits,U>& src) {
  return src.abs();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> sqrt(const SimpleFloat<T,W,bits,U>& src) {
  return src.sqrt();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> exp(const SimpleFloat<T,W,bits,U>& src) {
  return src.exp();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> log(const SimpleFloat<T,W,bits,U>& src) {
  return src.log();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> sin(const SimpleFloat<T,W,bits,U>& src) {
  return src.sin();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> cos(const SimpleFloat<T,W,bits,U>& src) {
  return src.cos();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> tan(const SimpleFloat<T,W,bits,U>& src) {
  return src.sin() / src.cos();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> atan(const SimpleFloat<T,W,bits,U>& src) {
  return src.atan();
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> atan2(const SimpleFloat<T,W,bits,U>& y, const SimpleFloat<T,W,bits,U>& x) {
  auto atan0(y.halfpi());
  if(! x && ! y)
    return x / y;
  else if(isfinite(x)) {
    if(! isfinite(y) )
      goto ensure;
    if(! x)
      goto ensure;
    const auto yoverx((y / x).abs());
    if(! isfinite(yoverx) )
      goto ensure;
    const auto atan00(yoverx.atan());
    if(! isfinite(atan00) )
      goto ensure;
    atan0 = atan00;
    goto ensure;
  } else if(isfinite(y)) {
    atan0 = x.zero();
    goto ensure;
  }
  return y;
 ensure:
  if(y.s & (1 << y.SIGN)) {
    if(x.s & (1 << x.SIGN))
      atan0 = - (x.pi() - atan0);
    else
      atan0 = - atan0;
  } else if(x.s & (1 << x.SIGN))
    atan0  = x.pi() - atan0;
  return atan0;
}

template <typename T, typename W, int bits, typename U> static inline SimpleFloat<T,W,bits,U> pow(const SimpleFloat<T,W,bits,U>& src, const SimpleFloat<T,W,bits,U>& dst) {
  if(! dst) {
    if(! src)
      throw "0^0";
    return dst.one();
  }
  return exp(log(src) * dst);
}


// class complex part:
template <typename T> class Complex {
public:
  inline Complex();
  inline Complex(const Complex<T>& s);
  inline Complex(Complex<T>&& s);
  inline Complex(const T& real, const T& imag = T(0));
  inline Complex(T&& real);
  inline Complex(T&& real, T&& imag);
  inline ~Complex();

  inline Complex<T>  operator ~  ()                    const;
  inline Complex<T>  operator -  ()                    const;
  inline Complex<T>  operator +  (const Complex<T>& s) const;
  inline Complex<T>& operator += (const Complex<T>& s);
  inline Complex<T>  operator -  (const Complex<T>& s) const;
  inline Complex<T>& operator -= (const Complex<T>& s);
  inline Complex<T>  operator *  (const T& s)          const;
  inline Complex<T>& operator *= (const T& s);
  inline Complex<T>  operator *  (const Complex<T>& s) const;
  inline Complex<T>& operator *= (const Complex<T>& s);
  inline Complex<T>  operator /  (const T& s)          const;
  inline Complex<T>& operator /= (const T& s);
  inline Complex<T>  operator /  (const Complex<T>& s) const;
  inline Complex<T>& operator /= (const Complex<T>& s);
  inline bool        operator == (const Complex<T>& s) const;
  inline bool        operator != (const Complex<T>& s) const;
  inline bool        operator !  ()                    const;
  inline Complex<T>  operator &  (const Complex<T>& s) const;
  inline Complex<T>& operator &= (const Complex<T>& s);
  inline Complex<T>  operator |  (const Complex<T>& s) const;
  inline Complex<T>& operator |= (const Complex<T>& s);
  inline Complex<T>  operator ^  (const Complex<T>& s) const;
  inline Complex<T>& operator ^= (const Complex<T>& s);
  inline bool        operator && (const Complex<T>& s) const;
  inline bool        operator || (const Complex<T>& s) const;
  inline Complex<T>& operator =  (const Complex<T>& s);
  inline Complex<T>& operator =  (Complex<T>&& s);
  inline T&          operator [] (const size_t& i);
  inline             operator bool () const;
  inline             operator T    () const;
  
  const Complex<T>& i() const;

  // friend std::ostream& operator << (std::ostream& os, const SimpleVector<T>& v);
  // friend std::istream& operator >> (std::istream& os, SimpleVector<T>& v);
  
  inline T  abs() const;
  inline T  arg() const;
  inline T& real();
  inline T& imag();
  inline const T& real() const;
  inline const T& imag() const;
  T _real;
  T _imag;
};

template <typename T> inline Complex<T>::Complex() {
  ;
}

template <typename T> inline Complex<T>::Complex(const Complex<T>& src) {
  *this = src;
}

template <typename T> inline Complex<T>::Complex(Complex<T>&& src) {
  *this = src;
}

template <typename T> inline Complex<T>::Complex(const T& real, const T& imag) {
  _real = real;
  _imag = imag;
  return;
}

template <typename T> inline Complex<T>::Complex(T&& real) {
  const static T zero(0);
  _real = move(real);
  _imag = zero;
  return;
}

template <typename T> inline Complex<T>::Complex(T&& real, T&& imag) {
  _real = move(real);
  _imag = move(imag);
  return;
}

template <typename T> inline Complex<T>::~Complex() {
  ;
}

template <typename T> inline Complex<T> Complex<T>::operator ~ () const {
  return Complex<T>(  _real, - _imag);
}

template <typename T> inline Complex<T> Complex<T>::operator - () const {
  return Complex<T>(- _real, - _imag);
}

template <typename T> inline Complex<T> Complex<T>::operator + (const Complex<T>& s) const {
  auto result(*this);
  return result += s;
}

template <typename T> inline Complex<T>& Complex<T>::operator += (const Complex<T>& s) {
  _real += s._real;
  _imag += s._imag;
  return *this;
}

template <typename T> inline Complex<T> Complex<T>::operator - (const Complex<T>& s) const {
  auto result(*this);
  return result -= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator -= (const Complex<T>& s) {
  _real -= s._real;
  _imag -= s._imag;
  return *this;
}

template <typename T> inline Complex<T>  Complex<T>::operator * (const T& s) const {
  auto result(*this);
  return result *= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator *= (const T& s) {
  _real *= s;
  _imag *= s;
  return *this;
}
  
template <typename T> inline Complex<T> Complex<T>::operator * (const Complex<T>& s) const {
  return Complex<T>(_real * s._real - _imag * s._imag,
                    _real * s._imag + _imag * s._real);
}
 
template <typename T> inline Complex<T>& Complex<T>::operator *= (const Complex<T>& s) {
  return (*this) = (*this) * s;
}

template <typename T> inline Complex<T> Complex<T>::operator / (const T& s) const {
  auto result(*this);
  return result /= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator /= (const T& s) {
  _real /= s;
  _imag /= s;
  return *this;
}

template <typename T> inline Complex<T> Complex<T>::operator / (const Complex<T>& s) const {
  return (*this * (~ s)) / (s._real * s._real + s._imag * s._imag);
}

template <typename T> inline Complex<T>& Complex<T>::operator /= (const Complex<T>& s) {
  return *this = *this / s;
}

template <typename T> inline bool Complex<T>::operator == (const Complex<T>& s) const {
  return !(*this != s);
}

template <typename T> inline bool Complex<T>::operator != (const Complex<T>& s) const {
  return (_real != s._real) || (_imag != s._imag);
}

template <typename T> inline bool Complex<T>::operator ! () const {
  return !_real && !_imag;
}

template <typename T> inline Complex<T> Complex<T>::operator & (const Complex<T>& s) const {
  auto result(*this);
  return result &= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator &= (const Complex<T>& s) {
  _real &= s._real;
  _imag &= s._imag;
  return *this;
}

template <typename T> inline Complex<T> Complex<T>::operator | (const Complex<T>& s) const {
  auto result(*this);
  return result |= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator |= (const Complex<T>& s) {
  _real |= s._real;
  _imag |= s._imag;
  return *this;
}

template <typename T> inline Complex<T> Complex<T>::operator ^ (const Complex<T>& s) const {
  auto result(*this);
  return result ^= s;
}

template <typename T> inline Complex<T>& Complex<T>::operator ^= (const Complex<T>& s) {
  _real ^= s._real;
  _imag ^= s._imag;
  return *this;
}

template <typename T> inline bool Complex<T>::operator && (const Complex<T>& s) const {
  return *this && s;
}

template <typename T> inline bool Complex<T>::operator || (const Complex<T>& s) const {
  return *this || s;
}

template <typename T> inline Complex<T>& Complex<T>::operator =  (const Complex<T>& s) {
  _real = s._real;
  _imag = s._imag;
  return *this;
}

template <typename T> inline Complex<T>& Complex<T>::operator =  (Complex<T>&& s) {
  _real = move(s._real);
  _imag = move(s._imag);
  return *this;
}

template <typename T> inline T& Complex<T>::operator [] (const size_t& i) {
  assert(0 <= i && i < 2);
  if(i)
    return _imag;
  return _real;
}

template <typename T> inline Complex<T>::operator bool () const {
  return ! (! *this);
}

template <typename T> inline Complex<T>::operator T () const {
  return this->_real;
}

template <typename T> const Complex<T>& Complex<T>::i() const {
  const static auto I(Complex<T>(T(0), T(1)));
  return I;
}

template <typename T> inline T Complex<T>::abs() const {
  return sqrt(_real * _real + _imag * _imag);
}

template <typename T> inline T Complex<T>::arg() const {
  return atan2(_imag, _real);
}

template <typename T> inline T& Complex<T>::real() {
  return _real;
}

template <typename T> inline T& Complex<T>::imag() {
  return _imag;
}

template <typename T> inline const T& Complex<T>::real() const {
  return _real;
}

template <typename T> inline const T& Complex<T>::imag() const {
  return _imag;
}

template <typename T> std::ostream& operator << (std::ostream& os, const Complex<T>& v) {
  return os << v.real() << "+i" << v.imag();
}

template <typename T> std::istream& operator >> (std::istream& is, Complex<T>& v) {
  is >> v._real;
  if('+' != is.get()) {
    is.unget();
    goto ensure;
  }
  if('i' != is.get()) {
    is.unget();
    is.unget();
    goto ensure;
  }
  is >> v._imag;
  return is;
 ensure:
  v._imag = T(0);
  return is;
}


template <typename T> static inline T abs(const Complex<T>& s) {
  return s.abs();
}

template <typename T> static inline T arg(const Complex<T>& s) {
  return s.arg();
}

template <typename T> static inline const T& real(const Complex<T>& s) {
  return s.real();
}

template <typename T> static inline const T& imag(const Complex<T>& s) {
  return s.imag();
}

template <typename T> static inline Complex<T> exp(const Complex<T>& s) {
  return Complex<T>(exp(s.real())) * Complex<T>(cos(s.imag()), sin(s.imag()));
}

template <typename T> static inline Complex<T> log(const Complex<T>& s) {
  // N.B. main branch
  return Complex<T>(log(abs(s)), arg(s));
}

template <typename T> static inline Complex<T> sqrt(const Complex<T>& s) {
  return exp(log(s) * Complex<T>(T(1) / T(2)));
}

template <typename T> static inline Complex<T> csin(const Complex<T>& s) {
  return (exp(Complex<T>(T(0), s)) - exp(Complex<T>(T(0), - s))) / Complex<T>(T(0), T(2));
}

template <typename T> static inline Complex<T> ccos(const Complex<T>& s) {
  return (exp(Complex<T>(T(0), s)) + exp(Complex<T>(T(0), - s))) / T(2);
}

template <typename T> static inline Complex<T> ctan(const Complex<T>& s) {
  return csin(s) / ccos(s);
}

template <typename T> static inline Complex<T> ccsc(const Complex<T>& s) {
  return Complex<T>(T(1)) / csin(s);
}

template <typename T> static inline Complex<T> csec(const Complex<T>& s) {
  return Complex<T>(T(1)) / ccos(s);
}

template <typename T> static inline T ccot(const T& s) {
  return Complex<T>(T(1)) / ctan(s);
}

template <typename T> using complex = Complex<T>;

# if _FLOAT_BITS_ == 8
  typedef uint8_t myuint;
  typedef int8_t  myint;
  typedef SimpleFloat<myuint, uint16_t, 8, int64_t> myfloat;
# elif _FLOAT_BITS_ == 16
  typedef uint16_t myuint;
  typedef int16_t  myint;
  typedef SimpleFloat<myuint, uint32_t, 16, int64_t> myfloat;
# elif _FLOAT_BITS_ == 32
  typedef uint32_t myuint;
  typedef int32_t  myint;
  typedef SimpleFloat<myuint, uint64_t, 32, int64_t> myfloat;
# elif _FLOAT_BITS_ == 64
  typedef uint64_t myuint;
  typedef int64_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 64>, 64, int64_t> myfloat;
# elif _FLOAT_BITS_ == 128
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef Signed<uint128_t, 128> int128_t;
  typedef uint128_t myuint;
  typedef int128_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 128>, 128, int64_t> myfloat;
# elif _FLOAT_BITS_ == 256
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef DUInt<uint128_t, 128> uint256_t;
  typedef Signed<uint256_t, 256> int256_t;
  typedef uint256_t myuint;
  typedef int256_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 256>, 256, int64_t> myfloat;
# elif _FLOAT_BITS_ == 512
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef DUInt<uint128_t, 128> uint256_t;
  typedef DUInt<uint256_t, 256> uint512_t;
  typedef Signed<uint512_t, 512> int512_t;
  typedef uint512_t myuint;
  typedef int512_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 512>, 512, int64_t> myfloat;
# elif _FLOAT_BITS_ == 1024
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef DUInt<uint128_t, 128> uint256_t;
  typedef DUInt<uint256_t, 256> uint512_t;
  typedef DUInt<uint512_t, 512> uint1024_t;
  typedef Signed<uint1024_t, 1024> int1024_t;
  typedef uint1024_t myuint;
  typedef int1024_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 1024>, 1024, int64_t> myfloat;
# else
#   error cannot handle float
# endif
#endif


// N.B. start simplelin.
template <typename T> class SimpleVector {
public:
  inline SimpleVector();
  inline SimpleVector(const int& size);
  inline SimpleVector(const SimpleVector<T>& other);
  inline SimpleVector(SimpleVector<T>&& other);
  inline ~SimpleVector();
  
  inline       SimpleVector<T>  operator -  () const;
  inline       SimpleVector<T>  operator +  (const SimpleVector<T>& other) const;
  inline       SimpleVector<T>& operator += (const SimpleVector<T>& other);
  inline       SimpleVector<T>  operator -  (const SimpleVector<T>& other) const;
  inline       SimpleVector<T>& operator -= (const SimpleVector<T>& other);
  inline       SimpleVector<T>  operator *  (const T& other) const;
  inline       SimpleVector<T>& operator *= (const T& other);
  inline       SimpleVector<T>  operator /  (const T& other) const;
  inline       SimpleVector<T>& operator /= (const T& other);
  inline       SimpleVector<T>& operator =  (const SimpleVector<T>& other);
  inline       SimpleVector<T>& operator =  (SimpleVector<T>&& other);
  inline       bool             operator == (const SimpleVector<T>& other) const;
  inline       bool             operator != (const SimpleVector<T>& other) const;
  inline       T                dot         (const SimpleVector<T>& other) const;
  inline       T&               operator [] (const int& idx);
  inline const T&               operator [] (const int& idx) const;
  template <typename U> inline SimpleVector<U> real() const;
  template <typename U> inline SimpleVector<U> imag() const;
  template <typename U> inline SimpleVector<U> cast() const;
  inline const int& size() const;
  inline       void resize(const int& size);
  inline       SimpleVector<T>  subVector(const int& i, const int& s) const;
  inline       SimpleVector<T>& setVector(const int& i, const SimpleVector<T>& d);
  inline       SimpleVector<T>& O();
  inline       SimpleVector<T>& I(const T& r = T(1));
  inline       SimpleVector<T>& ek(const int& i, const T& r = T(1));
  
  // friend std::ostream& operator << (std::ostream& os, const SimpleVector<T>& v);
  // friend std::istream& operator >> (std::istream& os, SimpleVector<T>& v);
  
  T*  entity;
  int esize;
};

template <typename T> inline SimpleVector<T>::SimpleVector() {
  entity = NULL;
  esize  = 0;
}

template <typename T> inline SimpleVector<T>::SimpleVector(const int& size) {
  assert(0 <= size);
  if(0 < size)
    this->entity = new T[size];
  else
    this->entity = NULL;
  this->esize  = size;
  return;
}

template <typename T> inline SimpleVector<T>::SimpleVector(const SimpleVector<T>& other) {
  entity = new T[other.esize];
  esize  = other.esize;
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] = other.entity[i];
  return;
}

template <typename T> inline SimpleVector<T>::SimpleVector(SimpleVector<T>&& other) {
  esize  = move(other.esize);
  entity = move(other.entity);
  other.entity = 0;
  return;
}

template <typename T> inline SimpleVector<T>::~SimpleVector() {
  if(entity)
    delete[] entity;
  entity = NULL;
  return;
}

template <typename T> inline SimpleVector<T> SimpleVector<T>::operator - () const {
  SimpleVector<T> res(esize);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    res.entity[i] = - entity[i];
  return res;
}

template <typename T> inline SimpleVector<T> SimpleVector<T>::operator + (const SimpleVector<T>& other) const {
  auto res(*this);
  return res += other;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::operator += (const SimpleVector<T>& other) {
  assert(esize == other.esize);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] += other.entity[i];
  return *this;
}

template <typename T> inline SimpleVector<T> SimpleVector<T>::operator - (const SimpleVector<T>& other) const {
  auto res(*this);
  return res -= other;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::operator -= (const SimpleVector<T>& other) {
  return *this += - other;
}

template <typename T> inline SimpleVector<T> SimpleVector<T>::operator * (const T& other) const {
  auto res(*this);
  return res *= other;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::operator *= (const T& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] *= other;
  return *this;
}

template <typename T> inline SimpleVector<T> SimpleVector<T>::operator / (const T& other) const {
  auto res(*this);
  return res /= other;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::operator = (const SimpleVector<T>& other) {
  if(entity == other.entity && esize == other.esize)
    return *this;
  if(esize != other.esize) {
    if(entity)
      delete[] entity;
    if(other.esize)
      entity = new T[other.esize];
    else
      entity = NULL;
  }
  esize = other.esize;
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] = other.entity[i];
  return *this;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::operator = (SimpleVector<T>&& other) {
  if(entity == other.entity && esize == other.esize)
    return *this;
  esize  = move(other.esize);
  if(entity)
    delete[] entity;
  entity = move(other.entity);
  other.esize  = 0;
  other.entity = NULL;
  return *this;
}

template <typename T> inline bool SimpleVector<T>::operator == (const SimpleVector<T>& other) const {
  return ! (*this != other);
}

template <typename T> inline bool SimpleVector<T>::operator != (const SimpleVector<T>& other) const {
  if(esize != other.esize)
    return true;
  for(int i = 0; i < esize; i ++)
    if(entity[i] != other.entity[i])
      return true;
  return false;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::operator /= (const T& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] /= other;
  return *this;
}

template <typename T> inline T SimpleVector<T>::dot(const SimpleVector<T>& other) const {
  assert(esize == other.esize);
  SimpleVector<T> work(other.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    work[i] = entity[i] * other.entity[i];
  auto res(work[0]);
  for(int i = 1; i < esize; i ++)
    res += work[i];
  return res;
}

template <typename T> inline T& SimpleVector<T>::operator [] (const int& idx) {
  assert(0 <= idx && idx < esize && entity);
  return entity[idx];
}

template <typename T> inline const T& SimpleVector<T>::operator [] (const int& idx) const {
  assert(0 <= idx && idx < esize && entity);
  return entity[idx];
}

template <typename T> template <typename U> inline SimpleVector<U> SimpleVector<T>::real() const {
  SimpleVector<U> result(esize);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    result.entity[i] = U(entity[i].real());
  return result;
}

template <typename T> template <typename U> inline SimpleVector<U> SimpleVector<T>::imag() const {
  SimpleVector<U> result(esize);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    result.entity[i] = U(entity[i].imag());
  return result;
}

template <typename T> template <typename U> inline SimpleVector<U> SimpleVector<T>::cast() const {
  SimpleVector<U> result(esize);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    result.entity[i] = U(entity[i]);
  return result;
}

template <typename T> inline const int& SimpleVector<T>::size() const {
  return esize;
}

template <typename T> inline void SimpleVector<T>::resize(const int& size) {
  assert(0 <= size);
  if(size != esize) {
    esize = size;
    if(entity)
      delete[] entity;
    if(0 < size)
      entity = new T[esize];
    else
      entity = NULL;
  }
  return;
}

template <typename T> inline SimpleVector<T>  SimpleVector<T>::subVector(const int& i, const int& s) const {
  assert(0 <= s && 0 <= i && i + s <= size());
  SimpleVector<T> res(s);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int ii = i; ii < i + s; ii ++)
    res[ii - i] = (*this)[ii];
  return res;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::setVector(const int& i, const SimpleVector<T>& d) {
  assert(0 <= i && i + d.size() <= size());
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int ii = i; ii < i + d.size(); ii ++)
    (*this)[ii] = d[ii - i];
  return *this;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::O() {
  const static T zero(0);
  return I(zero);
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::I(const T& r) {
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < size(); i ++)
    (*this)[i] = r;
  return *this;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::ek(const int& i, const T& r) {
  static const T zero(0);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int ii = 0; ii < size(); ii ++)
    (*this)[ii] = ii == i ? r : zero;
  return *this;
}

template <typename T> std::ostream& operator << (std::ostream& os, const SimpleVector<T>& v) {
  SimpleVector<string> buf(v.size());
  int M(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < v.size(); i ++) {
    stringstream ss;
    ss << v[i];
    buf[i] = ss.str();
#if defined(_OPENMP)
#pragma omp critical
#endif
    {
      M = max(int(buf[i].size()), M);
    }
  }
  os << v.size() << " : [";
  for(int i = 0; i < buf.size(); i ++) {
    for(int j = buf[i].size(); j <= M; j ++)
      os << " ";
    os << buf[i];
    if(i < buf.size() - 1) os << ", ";
  }
  os << "]" << endl;
  return os;
}

template <typename T> std::istream& operator >> (std::istream& is, SimpleVector<T>& v) {
  int s;
  is >> s;
  if(s <= 0) return is;
  v.resize(s);
  int i(0);
  for( ; i < v.size() && ! is.eof() && ! is.bad(); ) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == '[' || c == ',' || c == ']' || c == '\n') continue;
    is.unget();
    is >> v[i ++];
  }
  while(!is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == '\n') continue;
    else if(c == ']') break;
    is.unget();
    cerr << "XXX SimpleVector<T>::operator >> (\']\')" << flush;
    break;
  }
  while(!is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t') continue;
    else if(c == '\n') break;
    is.unget();
    cerr << "XXX SimpleVector<T>::operator >> (\'\\n\')" << flush;
    break;
  }
  if(i < v.size()) {
    cerr << "XXX SimpleVector<T>::operator >> (index)" << flush;
    for( ; i < v.size(); i ++)
      v[i] = T(0);
  }
  return is;
}


template <typename T> class SimpleMatrix {
public:
  inline SimpleMatrix();
  inline SimpleMatrix(const int& rows, const int& cols);
  inline SimpleMatrix(const SimpleMatrix<T>& other);
  inline SimpleMatrix(SimpleMatrix<T>&& other);
  inline ~SimpleMatrix();
  
  inline       SimpleMatrix<T>  operator -  () const;
  inline       SimpleMatrix<T>  operator +  (const SimpleMatrix<T>& other) const;
  inline       SimpleMatrix<T>& operator += (const SimpleMatrix<T>& other);
  inline       SimpleMatrix<T>  operator -  (const SimpleMatrix<T>& other) const;
  inline       SimpleMatrix<T>& operator -= (const SimpleMatrix<T>& other);
  inline       SimpleMatrix<T>  operator *  (const T& other) const;
  inline       SimpleMatrix<T>& operator *= (const T& other);
  inline       SimpleMatrix<T>  operator *  (const SimpleMatrix<T>& other) const;
  inline       SimpleMatrix<T>& operator *= (const SimpleMatrix<T>& other);
  inline       SimpleVector<T>  operator *  (const SimpleVector<T>& other) const;
  inline       SimpleMatrix<T>  operator /  (const T& other) const;
  inline       SimpleMatrix<T>& operator /= (const T& other);
  inline       SimpleMatrix<T>& operator =  (const SimpleMatrix<T>& other);
  inline       SimpleMatrix<T>& operator =  (SimpleMatrix<T>&& other);
  inline       bool             operator == (const SimpleMatrix<T>& other) const;
  inline       bool             operator != (const SimpleMatrix<T>& other) const;
  inline       T&               operator () (const int& y, const int& x);
  inline const T&               operator () (const int& y, const int& x) const;
  inline       SimpleVector<T>& row(const int& y);
  inline const SimpleVector<T>& row(const int& y) const;
  inline const SimpleVector<T>  col(const int& x) const;
  inline       void             setCol(const int& x, const SimpleVector<T>& other);
  // N.B. transpose : exhaust of the resource, so Eigen library handles better.
  inline       SimpleMatrix<T>  transpose() const;
  inline       SimpleMatrix<T>  subMatrix(const int& y, const int& x, const int& h, const int& w) const;
  inline       SimpleMatrix<T>& setMatrix(const int& y, const int& x, const SimpleMatrix<T>& d);
  inline       SimpleMatrix<T>& O();
  inline       SimpleMatrix<T>& I(const T& r = T(1));
  inline       T                determinant() const;
  inline       SimpleMatrix<T>  inverse() const;
  inline       SimpleVector<T>  solve(SimpleVector<T> other) const;
  inline       SimpleVector<T>  projectionPt(const SimpleVector<T>& other) const;
  inline       SimpleMatrix<T>& fillP(const vector<int>& idx);
  inline       SimpleMatrix<T>  QR() const;
  inline       SimpleMatrix<T>  SVD() const;
  inline       pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SVD(const SimpleMatrix<T>& src) const;
  inline       SimpleVector<T>  innerFix(const SimpleMatrix<T>& A, vector<pair<T, int> >& fidx);
  inline       SimpleVector<T>  inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const;
  template <typename U> inline SimpleMatrix<U> real() const;
  template <typename U> inline SimpleMatrix<U> imag() const;
  template <typename U> inline SimpleMatrix<U> cast() const;
  inline const int& rows() const;
  inline const int& cols() const;
  inline       void resize(const int& rows, const int& cols);
  myfloat      epsilon;

  // friend std::ostream& operator << (std::ostream& os, const SimpleVector<T>& v);
  // friend std::istream& operator >> (std::istream& os, SimpleVector<T>& v);

private:
  // this isn't better idea for faster calculations.
  SimpleVector<T>* entity;
  int              erows;
  int              ecols;
};

template <typename T> inline SimpleMatrix<T>::SimpleMatrix() {
  erows  = 0;
  ecols  = 0;
  entity = NULL;
#if defined(_FLOAT_BITS_)
  epsilon = myfloat(1) >> int64_t(_FLOAT_BITS_ - 1);
#else
  epsilon = std::numeric_limits<myfloat>::epsilon();
#endif
  return;
}

template <typename T> inline SimpleMatrix<T>::SimpleMatrix(const int& rows, const int& cols) {
  assert(0 <= rows && 0 <= cols);
  if(0 < rows)
    entity = new SimpleVector<T>[rows];
  else
    entity = NULL;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < rows; i ++)
    entity[i].resize(cols);
  erows = rows;
  ecols = cols;
#if defined(_FLOAT_BITS_)
  epsilon = myfloat(1) >> int64_t(_FLOAT_BITS_ - 1);
#else
  epsilon = std::numeric_limits<myfloat>::epsilon();
#endif
  return; 
}

template <typename T> inline SimpleMatrix<T>::SimpleMatrix(const SimpleMatrix<T>& other) {
  erows = other.erows;
  ecols = other.ecols;
  entity = new SimpleVector<T>[other.erows];
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] = other.entity[i];
  epsilon = other.epsilon;
  return;
}

template <typename T> inline SimpleMatrix<T>::SimpleMatrix(SimpleMatrix<T>&& other) {
  erows  = move(other.erows);
  ecols  = move(other.ecols);
  entity = move(other.entity);
  other.entity = NULL;
  epsilon = move(other.epsilon);
  return;
}

template <typename T> inline SimpleMatrix<T>::~SimpleMatrix() {
  if(entity)
    delete[] entity;
  entity = NULL;
  return;
}
  
template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::operator - () const {
  SimpleMatrix<T> res(erows, ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    res.entity[i] = - entity[i];
  return res;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::operator + (const SimpleMatrix<T>& other) const {
  auto res(*this);
  return res += other;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::operator += (const SimpleMatrix<T>& other) {
  assert(erows == other.erows && ecols == other.ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] += other.entity[i];
  return *this;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::operator - (const SimpleMatrix<T>& other) const {
  auto res(*this);
  return res -= other;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::operator -= (const SimpleMatrix<T>& other) {
  return *this += - other;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::operator * (const T& other) const {
  auto res(*this);
  return res *= other;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::operator *= (const T& other) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] *= other;
  return *this;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::operator * (const SimpleMatrix<T>& other) const {
  assert(ecols == other.erows && entity && other.entity);
  auto            derived(other.transpose());
  SimpleMatrix<T> res(erows, other.ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++) {
          SimpleVector<T>& resi(res.entity[i]);
    const SimpleVector<T>& ei(entity[i]);
    for(int j = 0; j < other.ecols; j ++)
      resi[j] = ei.dot(derived.entity[j]);
  }
  return res;

}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::operator *= (const SimpleMatrix<T>& other) {
  return *this = *this * other;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::operator * (const SimpleVector<T>& other) const {
  assert(ecols == other.size());
  SimpleVector<T> res(erows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    res[i] = entity[i].dot(other);
  return res;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::operator / (const T& other) const {
  auto res(*this);
  return res /= other;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::operator /= (const T& other) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] /= other;
  return *this;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::operator = (const SimpleMatrix<T>& other) {
  if(entity == other.entity && erows == other.erows && ecols == other.ecols)
    return *this;
  if(erows != other.erows || ecols != other.ecols) {
    if(erows != other.erows) {
      if(entity)
        delete[] entity;
      if(other.erows)
        entity = new SimpleVector<T>[other.erows];
      else
        entity = NULL;
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < other.erows; i ++)
      entity[i].resize(other.ecols);
  }
  erows = other.erows;
  ecols = other.ecols;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] = other.entity[i];
  return *this;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::operator = (SimpleMatrix<T>&& other) {
  if(entity == other.entity && erows == other.erows && ecols == other.ecols)
    return *this;
  erows  = move(other.erows);
  ecols  = move(other.ecols);
  if(entity)
    delete[] entity;
  entity = move(other.entity);
  other.erows  = 0;
  other.ecols  = 0;
  other.entity = NULL;
  return *this;
}

template <typename T> inline bool SimpleMatrix<T>::operator == (const SimpleMatrix<T>& other) const {
  return ! (*this != other);
}

template <typename T> inline bool SimpleMatrix<T>::operator != (const SimpleMatrix<T>& other) const {
  if(erows != other.erows || ecols != other.ecols)
    return true;
  for(int i = 0; i < erows; i ++)
    if(entity[i] != other.entity[i])
      return true;
  return false;
}

template <typename T> inline T& SimpleMatrix<T>::operator () (const int& y, const int& x) {
  assert(0 <= y && y < erows && entity);
  return entity[y][x];
}

template <typename T> inline const T& SimpleMatrix<T>::operator () (const int& y, const int& x) const {
  assert(0 <= y && y < erows && entity);
  return entity[y][x];
}

template <typename T> inline SimpleVector<T>& SimpleMatrix<T>::row(const int& y) {
  assert(0 <= y && y < erows && entity);
  return entity[y];
}

template <typename T> inline const SimpleVector<T>& SimpleMatrix<T>::row(const int& y) const {
  assert(0 <= y && y < erows && entity);
  return entity[y];
}

template <typename T> inline const SimpleVector<T> SimpleMatrix<T>::col(const int& x) const {
  assert(0 <= erows && 0 <= x && x < ecols && entity);
  SimpleVector<T> res(erows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    res[i] = entity[i][x];
  return res;
}

template <typename T> inline void SimpleMatrix<T>::setCol(const int& x, const SimpleVector<T>& other) {
  assert(0 <= x && x < ecols && other.size() == erows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    entity[i][x] = other[i];
  return;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::transpose() const {
  SimpleMatrix<T> res(ecols, erows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < ecols; i ++) {
    SimpleVector<T>& resi(res.entity[i]);
    for(int j = 0; j < erows; j ++)
      resi[j] = entity[j][i];
  }
  return res;
}

template <typename T> inline SimpleMatrix<T>  SimpleMatrix<T>::subMatrix(const int& y, const int& x, const int& h, const int& w) const {
  assert(0 <= h && 0 <= w && 0 <= y && y + h <= rows() && 0 <= x && x + w <= cols());
  SimpleMatrix<T> res(h, w);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = y; i < y + h; i ++)
    for(int j = x; j < x + w; j ++)
      res(i - y, j - x) = (*this)(i, j);
  return res;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::setMatrix(const int& y, const int& x, const SimpleMatrix<T>& d) {
  assert(0 <= y && y + d.rows() <= rows() && 0 <= x && x + d.cols() <= cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = y; i < y + d.rows(); i ++)
    for(int j = x; j < x + d.cols(); j ++)
      (*this)(i, j) = d(i - y, j - x);
  return *this;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::O() {
  const static T zero(0);
  return I(zero);
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::I(const T& r) {
  const static T zero(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < rows(); i ++)
    for(int j = 0; j < cols(); j ++)
      (*this)(i, j) = (i == j ? r : zero);
  return *this;
}

template <typename T> inline T SimpleMatrix<T>::determinant() const {
  assert(0 <= erows && 0 <= ecols && erows == ecols);
  T det(1);
  auto work(*this);
  for(int i = 0; i < erows; i ++) {
    int xchg = i;
    for(int j = i + 1; j < erows; j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    swap(work.entity[i], work.entity[xchg]);
    const auto& ei(work.entity[i]);
    const auto& eii(ei[i]);
    if(epsilon * ei.dot(ei) < eii * eii) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = i + 1; j < erows; j ++) {
        const auto ratio(work.entity[j][i] / eii);
        work.entity[j] -= ei * ratio;
      }
    }
    det *= eii;
  }
  return det;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::inverse() const {
  // XXX: extremely slow implementation.
  SimpleMatrix<T> result(erows, ecols);
  result.I();
  for(int i = 0; i < result.cols(); i ++)
    result.setCol(i, solve(result.col(i)));
  return result;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::solve(SimpleVector<T> other) const {
  assert(0 <= erows && 0 <= ecols && erows == ecols && entity && erows == other.size());
  auto work(*this);
  for(int i = 0; i < erows; i ++) {
    int xchg = i;
    for(int j = i + 1; j < erows; j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    swap(work.entity[i], work.entity[xchg]);
    swap(other[i], other[xchg]);
    const auto& ei(work.entity[i]);
    const auto& eii(ei[i]);
    if(epsilon * ei.dot(ei) < eii * eii) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = i + 1; j < erows; j ++) {
        const auto ratio(work.entity[j][i] / eii);
        work.entity[j] -= ei       * ratio;
        other[j]       -= other[i] * ratio;
      }
    }
  }
  for(int i = erows - 1; 0 <= i; i --) {
    if(work.entity[i][i] == T(0)) continue;
    const auto buf(other[i] / work.entity[i][i]);
    if(!isfinite(buf) || isnan(buf)) {
    //  assert(!isfinite(work.entity[i][i] / other[i]) || isnan(work.entity[i][i] / other[i]));
      continue;
    }
    other[i]    = buf;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = i - 1; 0 <= j; j --)
      other[j] -= other[i] * work.entity[j][i];
  }
  return other;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::projectionPt(const SimpleVector<T>& other) const {
  assert(0 < erows && 0 < ecols && ecols == other.size());
  // also needs class or this->transpose() * (*this) == I assertion is needed.
  SimpleMatrix<T> work(erows, ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < work.rows(); i ++)
    work.row(i) = entity[i] * entity[i].dot(other);
  SimpleVector<T> res(ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < other.size(); i ++) {
    res[i] = T(0);
    for(int j = 0; j < erows; j ++)
      res[i] += work(j, i);
  }
  return res;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::fillP(const vector<int>& idx) {
  int ii(0);
  for(int j = 0; j < cols() && ii < idx.size(); j ++) {
    SimpleVector<T> ek(cols());
    ek.ek(j);
    ek -= this->projectionPt(ek);
    const auto n2(ek.dot(ek));
    if(n2 <= epsilon) continue;
    assert(0 <= idx[ii] && idx[ii] < this->rows());
    this->row(idx[ii ++]) = ek / sqrt(n2);
  }
  assert(idx.size() <= ii);
  return *this;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::QR() const {
  assert(0 < this->cols() && this->cols() <= this->rows());
  SimpleMatrix<T> Q(this->cols(), this->rows());
  Q.O();
  vector<int> residue;
  residue.reserve(Q.rows());
  for(int i = 0; i < Q.rows(); i ++) {
    const auto Atrowi(this->col(i));
    const auto work(Atrowi - Q.projectionPt(Atrowi));
    const auto n2(work.dot(work));
    if(n2 <= epsilon) {
      residue.emplace_back(i);
      continue;
    }
    Q.row(i) = work / sqrt(n2);
  }
  return Q.fillP(residue);
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::SVD() const {
  // N.B. S S^t == L Q Q^t L^t == L L^t.
  //      eigen value on L is same on S, S is symmetric.
  const auto S(*this * transpose());
  const auto R(S.QR() * S);
        auto Ut(SimpleMatrix<T>(S.rows(), S.cols()).O());
  vector<int> fill;
  fill.reserve(Ut.rows());
  // XXX: we can boost this function with
  //   S = counterI L counterI L compilicated form,
  // in such case, we should repalce rows by reverse order first,
  // then, some transpose is needed.
  //  after all we get (L - lambda I) x == 0 whole.
  for(int i = 0; i < Ut.rows(); i ++) {
    // singular index is same as vanish index.
    auto SS(S);
    SS = S - SS.I(R(i, i));
    const auto work(SimpleMatrix<T>(Ut.cols() - 1, Ut.cols() - 1).O().setMatrix(0, 0, SS.subMatrix(0, 0, i, i)).setMatrix(i, 0, SS.subMatrix(i + 1, 0, SS.rows() - i - 1, i)).setMatrix(0, i, SS.subMatrix(0, i + 1, i, SS.cols() - i - 1)).setMatrix(i, i, SS.subMatrix(i + 1, i + 1, SS.rows() - i - 1, SS.cols() - i - 1)).solve(SimpleVector<T>(Ut.cols() - 1).O().setVector(0, SS.col(i).subVector(0, i)).setVector(i, SS.col(i).subVector(i + 1, SS.rows() - i - 1))));
    Ut.row(i).setVector(0, work.subVector(0, i)).setVector(i + 1, work.subVector(i, work.size() - i));
    const auto n2(Ut.row(i).dot(Ut.row(i)));
    if(n2 <= epsilon)
      fill.emplace_back(i);
    else
      Ut.row(i) /= sqrt(n2);
  }
  return Ut.fillP(fill);
}

template <typename T> inline pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SimpleMatrix<T>::SVD(const SimpleMatrix<T>& src) const {
  // refered from : https://en.wikipedia.org/wiki/Generalized_singular_value_decomposition .
  assert(this->cols() == src.cols());
  SimpleMatrix<T> C(this->rows() + src.rows(), this->cols());
  C.setMatrix(0, 0, *this);
  C.setMatrix(this->rows(), 0, src);
  const auto P(C.SVD());
  SimpleVector<T> d(this->cols());
        auto Qt(P * C);
  vector<int> fill;
  fill.reserve(d.size());
  for(int i = 0; i < d.size(); i ++) {
    d[i] = Qt.row(i).dot(Qt.row(i));
    // XXX ratio:
    if(d[i] <= epsilon) {
      fill.emplace_back(i);
      d[i] = sqrt(d[i]);
    } else
      Qt.row(i) /= (d[i] = sqrt(d[i]));
  }
  Qt.fillP(fill);
  const auto D(P * C * Qt.transpose());
  SimpleMatrix<T> P1(this->rows(), d.size());
  SimpleMatrix<T> P2(src.rows(), d.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < P1.rows(); i ++)
    P1.row(i) = P.col(i);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < P2.rows(); i ++)
    P2.row(i) = P.col(i + P1.rows());
  auto U1(P1.SVD());
  auto Wt(U1 * P1);
  fill = vector<int>();
  fill.reserve(Wt.rows());
  for(int i = 0; i < Wt.rows(); i ++) {
    const auto n2(Wt.row(i).dot(Wt.row(i)));
    if(n2 <= epsilon)
      fill.emplace_back(i);
    else
      Wt.row(i) /= sqrt(n2);
  }
  Wt.fillP(fill);
  auto U2(Wt * P2.transpose());
  fill = vector<int>();
  fill.reserve(U2.rows());
  for(int i = 0; i < U2.rows(); i ++) {
    const auto n2(U2.row(i).dot(U2.row(i)));
    if(n2 <= epsilon)
      fill.emplace_back(i);
    else
      U2.row(i) /= sqrt(n2);
  }
  return make_pair(make_pair(move(U1), move(U2.fillP(fill))), (Wt * D).transpose().QR() * Qt);
}

template <typename T> template <typename U> inline SimpleMatrix<U> SimpleMatrix<T>::real() const {
  assert(0 < erows && 0 < ecols);
  SimpleMatrix<U> res(erows, ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    for(int j = 0; j < ecols; j ++)
      res(i, j) = U(entity[i][j].real());
  return res;
}

template <typename T> template <typename U> inline SimpleMatrix<U> SimpleMatrix<T>::imag() const {
  assert(0 < erows && 0 < ecols);
  SimpleMatrix<U> res(erows, ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    for(int j = 0; j < ecols; j ++)
      res(i, j) = U(entity[i][j].imag());
  return res;
}

template <typename T> template <typename U> inline SimpleMatrix<U> SimpleMatrix<T>::cast() const {
  assert(0 < erows && 0 < ecols);
  SimpleMatrix<U> res(erows, ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    for(int j = 0; j < ecols; j ++)
      res(i, j) = U(entity[i][j]);
  return res;
}

template <typename T> inline const int& SimpleMatrix<T>::rows() const {
  return erows;
}

template <typename T> inline const int& SimpleMatrix<T>::cols() const {
  return ecols;
}

template <typename T> inline void SimpleMatrix<T>::resize(const int& rows, const int& cols) {
  assert(0 <= rows && 0 <= cols);
  if(rows != erows) {
    erows = rows;
    if(entity)
      delete[] entity;
    if(0 < rows)
      entity = new SimpleVector<T>[erows];
    else
      entity = NULL;
    ecols = 0;
  }
  if(cols != ecols) {
    ecols = cols;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < erows; i ++)
      entity[i].resize(ecols);
  }
  return;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::innerFix(const SimpleMatrix<T>& A, vector<pair<T, int> >& fidx) {
  // N.B. we now have |[A -bb] [x t]| <= 1 condition.
  // N.B. there's no difference |[A - bb] [x t]|^2 <= 1 condition in this.
  //      but not with mixed condition.
  const auto R((*this) * A);
  SimpleVector<T> one(this->cols());
  one.I(T(1));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < fidx.size(); i ++)
    one[fidx[i].second] = T(0);
  // we now have: Q [R [x t] ] <= {0, 1}^m cond.
  const auto on(projectionPt(one));
  fidx.reserve(fidx.size() + this->cols());
  for(int i = 0; i < this->cols(); i ++)
    fidx.emplace_back(make_pair(abs(on[i]), i));
  sort(fidx.begin(), fidx.end());
  // sort by: |<Q^t(1), q_k>|, we subject to minimize each, to do this,
  //   maximize minimum q_k orthogonality.
  for(int i = 0, idx = 0; i < this->rows() - 1 && idx < fidx.size(); idx ++) {
    const auto& iidx(fidx[idx].second);
    const auto  orth(this->col(iidx));
    const auto  n2(orth.dot(orth));
    if(n2 <= epsilon)
      continue;
    // N.B. O(mn) can be writed into O(lg m + lg n) in many core cond.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < this->cols(); j ++)
      this->setCol(j, this->col(j) - orth * this->col(j).dot(orth) / n2);
    i ++;
  }
  // N.B. now we have fix indexes that to be P R [x 1] * t == 0.
  return R.solve((*this) * one);
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const {
  assert(this->rows() == bl.size() && this->rows() == bu.size() &&
         0 < this->cols() && 0 < this->rows() && this->cols() < this->rows());
  // bu - bb == A, bl - bb == - A <=> bu - bl == 2 A. 
  const auto bb((bu + bl) / T(2));
  const auto upper(bu - bb);
  SimpleMatrix<T> A(*this);
  vector<pair<T, int> > fidx;
  fidx.reserve(this->rows());
  for(int i = 0; i < this->rows(); i ++) {
    assert(T(0) <= upper[i]);
    if(upper[i] == T(0)) {
      const auto n2(A.row(i).dot(A.row(i)));
      if(n2 != T(0)) {
        fidx.emplace_back(make_pair(- T(1), i));
        A.row(i) /= sqrt(n2);
      } else
        A.row(i) *= n2;
    } else
      A.row(i) /= upper[i];
    assert(isfinite(A.row(i).dot(A.row(i))));
  }
  // N.B. we now have |[A -bb] [x t]| <= 1 condition.
  // N.B. there's no difference |[A -bb] [x t]|^2 <= 1 condition in this.
  //      but not with mixed condition.
  // N.B. we don't matter signs on [A -bb] [x t].
  return A.QR().innerFix(A, fidx);
}

template <typename T> std::ostream& operator << (std::ostream& os, const SimpleMatrix<T>& v) {
  SimpleMatrix<string> buf(v.rows(), v.cols());
  int M(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < v.rows(); i ++)
    for(int j = 0; j < v.cols(); j ++) {
      stringstream ss;
      ss << v(i, j);
      buf(i, j) = ss.str();
#if defined(_OPENMP)
#pragma omp critical
#endif
      {
        M = max(int(buf(i, j).size()), M);
      }
    }
  os << "(" << buf.rows() << ", " << buf.cols() << ")" << "[" << endl;
  for(int i = 0; i < buf.rows(); i ++) {
    os << "[";
    for(int j = 0; j < buf.cols(); j ++) {
      for(int k = buf(i, j).size(); k <= M; k ++)
        os << " ";
      os << buf(i, j);
      if(j < buf.cols() - 1) os << ", ";
    }
    os << "]";
    if(i < buf.rows() - 1) os << ", ";
    os << endl;
  }
  os << "]" << endl;
  return os;
}

template <typename T> std::istream& operator >> (std::istream& is, SimpleMatrix<T>& v) {
  while(! is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == '\n') continue;
    else if(c == '(') break;
    is.unget();
    break;
  }
  int r, c;
  is >> r;
  if(r <= 0) return is;
  while(! is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == '\n') continue;
    else if(c == ',') break;
    is.unget();
    break;
  }
  is >> c;
  while(! is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == '\n') continue;
    else if(c == ')') break;
    is.unget();
    break;
  }
  if(c <= 0) return is;
  v.resize(r, c);
  int i(0);
  int j(0);
  while(! is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == '\n') continue;
    else if(c == '[') break;
    is.unget();
    break;
  }
  for( ; i < v.rows() && ! is.eof() && ! is.bad(); i ++) {
    while(! is.eof() && ! is.bad()) {
      const auto c(is.get());
      if(c == ' ' || c == '\t' || c == '\n') continue;
      else if(c == '[') break;
      is.unget();
      break;
    }
    for(j = 0; j < v.cols() && ! is.eof() && ! is.bad(); j ++) {
      is >> v(i, j);
      if(v.cols() - 1 <= j) {
        j ++;
        break;
      }
      while(! is.eof() && ! is.bad()) {
        const auto c(is.get());
        if(c == ' ' || c == '\t' || c == '\n') continue;
        else if(c == ',') break;
        is.unget();
        break;
      }
    }
    while(! is.eof() && ! is.bad()) {
      const auto c(is.get());
      if(c == ' ' || c == '\t' || c == '\n') continue;
      else if(c == ']') break;
      is.unget();
      break;
    }
    if(v.rows() - 1 <= i) {
      i ++;
      break;
    }
    while(! is.eof() && ! is.bad()) {
      const auto c(is.get());
      if(c == ' ' || c == '\t' || c == '\n') continue;
      else if(c == ',') break;
      is.unget();
      break;
    }
  }
  while(! is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c =='\n') continue;
    else if(c == ']') break;
    cerr << "XXX SimpleMatrix<T>::operator >> (\']\')" << flush;
    is.unget();
    break;
  }
  while(! is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t') continue;
    else if(c == '\n') break;
    cerr << "XXX SimpleMatrix<T>::operator >> (\'\\n\')" << flush;
    is.unget();
    break;
  }
  if(i < v.rows() || j < v.cols()) {
    cerr << "XXX SimpleMatrix<T>::operator >> (index)" << flush;
    for( ; i < v.rows(); i ++)
      for( ; j < v.cols(); j ++)
        v(i, j) = T(0);
  }
  return is;
}


template <typename T> SimpleMatrix<complex<T> > dft(const int& size0) {
  const auto size(abs(size0));
  if(! size) {
    const static SimpleMatrix<complex<T> > m0;
    return m0;
  }
  SimpleMatrix<complex<T> > edft( size, size);
  SimpleMatrix<complex<T> > eidft(size, size);
  const auto file(std::string("./.cache/lieonn/dft-") + std::to_string(size) +
#if defined(_FLOAT_BITS_)
    std::string("-") + std::to_string(_FLOAT_BITS_)
#else
    std::string("-ld")
#endif
  );
  ifstream cache(file.c_str());
  if(cache.is_open()) {
    cache >> edft;
    cache >> eidft;
    cache.close();
  } else {
    static const auto Pi(T(4) * atan2(T(1), T(1)));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < edft.rows(); i ++) {
      for(int j = 0; j < edft.cols(); j ++) {
        const auto theta(- T(2) * Pi * T(i) * T(j) / T(edft.rows()));
        const auto c(cos(theta));
        const auto s(sin(theta));
        edft( i, j) = complex<T>(c,   s);
        eidft(i, j) = complex<T>(c, - s) / complex<T>(T(size));
      }
    }
    ofstream ocache(file.c_str());
    ocache << edft;
    ocache << eidft;
    ocache.close();
  }
  return size0 < 0 ? eidft : edft;
}

template <typename T> SimpleMatrix<T> diff(const int& size0) {
  const auto size(abs(size0));
  if(! size) {
    static const SimpleMatrix<T> m0;
    return m0;
  }
  SimpleMatrix<T> dd;
  SimpleMatrix<T> ii;
  const auto file(std::string("./.cache/lieonn/diff-") + std::to_string(size) +
#if defined(_FLOAT_BITS_)
    std::string("-") + std::to_string(_FLOAT_BITS_)
#else
    std::string("-ld")
#endif
  );
  ifstream cache(file.c_str());
  if(cache.is_open()) {
    cache >> dd;
    cache >> ii;
    cache.close();
  } else {
    if(2 < size) {
      const auto dd0(diff<T>(  (size - 1)) * T(size - 1));
      const auto ii0(diff<T>(- (size - 1)) * T(size - 1));
      dd = SimpleMatrix<T>(size, size).O();
      ii = SimpleMatrix<T>(size, size).O();
            auto ddo(SimpleMatrix<T>(size, size).O().setMatrix(0, 0, dd0));
            auto iio(SimpleMatrix<T>(size, size).O().setMatrix(0, 0, ii0));
      ddo.setMatrix(1, 1, dd.subMatrix(1, 1, size - 1, size - 1) + dd0);
      iio.setMatrix(1, 1, ii.subMatrix(1, 1, size - 1, size - 1) + ii0);
      ddo.row(0) *= T(2);
      ddo.row(ddo.rows() - 1) *= T(2);
      iio.row(0) *= T(2);
      iio.row(iio.rows() - 1) *= T(2);
      dd = ddo / T(2);
      ii = iio / T(2);
    } else {
      dd = SimpleMatrix<T>(size, size).O();
      ii = SimpleMatrix<T>(size, size).O();
    }
    auto DD(dft<T>(size));
    auto II(dft<T>(size));
    static const auto Pi(T(4) * atan2(T(1), T(1)));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < DD.rows(); i ++)
      DD.row(i) *=   complex<T>(T(0), - T(2) * Pi * T(i) / T(DD.rows()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 1; i < II.rows(); i ++)
      II.row(i) /= - complex<T>(T(0), - T(2) * Pi * T(i) / T(DD.rows()));
    // N.B. in continuous function, we don't divide dd by &pi;.
    //      (d/dx sum exp(2 Pi i x theta / N) exp(- 2 Pi i y theta / N) f(y)).
    // XXX: from some numerical test, sign on DD, II is reverse side.
    // XXX: In discrete function, we had be chosen ||dd|| := 1, and not now.
    //      (because the (left or right) differential itself cannot be
    //       larger than |x_{k+1}-x_k| < ||x|| in discrete,
    //       sum_0^1 - 2 Pi i (theta/n)^2/2 -&gt; Pi)
    // XXX: if we make plain differential with no error on cosine curve,
    //      it causes constant 0 vector.
    dd += (dft<T>(- size) * DD).template real<T>();
    ii += (dft<T>(- size) * II).template real<T>();
    if(3 < size) {
      dd /= T(size);
      ii /= T(size);
    }
    ofstream ocache(file.c_str());
    ocache << dd;
    ocache << ii;
    ocache.close();
    cerr << "." << flush;
  }
  return size0 < 0 ? ii : dd;
}

template <typename T> inline SimpleVector<T> taylor(const int& size, const T& step) {
  const int  step00(max(int(0), min(size - 1, int(floor(step)))));
  const auto residue0(step - T(step00));
  const auto step0(step00 == size - 1 || abs(residue0) <= T(1) / T(2) ? step00 : step00 + 1);
  const auto residue(step - T(step0));
  SimpleVector<T> res(size);
  res.ek(step0);
  if(residue == T(0)) return res;
  const auto Dt(diff<T>(size).transpose());
        auto dt(Dt.col(step0) * residue);
  // N.B.
  // if we deal with (D *= r, residue /= r), it is identical with (D, residue)
  // So ||D^n * residue^n|| / T(n!) < 1 case, this loop converges.
  // but with n^n v.s. n!, differential of n! is faster than n^n.
  // (n! < n^n but a^n < n! somewhere).
  // And, we treat D * residue as a block, so Readme.md's condition 1/x^k needs
  // to be in the series in this.
  for(int i = 2; ; i ++) {
    const auto last(res);
    res += dt;
    if(last == res) break;
    dt   = Dt * dt * residue / T(i);
  }
  return res;
}

template <typename T> SimpleVector<T> linearInvariant(const vector<SimpleVector<T> >& in, const bool& nonzero = false) {
  SimpleMatrix<T> A(in.size(), in[0].size());
  assert(in[0].size() <= in.size());
  SimpleVector<T> fvec(A.cols());
  SimpleVector<T> one(A.rows());
  fvec.O();
  one.I();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < in.size(); i ++) {
    assert(in[i].size() == in[0].size());
    for(int j = 0; j < in[0].size(); j ++)
      A(i, j) = in[i][j];
    assert(isfinite(A.row(i).dot(A.row(i))));
  }
  vector<pair<T, int> > sute;
  return A.QR().innerFix(A, sute);
}

// N.B. please refer bitsofcotton/randtools.
template <typename T> inline pair<SimpleVector<T>, T> makeProgramInvariant(const SimpleVector<T>& in, const T& index = - T(1)) {
  SimpleVector<T> res(in.size() + (T(0) <= index ? 3 : 2));
  res.setVector(0, in);
  res[in.size()] = T(1);
  if(T(0) <= index)
    res[in.size() + 1] = T(index);
  res[res.size() - 1] = T(0);
  T   lsum(0);
  for(int i = 0; i < res.size() - 1; i ++) {
    assert(- T(1) <= res[i] && res[i] <= T(1));
    res[i] += T(1);
    if(res[i] != T(0)) lsum += log(res[i]);
  }
  T ratio(1);
  // N.B. x_1 ... x_n == 1.
  // <=> x_1 / (x_1 ... x_n)^(1/n) ... == 1.
  if(lsum != T(0)) res /= ratio = exp(lsum / T(res.size() - 1));
  res[res.size() - 1] = T(1);
  return make_pair(res, ratio);
}

template <typename T> inline T revertProgramInvariant(const pair<T, T>& in) {
  return in.first * in.second - T(1);
}

template <typename T> class linearFeeder {
public:
  inline linearFeeder() { t = 0; full = false; }
  inline linearFeeder(const int& size) {
    res.resize(size);
    for(int i = 0; i < res.size(); i ++)
      res[i] = T(0);
    t = 0;
    full = false;
  }
  inline ~linearFeeder() { ; }
  inline const SimpleVector<T>& next(const T& in) {
    for(int i = 0; i < res.size() - 1; i ++)
      res[i] = move(res[i + 1]);
    res[res.size() - 1] = in;
    if(res.size() < t ++) full = true;
    return res;
  }
  SimpleVector<T> res;
  bool full;
private:
  int t;
};

template <typename T> class arctanFeeder {
public:
  inline arctanFeeder() { t = 0; full = false; }
  inline arctanFeeder(const int& size) {
    res.resize(size);
    buf.resize(1 + int(ceil(T(1) / tan(T(1) * atan(T(1)) / T(res.size() - 1)))));
    for(int i = 0; i < buf.size(); i ++)
      buf[i] = T(0);
    for(int i = 0; i < res.size(); i ++)
      res[i] = T(0);
    t = 0;
    full = false;
  }
  inline ~arctanFeeder() { ; }
  inline const SimpleVector<T>& next(const T& in) {
    for(int i = 0; i < buf.size() - 1; i ++)
      buf[i] = move(buf[i + 1]);
    buf[buf.size() - 1] = in;
    if(buf.size() < t ++) full = true;
    for(int i = 0; i < res.size(); i ++)
      res[res.size() - i - 1] = buf[buf.size() - 1 - int(tan(T(i) * atan(T(1)) / T(res.size() - 1)) / tan(T(1) * atan(T(1)) / T(res.size() - 1)))];
    return res;
  }
  SimpleVector<T> res;
  bool full;
private:
  SimpleVector<T> buf;
  int t;
};

template <typename T> class SimpleSparseVector {
public:
  inline SimpleSparseVector();
  inline SimpleSparseVector(const int& sute);
  inline SimpleSparseVector(const SimpleSparseVector<T>& other);
  inline SimpleSparseVector(SimpleSparseVector<T>&& other);
  inline ~SimpleSparseVector();
  
  inline SimpleSparseVector<T>  operator -  () const;
  inline SimpleSparseVector<T>  operator +  (const SimpleSparseVector<T>& other) const;
  inline SimpleSparseVector<T>& operator += (const SimpleSparseVector<T>& other);
  inline SimpleSparseVector<T>  operator -  (const SimpleSparseVector<T>& other) const;
  inline SimpleSparseVector<T>& operator -= (const SimpleSparseVector<T>& other);
  template <typename U> inline SimpleSparseVector<T>  operator *  (const U& other) const;
  template <typename U> inline SimpleSparseVector<T>& operator *= (const U& other);
  template <typename U> inline SimpleSparseVector<T>  operator /  (const U& other) const;
  template <typename U> inline SimpleSparseVector<T>& operator /= (const U& other);
  inline       SimpleSparseVector<T>& operator =  (const SimpleSparseVector<T>& other);
  inline       SimpleSparseVector<T>& operator =  (SimpleSparseVector<T>&& other);
  inline       bool                   operator != (const SimpleSparseVector<T>& other) const;
  inline       bool                   operator == (const SimpleSparseVector<T>& other) const;
  inline       T  dot         (const SimpleSparseVector<T>& other) const;
  inline       T& operator [] (const int& idx);
  inline const T& operator [] (const int& idx) const;
  inline void     clear();
  inline       map<int, T>& iter();
  inline const map<int, T>& iter() const;
private:
  map<int, T>  entity;
};

template <typename T> inline SimpleSparseVector<T>::SimpleSparseVector() {
  return;
}

template <typename T> inline SimpleSparseVector<T>::SimpleSparseVector(const int& sute) {
  assert(sute == 0);
  return;
}

template <typename T> inline SimpleSparseVector<T>::SimpleSparseVector(const SimpleSparseVector<T>& other) {
  *this = other;
}

template <typename T> inline SimpleSparseVector<T>::SimpleSparseVector(SimpleSparseVector<T>&& other) {
  *this = other;
}

template <typename T> inline SimpleSparseVector<T>::~SimpleSparseVector() {
  return;
}

template <typename T> inline SimpleSparseVector<T> SimpleSparseVector<T>::operator - () const {
  auto res(*this);
  for(auto itr(res.entity.begin()); itr != res.entity.end(); ++ itr)
    itr->second = - itr->second;
  return res;
}

template <typename T> inline SimpleSparseVector<T> SimpleSparseVector<T>::operator + (const SimpleSparseVector<T>& other) const {
  auto res(*this);
  return res += other;
}

template <typename T> inline SimpleSparseVector<T>& SimpleSparseVector<T>::operator += (const SimpleSparseVector<T>& other) {
  for(auto itr(other.entity.begin()); itr != other.entity.end(); ++ itr) {
    if(itr->second == T(0)) continue;
    auto search(entity.lower_bound(itr->first));
    if(search == entity.end() || search->first != itr->first)
      (*this)[itr->first] = itr->second;
    else
      search->second += itr->second;
  }
  return *this;
}

template <typename T> inline SimpleSparseVector<T> SimpleSparseVector<T>::operator - (const SimpleSparseVector<T>& other) const {
  auto res(*this);
  return res -= other;
}

template <typename T> inline SimpleSparseVector<T>& SimpleSparseVector<T>::operator -= (const SimpleSparseVector<T>& other) {
  return *this += - other;
}

template <typename T> template <typename U> inline SimpleSparseVector<T> SimpleSparseVector<T>::operator * (const U& other) const {
  auto res(*this);
  return res *= other;
}

template <typename T> template <typename U> inline SimpleSparseVector<T>& SimpleSparseVector<T>::operator *= (const U& other) {
  for(auto itr(entity.begin()); itr != entity.end(); ++ itr)
    itr->second *= other;
  return *this;
}

template <typename T> template <typename U> inline SimpleSparseVector<T> SimpleSparseVector<T>::operator / (const U& other) const {
  auto res(*this);
  return res /= other;
}

template <typename T> inline SimpleSparseVector<T>& SimpleSparseVector<T>::operator = (const SimpleSparseVector<T>& other) {
  entity = other.entity;
  return *this;
}

template <typename T> inline SimpleSparseVector<T>& SimpleSparseVector<T>::operator = (SimpleSparseVector<T>&& other) {
  entity = move(other.entity);
  return *this;
}

template <typename T> inline bool SimpleSparseVector<T>::operator != (const SimpleSparseVector<T>& other) const {
  for(auto itr(entity.begin()); itr != entity.end(); ++ itr)
    if(itr->second != other[itr->first])
      return true;
  for(auto itr(other.entity.begin()); itr != other.entity.end(); ++ itr)
    if(itr->second != const_cast<const SimpleSparseVector<T>&>(*this)[itr->first])
      return true;
  return false;
}

template <typename T> inline bool SimpleSparseVector<T>::operator == (const SimpleSparseVector<T>& other) const {
  return ! (*this != other);
}

template <typename T> template <typename U> inline SimpleSparseVector<T>& SimpleSparseVector<T>::operator /= (const U& other) {
  for(auto itr(entity.begin()); itr != entity.end(); ++ itr)
    itr->second /= other;
  return *this;
}

template <typename T> inline T SimpleSparseVector<T>::dot(const SimpleSparseVector<T>& other) const {
  T res(0);
  for(auto itr(other.entity.begin()); itr < other.entity.end(); ++ itr) {
    auto search(entity.lower_bound(itr->first));
    if(search != entity.end() && search->first == itr->first)
      res += search->second * itr->second;
  }
  return res;
}

template <typename T> inline T& SimpleSparseVector<T>::operator [] (const int& idx) {
  assert(0 <= idx);
  const auto search(entity.lower_bound(idx));
  if(search != entity.end() && search->first == idx)
    return search->second;
  else
    entity[idx] = T(0);
  const auto search2(entity.lower_bound(idx));
  assert(search2 != entity.end() && search2->first == idx);
  return search2->second;
}

template <typename T> inline const T& SimpleSparseVector<T>::operator [] (const int& idx) const {
  assert(0 <= idx);
  if(entity.size()) {
    const auto search(entity.lower_bound(idx));
    if(search != entity.end() && search->first == idx)
      return search->second;
  }
  const static T zero(0);
  return zero;
}

template <typename T> inline void SimpleSparseVector<T>::clear() {
  entity = map<int, T>();
  return;
}

template <typename T> inline map<int, T>& SimpleSparseVector<T>::iter() {
  return entity;
}

template <typename T> inline const map<int, T>& SimpleSparseVector<T>::iter() const {
  return entity;
}

template <typename T> using SimpleSparseMatrix = SimpleSparseVector<SimpleSparseVector<T> >;
template <typename T> using SimpleSparseTensor = SimpleSparseVector<SimpleSparseVector<SimpleSparseVector<T> > >;

#define _SIMPLELIN_
#endif


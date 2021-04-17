/*
BSD 3-Clause License

Copyright (c) 2019-2021, bitsofcotton
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

#if !defined(_INTEGER_FLOAT_)

#if !defined(_FLOAT_BITS_)
  #include <complex>
  #include <cmath>
  using namespace std;
  typedef uint64_t myuint;
  typedef int64_t  myint;
  typedef long double myfloat;
  #define mybits 64
#else

using std::move;
using std::max;
using std::min;
using std::vector;

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
  inline DUInt<T,bits>  operator ++ (int);
  inline DUInt<T,bits>& operator -- ();
  inline DUInt<T,bits>  operator -- (int);
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
  inline DUInt<T,bits>& operator =  (const int& src);
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
  e[0]   = abssrc;
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

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator ++ (int) {
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

template <typename T, int bits> inline DUInt<T,bits>  DUInt<T,bits>::operator -- (int) {
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
  else if(b > bits * 2)
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
  else if(b > bits * 2)
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

template <typename T, int bits> inline DUInt<T,bits>& DUInt<T,bits>::operator =  (const int& src) {
  e[0]  = src;
  e[1] ^= e[1];
  return *this;
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
    } else
      goto ensure;
  }
 ensure:
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

  const vector<SimpleFloat<T,W,bits,U> >& exparray()    const;
  const vector<SimpleFloat<T,W,bits,U> >& invexparray() const;
};

template <typename T, typename W, int bits, typename U> inline SimpleFloat<T,W,bits,U>::SimpleFloat() {
  assert(0 < bits && ! (bits & 1));
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

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>
& SimpleFloat<T,W,bits,U>::halfpi() const {
  const static auto vhalfpi(quatpi() << U(1));
  return vhalfpi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>
& SimpleFloat<T,W,bits,U>::quatpi() const {
  const static auto vquatpi(one().atan());
  return vquatpi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>
& SimpleFloat<T,W,bits,U>::twopi() const {
  const static auto vtwopi(quatpi() << U(3));
  return vtwopi;
}

template <typename T, typename W, int bits, typename U> const SimpleFloat<T,W,bits,U>
& SimpleFloat<T,W,bits,U>::sqrt2() const {
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
    return os << "NaN";
  if(isinf(v))
    return os << (const char*)(v.s & (1 << v.SIGN) ? "-" : "") << "Inf";
  return os << (const char*)(v.s & (1 << v.SIGN) ? "-" : "") << T(v.m) << "*2^" << v.e;
}

template <typename T, typename W, int bits, typename U> std::istream& operator >> (std::istream& is, SimpleFloat<T,W,bits,U>& v) {
  const static SimpleFloat<T,W,bits,U> ten(10);
               SimpleFloat<T,W,bits,U> e(0);
  bool mode(false);
  bool sign(false);
  bool fsign(false);
  v = SimpleFloat<T,W,bits,U>(0);
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
    switch(buf) {
    case '-':
      sign  = true;
    case '+':
      if(fsign)
        throw "Wrong input";
      fsign = true;
      break;
    case 'e':
      if(mode)
        goto ensure;
      if(sign)
        v   = - v;
      mode  = true;
      sign  = false;
      fsign = false;
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
  v *= pow(ten, e);
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
  #define mybits 8
# elif _FLOAT_BITS_ == 16
  typedef uint16_t myuint;
  typedef int16_t  myint;
  typedef SimpleFloat<myuint, uint32_t, 16, int64_t> myfloat;
  #define mybits 16
# elif _FLOAT_BITS_ == 32
  typedef uint32_t myuint;
  typedef int32_t  myint;
  typedef SimpleFloat<myuint, uint64_t, 32, int64_t> myfloat;
  #define mybits 32
# elif _FLOAT_BITS_ == 64
  typedef uint64_t myuint;
  typedef int64_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 64>, 64, int64_t> myfloat;
  #define mybits 64
# elif _FLOAT_BITS_ == 128
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef Signed<uint128_t, 128> int128_t;
  typedef uint128_t myuint;
  typedef int128_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 128>, 128, int64_t> myfloat;
  #define mybits 128
# elif _FLOAT_BITS_ == 256
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef DUInt<uint128_t, 128> uint256_t;
  typedef Signed<uint256_t, 128> int256_t;
  typedef uint256_t myuint;
  typedef int256_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 256>, 256, int64_t> myfloat;
  #define mybits 256
# elif _FLOAT_BITS_ == 512
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef DUInt<uint128_t, 128> uint256_t;
  typedef DUInt<uint256_t, 128> int256_t;
  typedef Signed<uint512_t, 128> int1024_t;
  typedef uint512_t myuint;
  typedef int512_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 512>, 512, int64_t> myfloat;
  #define mybits 256
# else
#   error cannot handle float
# endif
#endif

#define _INTEGER_FLOAT_
#endif


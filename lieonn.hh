/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/*
BSD 3-Clause License

Copyright (c) 2013-2023, bitsofcotton (kazunobu watatsu)
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
using std::vector;
using std::max;
using std::min;
using std::swap;
using std::sort;
using std::binary_search;
using std::pair;
using std::make_pair;
using std::map;

using std::stringstream;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;

using std::string;
using std::to_string;
using std::cerr;
using std::endl;
using std::flush;
using std::getline;


// Double int to new int class.
template <typename T, int bits> class DUInt {
public:
  inline DUInt(const int& src = 0) {
    assert(0 < bits && ! (bits & 3));
    const auto abssrc(src < 0 ? - src : src);
    e[0]   = T(abssrc);
    e[1]  ^= e[1];
    if(abssrc != src)
      *this = - *this;
  }
  inline DUInt(const T& src) {
    const auto abssrc(src < T(int(0)) ? - src : src);
    e[0]   = abssrc;
    e[1]  ^= e[1];
    if(abssrc != src)
      *this = - *this;
  }
  inline DUInt(const DUInt<T,bits>& src) { *this = src; }
  inline DUInt(const DUInt<DUInt<T,bits>,bits*2>& src) { *this = src; }
  inline DUInt(DUInt<T,bits>&& src) { *this = src; }
  inline ~DUInt() { ; }
  
  inline DUInt<T,bits>& operator ++ () {
    ++ e[0];
    if(!e[0]) ++ e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator ++ (int32_t) {
    const auto work(*this);
    ++ *this;
    return work;
  }
  inline DUInt<T,bits>& operator -- () {
    if(!e[0]) -- e[1];
    -- e[0];
    return *this;
  }
  inline DUInt<T,bits>  operator -- (int32_t) {
    const auto work(*this);
    -- *this;
    return work;
  }
  inline DUInt<T,bits>  operator -  () const {
    auto work(~ *this);
    return ++ work;
  }
  inline DUInt<T,bits>  operator +  (const DUInt<T,bits>& src) const {
    auto work(*this);
    return work += src;
  }
  inline DUInt<T,bits>& operator += (const DUInt<T,bits>& src) {
    // N.B. assembler can boost dramatically this code. but not here.
    const auto e0(max(e[0], src.e[0]));
    e[0] += src.e[0];
    if(e[0] < e0)
      e[1] ++;
    e[1] += src.e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator -  (const DUInt<T,bits>& src) const {
    auto work(*this);
    return work -= src;
  }
  inline DUInt<T,bits>& operator -= (const DUInt<T,bits>& src) {
    return *this += - src;
  }
  inline DUInt<T,bits>  operator *  (const DUInt<T,bits>& src) const {
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
  inline DUInt<T,bits>& operator *= (const DUInt<T,bits>& src) {
    return *this = *this * src;
  }
  inline DUInt<T,bits>  operator /  (const DUInt<T,bits>& src) const {
    auto work(*this);
    return work /= src;
  }
  inline DUInt<T,bits>& operator /= (const DUInt<T,bits>& src) {
    const static DUInt<T,bits> one(int(1));
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
  inline DUInt<T,bits>  operator %  (const DUInt<T,bits>& src) const {
    return *this - ((*this / src) * src);
  }
  inline DUInt<T,bits>& operator %= (const DUInt<T,bits>& src) {
    return *this = *this % src;
  }
  inline DUInt<T,bits>  operator << ( const int& b)            const {
    auto work(*this);
    return work <<= b;
  }
  inline DUInt<T,bits>& operator <<= (const int& b) {
    if(! b)                return *this;
    else if(b < 0)         return *this >>= (- b);
    else if(b >= bits * 2) return *this ^= *this;
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
  inline DUInt<T,bits>  operator >> ( const int& b)            const {
    auto work(*this);
    return work >>= b;
  }
  inline DUInt<T,bits>& operator >>= (const int& b) {
    if(! b)                return *this;
    else if(b < 0)         return *this <<= (- b);
    else if(b >= bits * 2) return *this ^= *this;
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
  inline DUInt<T,bits>  operator &  (const DUInt<T,bits>& src) const {
    auto work(*this);
    return work &= src;
  }
  inline DUInt<T,bits>& operator &= (const DUInt<T,bits>& src) {
    e[0] &= src.e[0]; e[1] &= src.e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator |  (const DUInt<T,bits>& src) const {
    auto work(*this);
    return work |= src;
  }
  inline DUInt<T,bits>& operator |= (const DUInt<T,bits>& src) {
    e[0] |= src.e[0]; e[1] |= src.e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator ^  (const DUInt<T,bits>& src) const {
    auto work(*this);
    return work ^= src;
  }
  inline DUInt<T,bits>& operator ^= (const DUInt<T,bits>& src) {
    e[0] ^= src.e[0]; e[1] ^= src.e[1];
    return *this;
  }
  inline DUInt<T,bits>  operator ~  ()                         const {
    DUInt<T,bits> work;
    work.e[0] = ~ e[0]; work.e[1] = ~ e[1];
    return work;
  }
  inline DUInt<T,bits>& operator =  (const DUInt<T,bits>& src) {
    e[0] = src.e[0]; e[1] = src.e[1];
    return *this;
  }
  inline DUInt<T,bits>& operator =  (const DUInt<DUInt<T,bits>,bits*2>& src) {
    return *this = src.e[0];
  }
  inline DUInt<T,bits>& operator =  (DUInt<T,bits>&& src) {
    e[0] = move(src.e[0]); e[1] = move(src.e[1]);
    return *this;
  }
  inline bool           operator <  (const DUInt<T,bits>& src) const {
    if(e[1]) return e[1] != src.e[1] ? e[1] < src.e[1] : e[0] < src.e[0];
    return bool(src.e[1]) || e[0] < src.e[0];
  }
  inline bool           operator <= (const DUInt<T,bits>& src) const {
    return *this < src || *this == src;
  }
  inline bool           operator >  (const DUInt<T,bits>& src) const {
    return ! (*this <= src);
  }
  inline bool           operator >= (const DUInt<T,bits>& src) const {
    return ! (*this < src);
  }
  inline bool           operator == (const DUInt<T,bits>& src) const {
    return ! (*this != src);
  }
  inline bool           operator != (const DUInt<T,bits>& src) const {
    return (*this ^ src).operator bool();
  }
  inline bool           operator && (const DUInt<T,bits>& src) const {
    return this->operator bool() && src.operator bool();
  }
  inline bool           operator || (const DUInt<T,bits>& src) const {
    return this->operator bool() || src.operator bool();
  }
  inline bool           operator !    () const {
    return ! this->operator bool();
  }
  inline                operator bool () const {
    return e[0] || e[1];
  }
  inline                operator int  () const {
    return int(e[0]);
  }
  inline                operator T    () const {
    return e[0];
  }
  inline                operator DUInt<T,bits> () const {
    return *this;
  }

  T e[2];
};

template <typename T, int bits> ostream&  operator << (ostream& os, DUInt<T,bits> v) {
  static const char* table = "0123456789abcdef";
  vector<char> buf;
  while(v) {
    buf.emplace_back(table[int(v) & 0x0f]);
    v >>= 4;
  }
  if(buf.size()) {
    for(int i = 0; 0 <= i && i < buf.size(); i ++)
      os << char(buf[buf.size() - 1 - i]);
    return os;
  }
  return os << '0';
}

template <typename T, int bits> istream&  operator >> (istream& is, DUInt<T,bits>& v) {
  v ^= v;
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
      v <<= 4;
      v |= DUInt<T,bits>(int(buf - '0'));
    } else if('a' <= buf && buf <= 'f') {
      v <<= 4;
      v |= DUInt<T,bits>(int(buf - 'a' + 10));
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
  inline Signed() { ; }
  inline Signed(const int& src) {
    T tsrc(src);
    *this = reinterpret_cast<const Signed<T,bits>&>(tsrc);
  }
  inline Signed(const T& src) {
    *this = reinterpret_cast<const Signed<T,bits>&>(src);
  }
  inline Signed(const Signed<T,bits>& src) {
    *this = src;
  }
  inline bool operator <  (const Signed<T,bits>& src) const {
    const auto mthis(int(*this >> (bits - 1)));
    const auto msrc( int(src   >> (bits - 1)));
    if(mthis ^ msrc) return mthis;
    if(mthis)
      return - dynamic_cast<const T&>(src) < - dynamic_cast<const T&>(*this);
    return dynamic_cast<const T&>(*this) < dynamic_cast<const T&>(src);
  }
  inline bool operator <= (const Signed<T,bits>& src) const {
    return ! (*this > src);
  }
  inline bool operator >  (const Signed<T,bits>& src) const {
    return ! (*this < src) && *this != src;
  }
  inline bool operator >= (const Signed<T,bits>& src) const {
    return ! (*this < src);
  }
};

template <typename T, int bits> ostream& operator << (ostream& os, Signed<T,bits> v) {
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
  inline SimpleFloat() {
    assert(0 < bits && ! (bits & 1));
    s |= (1 << NaN) | (1 << INF);
  }
  template <typename V> inline SimpleFloat(const V& src) {
    const static V vzero(0);
    s ^= s;
    m  = T(int(src < vzero ? - src : src));
    e ^= e;
    s |= safeAdd(e, normalize(m));
    if(src < vzero)
      s |= 1 << SIGN;
    ensureFlag();
  }
  inline SimpleFloat(const SimpleFloat<T,W,bits,U>& src) { *this = src; }
  inline SimpleFloat(SimpleFloat<T,W,bits,U>&& src) { *this = src; }
  inline ~SimpleFloat() { ; }
  
  inline SimpleFloat<T,W,bits,U>  operator -  () const {
    auto work(*this);
    work.s ^= 1 << SIGN;
    return work;
  }
  inline SimpleFloat<T,W,bits,U>  operator +  (const SimpleFloat<T,W,bits,U>& src) const {
    auto work(*this);
    return work += src;
  }
         SimpleFloat<T,W,bits,U>& operator += (const SimpleFloat<T,W,bits,U>& src) {
    if((s |= src.s & (1 << NaN)) & (1 << NaN)) return *this;
    if(s & (1 << INF)) {
      if((src.s & (1 << INF)) && (s ^ src.s) & (1 << SIGN)) s |= 1 << NaN;
      return *this;
    }
    if(src.s & (1 << INF)) return *this = src;
    if(! m) return *this = src;
    if(! src.m) return *this;
    if(! ((s ^ src.s) & (1 << SIGN))) {
      if(e >= src.e) {
        m >>= 1;
        s |= safeAdd(e, 1);
        U se(e);
        if(! safeAdd(se, - src.e) && se < U(bits))
          m += se ? src.m >> int(se) : src.m;
      } else
        return *this = src + *this;
    } else {
      if(e > src.e) {
        U se(e);
        if(! safeAdd(se, - src.e) && se < U(bits))
          m -= se ? src.m >> int(se) : src.m;
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
  inline SimpleFloat<T,W,bits,U>  operator -  (const SimpleFloat<T,W,bits,U>& src) const {
    auto work(*this);
    return work -= src;
  }
  inline SimpleFloat<T,W,bits,U>& operator -= (const SimpleFloat<T,W,bits,U>& src) {
    s ^= 1 << SIGN;
    *this += src;
    s ^= 1 << SIGN;
    return *this;
  }
  inline SimpleFloat<T,W,bits,U>  operator *  (const SimpleFloat<T,W,bits,U>& src) const {
    auto work(*this);
    return work *= src;
  }
         SimpleFloat<T,W,bits,U>& operator *= (const SimpleFloat<T,W,bits,U>& src) {
    s ^= src.s & (1 << SIGN);
    if((s |= src.s & (1 << NaN)) & (1 << NaN)) return *this;
    if((! m) || (! src.m)) {
      s |= 1 << DWRK;
      return ensureFlag();
    }
    if((s |= src.s & (1 << INF)) & (1 << INF)) return *this;
    auto mm(W(m) * W(src.m));
    s |= safeAdd(e, src.e);
    s |= safeAdd(e, normalize(mm));
    s |= safeAdd(e, U(bits));
    m  = T(mm >> bits);
    return ensureFlag();
  }
  inline SimpleFloat<T,W,bits,U>  operator /  (const SimpleFloat<T,W,bits,U>& src) const {
    auto work(*this);
    return work /= src;
  }
         SimpleFloat<T,W,bits,U>& operator /= (const SimpleFloat<T,W,bits,U>& src) {
    s ^= src.s & (1 << SIGN);
    if((s |= src.s & (1 << NaN)) & (1 << NaN)) return *this;
    if(! (s & (1 << INF)) && (src.s & (1 << INF))) {
      s |= 1 << DWRK;
      return ensureFlag();
    }
    if(s & (1 << INF)) {
      if(src.s & (1 << INF)) s |= 1 << NaN;
      return *this;
    }
    if(! src.m) {
      throw "Zero division";
      s |= 1 << NaN;
      return *this;
    }
    if(! m) return *this;
    auto mm((W(m) << bits) / W(src.m));
    s |= safeAdd(e, - src.e);
    s |= safeAdd(e, normalize(mm));
    m  = T(mm >> bits);
    return ensureFlag();
  }

  inline SimpleFloat<T,W,bits,U>  operator %  (const SimpleFloat<T,W,bits,U>& src) const {
    return *this - (*this / src).floor() * src;
  }
  inline SimpleFloat<T,W,bits,U>& operator %= (const SimpleFloat<T,W,bits,U>& src) {
    return *this = *this % src;
  }
  inline SimpleFloat<T,W,bits,U>  operator <<  (const U& b) const {
    auto work(*this);
    return work <<= b;
  }
  inline SimpleFloat<T,W,bits,U>& operator <<= (const U& b) {
    if(s & ((1 << INF) | (1 << NaN))) return *this;
    s |= safeAdd(e, b);
    return ensureFlag();
  }
  inline SimpleFloat<T,W,bits,U>  operator >>  (const U& b) const {
    auto work(*this);
    return work >>= b;
  }
  inline SimpleFloat<T,W,bits,U>& operator >>= (const U& b) {
    if(s & ((1 << INF) | (1 << NaN))) return *this;
    s |= safeAdd(e, - b);
    return ensureFlag();
  }
  inline SimpleFloat<T,W,bits,U>& operator =  (const SimpleFloat<T,W,bits,U>& src) {
    s = src.s;
    e = src.e;
    m = src.m;
    return *this;
  }
  inline SimpleFloat<T,W,bits,U>& operator =  (SimpleFloat<T,W,bits,U>&& src) {
    s = move(src.s);
    e = move(src.e);
    m = move(src.m);
    return *this;
  }
  inline bool             operator == (const SimpleFloat<T,W,bits,U>& src) const {
    return ! (*this != src);
  }
  inline bool             operator != (const SimpleFloat<T,W,bits,U>& src) const {
    return (((s | src.s) & ((1 << INF) | (1 << NaN))) ||
             (s != src.s || e != src.e || m != src.m)) &&
           ! (! m && ! src.m);
  }
  inline bool             operator <  (const SimpleFloat<T,W,bits,U>& src) const {
    if((s | src.s) & (1 << NaN)) throw "compair NaN";
    const auto s_is_minus(s & (1 << SIGN));
    if(s_is_minus ^ (src.s & (1 << SIGN))) return s_is_minus;
    if(s & (1 << INF)) {
      if(src.s & (1 << INF)) throw "compair INF";
      return s_is_minus;
    }
    if(src.s & (1 << INF)) return ! s_is_minus;
    if(m && src.m) {
      if(e < src.e) return ! s_is_minus;
      if(e == src.e) return s_is_minus ? src.m < m : m < src.m;
      return s_is_minus;
    }
    return !m ? (! src.m ? m != src.m : ! (src.s & (1 << SIGN))) : s_is_minus;
  }
  inline bool             operator <= (const SimpleFloat<T,W,bits,U>& src) const {
    return *this < src || *this == src;
  }
  inline bool             operator >  (const SimpleFloat<T,W,bits,U>& src) const {
    return ! (*this <= src);
  }
  inline bool             operator >= (const SimpleFloat<T,W,bits,U>& src) const {
    return ! (*this < src);
  }
  inline bool             operator !  () const {
    return ! m && isfinite(*this);
  }
  inline                  operator bool () const {
    return ! (!*this);
  }
  inline                  operator int  () const {
    return int(this->operator T());
  }
  inline                  operator T    () const {
    auto deci(*this);
    if(deci.s & (1 << INF)) throw "Inf to convert int";
    if(deci.s & (1 << NaN)) throw "NaN to convert int";
    if(! deci.m) return T(int(0));
    if(U(bits) <= deci.e || (uzero() < deci.e && (deci.m << int(deci.e)) >> int(deci.e) != deci.m))
      throw "Overflow to convert int.";
    if(deci.e <= - U(bits)) return T(int(0));
    if(     deci.e < uzero()) deci.m >>= - int(deci.e);
    else if(uzero() < deci.e) deci.m <<=   int(deci.e);
    return s & (1 << SIGN) ? - deci.m : deci.m;
  }
  inline                  operator SimpleFloat<T,W,bits,U> () const {
    return *this;
  }
  // XXX: absfloor, absceil implementation.
  inline SimpleFloat<T,W,bits,U>  floor() const {
    if(uzero() <= e) return *this;
    if(e <= - U(bits)) return zero();
    auto deci(*this);
    deci.m >>= - int(deci.e);
    deci.m <<= - int(deci.e);
    return deci;
  }
  inline SimpleFloat<T,W,bits,U>  ceil() const {
    const auto fl(this->floor());
    if(*this - fl) {
      auto pmone(one());
      pmone.s |= s & (1 << SIGN);
      return fl + pmone;
    }
    return fl;
  }

  inline SimpleFloat<T,W,bits,U>  abs()  const {
    auto work(*this);
    work.s &= ~ (1 << SIGN);
    return work;
  }
         SimpleFloat<T,W,bits,U>  log()  const;
         SimpleFloat<T,W,bits,U>  exp()  const;
         SimpleFloat<T,W,bits,U>  sin()  const;
         SimpleFloat<T,W,bits,U>  cos()  const;
         SimpleFloat<T,W,bits,U>  atan() const;
  inline SimpleFloat<T,W,bits,U>  sqrt() const;
  
  unsigned char s;
  typedef enum {
    INF = 0,
    NaN = 1,
    SIGN = 2,
    DWRK = 3
  } state_t;
  T m;
  U e;
  const U& uzero() const {
    const static U vuzero(0);
    return vuzero;
  }
  const SimpleFloat<T,W,bits,U>& zero()   const {
    const static SimpleFloat<T,W,bits,U> vzero(0);
    return vzero;
  }
  const SimpleFloat<T,W,bits,U>& one()    const {
    const static SimpleFloat<T,W,bits,U> vone(1);
    return vone;
  }
  const SimpleFloat<T,W,bits,U>& two()    const {
    const static auto vtwo(one() << U(1));
    return vtwo;
  }
  const SimpleFloat<T,W,bits,U>& pi()     const {
    const static auto vpi(quatpi() << U(2));
    return vpi;
  }
  const SimpleFloat<T,W,bits,U>& halfpi() const {
    const static auto vhalfpi(quatpi() << U(1));
    return vhalfpi;
  }
  const SimpleFloat<T,W,bits,U>& quatpi() const {
    const static auto vquatpi(one().atan());
    return vquatpi;
  }
  const SimpleFloat<T,W,bits,U>& twopi()  const {
    const static auto vtwopi(quatpi() << U(3));
    return vtwopi;
  }
  const SimpleFloat<T,W,bits,U>& sqrt2()  const {
    const static auto vsqrt2((one() << U(1)).sqrt());
    return vsqrt2;
  }
private:
  template <typename V> inline U normalize(V& src) const {
    V   bt(int(1));
    int b(0);
    int tb(0);
    for( ; bt; tb ++) {
      if(src & bt) b = tb;
      bt <<= 1;
    }
    const auto shift(tb - b - 1);
    assert(0 <= shift);
    if(shift) src <<= shift;
    return - U(shift);
  }
  inline SimpleFloat<T,W,bits,U>& ensureFlag() {
    if(s & (1 << INF))
      s &= ~ (1 << DWRK);
    else if(! m || (s & (1 << DWRK))) {
      e ^= e;
      m ^= m;
      s &= ~ ((1 << DWRK) | (1 << INF));
    }
    return * this;
  }
  inline unsigned char safeAdd(U& dst, const U& src) {
    const auto dst0(dst);
    dst += src;
    if((dst0 > uzero() && src > uzero() && dst <= uzero()) ||
       (dst0 < uzero() && src < uzero() && dst >= uzero()))
      return 1 << (dst0 < uzero() ? DWRK : INF);
    return 0;
  }
  inline char residue2() const {
    if(uzero() < e || U(bits) <= - e) return 0;
    if(! e) return char(int(m) & 1);
    return char(int(m >> - int(e)) & 1);
  }

  // XXX: these are NOT threadsafe on first call.
  const vector<SimpleFloat<T,W,bits,U> >& exparray()    const;
  const vector<SimpleFloat<T,W,bits,U> >& invexparray() const;
/*
friend:
  ostream&    operator << (ostream& os, const SimpleFloat<T,W,bits,U>& v);
  istream&    operator >> (istream& is, SimpleFloat<T,W,bits,U>& v);
*/
};

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
  if(s & (1 << SIGN))
    return - (- *this).atan();
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
    const auto v(five * *this / (four + (*this << U(1))) - half);
    assert(v < *this);
    return atanhalf + v.atan();
  }
  // N.B.
  //    in u = v case,
  //  2 atan(u) = atan(2 * u / (1 - u * u))
  //  2 atan(u) = atan(2 * u / (1 - u) / (1 + u))
  //            = atan(2 / (1 / u - u))
  //    in 2Y := 1 / u - u case,
  //            = atan(1 / Y),
  //  u^2 + 2Yu - 1 == 0, u = - Y \pm sqrt(Y^2 + 1)
  const auto Y(one() / (*this));
  const auto u((Y * Y + one()).sqrt() - Y);
  assert(- *this < u && u < *this);
  return u.atan() << U(1);
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

template <typename T, typename W, int bits, typename U> ostream& operator << (ostream& os, const SimpleFloat<T,W,bits,U>& v) {
  static const U uzero(int(0));
  if(isnan(v))
    return os << "NaN ";
  if(isinf(v))
    return os << (const char*)(v.s & (1 << v.SIGN) ? "-" : "") << "Inf ";
  return os << (const char*)(v.s & (1 << v.SIGN) ? "-" : "") << std::hex << T(v.m) << "*2^" << (const char*)(v.e < uzero ? "-" : "") << (v.e < uzero ? U(- v.e) : v.e) << " " << std::dec;
}

template <typename T, typename W, int bits, typename U> istream& operator >> (istream& is, SimpleFloat<T,W,bits,U>& v) {
  const static SimpleFloat<T,W,bits,U> two(2);
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
        e <<= U(int(4));
        e  += SimpleFloat<T,W,bits,U>(int(buf - '0'));
      } else {
        v <<= U(int(4));
        v  += SimpleFloat<T,W,bits,U>(int(buf - '0'));
      }
      fsign = true;
      break;
    case 'a': case'b': case 'c': case 'd': case 'e': case 'f':
      if(mode) {
        e <<= U(int(4));
        e  += SimpleFloat<T,W,bits,U>(int(buf - 'a' + 10));
      } else {
        v <<= U(int(4));
        v  += SimpleFloat<T,W,bits,U>(int(buf - 'a' + 10));
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
  inline Complex() { ; }
  inline Complex(const T& real, const T& imag = T(int(0))) {
    _real = real; _imag = imag;
  }
  inline Complex(const Complex<T>& s) { *this = s; }
  inline Complex(Complex<T>&& s) { *this = s; }
  inline Complex(T&& real) {
    const static T zero(0);
    _real = move(real);
    _imag = zero;
  }
  inline Complex(T&& real, T&& imag) {
    _real = move(real);
    _imag = move(imag);
  }
  inline ~Complex() { ; }

  inline Complex<T>  operator ~  ()                    const {
    return Complex<T>(  _real, - _imag);
  }
  inline Complex<T>  operator -  ()                    const {
    return Complex<T>(- _real, - _imag);
  }
  inline Complex<T>  operator +  (const Complex<T>& s) const {
    auto result(*this);
    return result += s;
  }
  inline Complex<T>& operator += (const Complex<T>& s) {
    _real += s._real;
    _imag += s._imag;
    return *this;
  }
  inline Complex<T>  operator -  (const Complex<T>& s) const {
    auto result(*this);
    return result -= s;
  }
  inline Complex<T>& operator -= (const Complex<T>& s) {
    _real -= s._real;
    _imag -= s._imag;
    return *this;
  }
  inline Complex<T>  operator *  (const T& s)          const {
    auto result(*this);
    return result *= s;
  }
  inline Complex<T>& operator *= (const T& s) {
    _real *= s;
    _imag *= s;
    return *this;
  }
  inline Complex<T>  operator *  (const Complex<T>& s) const {
    return Complex<T>(_real * s._real - _imag * s._imag,
                      _real * s._imag + _imag * s._real);
  }
  inline Complex<T>& operator *= (const Complex<T>& s) {
    return (*this) = (*this) * s;
  }
  inline Complex<T>  operator /  (const T& s)          const {
    auto result(*this);
    return result /= s;
  }
  inline Complex<T>& operator /= (const T& s) {
    _real /= s;
    _imag /= s;
    return *this;
  }
  inline Complex<T>  operator /  (const Complex<T>& s) const {
    return (*this * (~ s)) / (s._real * s._real + s._imag * s._imag);
  }
  inline Complex<T>& operator /= (const Complex<T>& s) {
    return *this = *this / s;
  }
  inline bool        operator == (const Complex<T>& s) const {
    return !(*this != s);
  }
  inline bool        operator != (const Complex<T>& s) const {
    return (_real != s._real) || (_imag != s._imag);
  }
  inline bool        operator !  ()                    const {
    return !_real && !_imag;
  }
  inline Complex<T>  operator &  (const Complex<T>& s) const {
    auto result(*this);
    return result &= s;
  }
  inline Complex<T>& operator &= (const Complex<T>& s) {
    _real &= s._real;
    _imag &= s._imag;
    return *this;
  }
  inline Complex<T>  operator |  (const Complex<T>& s) const {
    auto result(*this);
    return result |= s;
  }
  inline Complex<T>& operator |= (const Complex<T>& s) {
    _real |= s._real;
    _imag |= s._imag;
    return *this;
  }
  inline Complex<T>  operator ^  (const Complex<T>& s) const {
    auto result(*this);
    return result ^= s;
  }
  inline Complex<T>& operator ^= (const Complex<T>& s) {
    _real ^= s._real;
    _imag ^= s._imag;
    return *this;
  }
  inline bool        operator && (const Complex<T>& s) const {
    return *this && s;
  }
  inline bool        operator || (const Complex<T>& s) const {
    return *this || s;
  }
  inline Complex<T>& operator =  (const Complex<T>& s) {
    _real = s._real;
    _imag = s._imag;
    return *this;
  }
  inline Complex<T>& operator =  (Complex<T>&& s) {
    _real = move(s._real);
    _imag = move(s._imag);
    return *this;
  }
  inline T&          operator [] (const size_t& i) {
    assert(0 <= i && i < 2);
    if(i) return _imag;
    return _real;
  }
  inline             operator bool () const {
    return ! (! *this);
  }
  inline             operator T    () const {
    return this->_real;
  }
  
  const Complex<T>& i() const {
    const static auto I(Complex<T>(T(int(0)), T(int(1))));
    return I;
  }
  inline T  abs() const {
    return sqrt(_real * _real + _imag * _imag);
  }
  inline T  arg() const {
    return atan2(_imag, _real);
  }
  inline T& real() {
    return _real;
  }
  inline T& imag() {
    return _imag;
  }
  inline const T& real() const {
    return _real;
  }
  inline const T& imag() const {
    return _imag;
  }
  T _real;
  T _imag;
};

template <typename T> ostream& operator << (ostream& os, const Complex<T>& v) {
  return os << v.real() << "+i" << v.imag();
}

template <typename T> istream& operator >> (istream& is, Complex<T>& v) {
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
  v._imag = T(int(0));
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
  return exp(log(s) * Complex<T>(T(int(1)) / T(int(2))));
}

template <typename T> static inline Complex<T> csin(const Complex<T>& s) {
  return (exp(Complex<T>(T(int(0)), s)) - exp(Complex<T>(T(int(0)), - s))) / Complex<T>(T(int(0)), T(int(2)));
}

template <typename T> static inline Complex<T> ccos(const Complex<T>& s) {
  return (exp(Complex<T>(T(int(0)), s)) + exp(Complex<T>(T(int(0)), - s))) / T(int(2));
}

template <typename T> static inline Complex<T> ctan(const Complex<T>& s) {
  return csin(s) / ccos(s);
}

template <typename T> static inline Complex<T> ccsc(const Complex<T>& s) {
  return Complex<T>(T(int(1))) / csin(s);
}

template <typename T> static inline Complex<T> csec(const Complex<T>& s) {
  return Complex<T>(T(int(1))) / ccos(s);
}

template <typename T> static inline T ccot(const T& s) {
  return Complex<T>(T(int(1))) / ctan(s);
}

template <typename T> using complex = Complex<T>;

#if !defined(_FLOAT_BITS_)
  #include <cmath>
  using namespace std;
  typedef uint64_t myuint;
  typedef int64_t  myint;
  // XXX:
  // typedef long double myfloat;
  typedef double myfloat;
#elif _FLOAT_BITS_ == 8
  typedef uint8_t myuint;
  typedef int8_t  myint;
  typedef SimpleFloat<myuint, uint16_t, 8, myint> myfloat;
#elif _FLOAT_BITS_ == 16
  typedef uint16_t myuint;
  typedef int16_t  myint;
  typedef SimpleFloat<myuint, uint32_t, 16, myint> myfloat;
#elif _FLOAT_BITS_ == 32
  typedef uint32_t myuint;
  typedef int32_t  myint;
  typedef SimpleFloat<myuint, uint64_t, 32, myint> myfloat;
#elif _FLOAT_BITS_ == 64
  typedef uint64_t myuint;
  typedef int64_t  myint;
  typedef SimpleFloat<myuint, unsigned __int128, 64, myint> myfloat;
#elif _FLOAT_BITS_ == 128
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef Signed<uint128_t, 128> int128_t;
  typedef uint128_t myuint;
  typedef int128_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 128>, 128, myint> myfloat;
#elif _FLOAT_BITS_ == 256
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef DUInt<uint128_t, 128> uint256_t;
  typedef Signed<uint256_t, 256> int256_t;
  typedef uint256_t myuint;
  typedef int256_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 256>, 256, myint> myfloat;
#elif _FLOAT_BITS_ == 512
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef DUInt<uint128_t, 128> uint256_t;
  typedef DUInt<uint256_t, 256> uint512_t;
  typedef Signed<uint512_t, 512> int512_t;
  typedef uint512_t myuint;
  typedef int512_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 512>, 512, myint> myfloat;
#elif _FLOAT_BITS_ == 1024
  typedef DUInt<uint64_t, 64> uint128_t;
  typedef DUInt<uint128_t, 128> uint256_t;
  typedef DUInt<uint256_t, 256> uint512_t;
  typedef DUInt<uint512_t, 512> uint1024_t;
  typedef Signed<uint1024_t, 1024> int1024_t;
  typedef uint1024_t myuint;
  typedef int1024_t  myint;
  typedef SimpleFloat<myuint, DUInt<myuint, 1024>, 1024, myint> myfloat;
#else
# error cannot handle float
#endif


// N.B. start simplelin.
template <typename T> class SimpleVector {
public:
  inline SimpleVector() { ; }
  inline SimpleVector(const int& size) {
    assert(0 <= size);
    this->entity.resize(size);
  }
  inline SimpleVector(const SimpleVector<T>& other) { *this = other; }
  inline SimpleVector(SimpleVector<T>&& other) { *this = other; }
  inline ~SimpleVector() { ; }
  
  inline       SimpleVector<T>  operator -  () const {
    SimpleVector<T> res(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      res.entity[i] = - entity[i];
    return res;
  }
  inline       SimpleVector<T>  operator +  (const SimpleVector<T>& other) const {
    auto res(*this);
    return res += other;
  }
  inline       SimpleVector<T>& operator += (const SimpleVector<T>& other) {
    assert(entity.size() == other.entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] += other.entity[i];
    return *this;
  }
  inline       SimpleVector<T>  operator -  (const SimpleVector<T>& other) const {
    auto res(*this);
    return res -= other;
  }
  inline       SimpleVector<T>& operator -= (const SimpleVector<T>& other) {
    return *this += - other;
  }
  inline       SimpleVector<T>  operator *  (const T& other) const {
    auto res(*this);
    return res *= other;
  }
  inline       SimpleVector<T>& operator *= (const T& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] *= other;
    return *this;
  }
  inline       SimpleVector<T>  operator /  (const T& other) const {
    auto res(*this);
    return res /= other;
  }
  inline       SimpleVector<T>& operator /= (const T& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] /= other;
    return *this;
  }
  inline       SimpleVector<T>& operator =  (const SimpleVector<T>& other) {
    entity = other.entity;
    return *this;
  }
  inline       SimpleVector<T>& operator =  (SimpleVector<T>&& other) {
    entity = move(other.entity);
    return *this;
  }
  inline       bool             operator == (const SimpleVector<T>& other) const {
    return ! (*this != other);
  }
  inline       bool             operator != (const SimpleVector<T>& other) const {
    if(entity.size() != other.entity.size()) return true;
    for(int i = 0; i < entity.size(); i ++)
      if(entity[i] != other.entity[i]) return true;
    return false;
  }
  template <typename U> inline T dot(const SimpleVector<U>& other) const {
    assert(entity.size() == other.entity.size());
    SimpleVector<T> work(other.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      work[i] = entity[i] * other.entity[i];
    auto res(work[0]);
    for(int i = 1; i < entity.size(); i ++)
      res += work[i];
    return res;
  }
  inline       T&               operator [] (const int& idx) {
    assert(0 <= idx && idx < entity.size());
    return entity[idx];
  }
  inline const T&               operator [] (const int& idx) const {
    assert(0 <= idx && idx < entity.size());
    return entity[idx];
  }
  template <typename U> inline SimpleVector<U> real() const {
    SimpleVector<U> result(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      result.entity[i] = U(entity[i].real());
    return result;
  }
  template <typename U> inline SimpleVector<U> imag() const {
    SimpleVector<U> result(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      result.entity[i] = U(entity[i].imag());
    return result;
  }
  template <typename U> inline SimpleVector<U> cast() const {
    SimpleVector<U> result(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < entity.size(); i ++)
      result.entity[i] = U(entity[i]);
    return result;
  }
  inline const int size() const {
    return entity.size();
  }
  inline       void resize(const int& size) {
    assert(0 <= size);
    entity.resize(size);
    return;
  }
  inline       SimpleVector<T>  subVector(const int& i, const int& s) const {
    assert(0 <= s && 0 <= i && i + s <= size());
    SimpleVector<T> res(s);
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int ii = i; ii < i + s; ii ++)
      res[ii - i] = (*this)[ii];
    return res;
  }
  inline       SimpleVector<T>& setVector(const int& i, const SimpleVector<T>& d) {
    assert(0 <= i && i + d.size() <= size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int ii = i; ii < i + d.size(); ii ++)
      (*this)[ii] = d[ii - i];
    return *this;
  }
  inline       SimpleVector<T>& O(const T& r = T(int(0))) {
    return I(r);
  }
  inline       SimpleVector<T>& I(const T& r = T(int(1))) {
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < size(); i ++)
      (*this)[i] = r;
    return *this;
  }
  inline       SimpleVector<T>& ek(const int& i, const T& r = T(int(1))) {
    static const T zero(0);
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int ii = 0; ii < size(); ii ++)
      (*this)[ii] = ii == i ? r : zero;
    return *this;
  }
  inline       SimpleVector<T>  reverse() {
    SimpleVector<T> res(entity.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
    for(int i = 0; i < res.size(); i ++)
      res[i] = entity[entity.size() - 1 - i];
    return res;
  }
  
  vector<T> entity;
};

template <typename T> ostream& operator << (ostream& os, const SimpleVector<T>& v) {
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

template <typename T> istream& operator >> (istream& is, SimpleVector<T>& v) {
  int s;
  is >> s;
  if(s <= 0) return is;
  v.resize(s);
  int i(0);
  for( ; i < v.size() && ! is.eof() && ! is.bad(); ) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == ':' || c == ',' || c == '[' || c == '\n') continue;
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
      v[i] = T(int(0));
  }
  return is;
}


template <typename T> class SimpleMatrix {
public:
  inline SimpleMatrix() { ecols = 0; }
  inline SimpleMatrix(const int& rows, const int& cols) {
    assert(0 <= rows && 0 <= cols);
    entity.resize(rows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i].resize(cols);
    ecols = cols;
  }
  inline SimpleMatrix(const SimpleMatrix<T>& other) { *this = other; }
  inline SimpleMatrix(SimpleMatrix<T>&& other) { *this = other; }
  inline ~SimpleMatrix() { ; }
  
  inline       SimpleMatrix<T>  operator -  () const {
    SimpleMatrix<T> res(entity.size(), ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      res.entity[i] = - entity[i];
    return res;
  }
  inline       SimpleMatrix<T>  operator +  (const SimpleMatrix<T>& other) const {
    auto res(*this);
    return res += other;
  }
  inline       SimpleMatrix<T>& operator += (const SimpleMatrix<T>& other) {
    assert(entity.size() == other.entity.size() && ecols == other.ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] += other.entity[i];
    return *this;
  }
  inline       SimpleMatrix<T>  operator -  (const SimpleMatrix<T>& other) const {
    auto res(*this);
    return res -= other;
  }
  inline       SimpleMatrix<T>& operator -= (const SimpleMatrix<T>& other) {
    return *this += - other;
  }
  inline       SimpleMatrix<T>  operator *  (const T& other) const {
    auto res(*this);
    return res *= other;
  }
  inline       SimpleMatrix<T>& operator *= (const T& other) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] *= other;
    return *this;
  }
  inline       SimpleMatrix<T>  operator *  (const SimpleMatrix<T>& other) const {
    assert(ecols == other.entity.size() && entity.size() && other.entity.size());
    auto            derived(other.transpose());
    SimpleMatrix<T> res(entity.size(), other.ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++) {
            SimpleVector<T>& resi(res.entity[i]);
      const SimpleVector<T>& ei(entity[i]);
      for(int j = 0; j < other.ecols; j ++)
        resi[j] = ei.dot(derived.entity[j]);
    }
    return res;
  }
  inline       SimpleMatrix<T>& operator *= (const SimpleMatrix<T>& other) {
    return *this = *this * other;
  }
  inline       SimpleVector<T>  operator *  (const SimpleVector<T>& other) const {
    assert(ecols == other.size());
    SimpleVector<T> res(entity.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      res[i] = entity[i].dot(other);
    return res;
  }
  inline       SimpleMatrix<T>  operator /  (const T& other) const {
    auto res(*this);
    return res /= other;
  }
  inline       SimpleMatrix<T>& operator /= (const T& other) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i] /= other;
    return *this;
  }
  inline       SimpleMatrix<T>& operator =  (const SimpleMatrix<T>& other) {
    ecols  = other.ecols;
    entity = other.entity;
    return *this;
  }
  inline       SimpleMatrix<T>& operator =  (SimpleMatrix<T>&& other) {
    ecols  = move(other.ecols);
    entity = move(other.entity);
    return *this;
  }
  inline       bool             operator == (const SimpleMatrix<T>& other) const {
    return ! (*this != other);
  }
  inline       bool             operator != (const SimpleMatrix<T>& other) const {
    if(entity.size() != other.entity.size() || ecols != other.ecols)
      return true;
    for(int i = 0; i < entity.size(); i ++)
      if(entity[i] != other.entity[i]) return true;
    return false;
  }
  inline       T&               operator () (const int& y, const int& x) {
    assert(0 <= y && y < entity.size());
    return entity[y][x];
  }
  inline const T&               operator () (const int& y, const int& x) const {
    assert(0 <= y && y < entity.size());
    return entity[y][x];
  }
  inline       SimpleVector<T>& row(const int& y) {
    assert(0 <= y && y < entity.size());
    return entity[y];
  }
  inline const SimpleVector<T>& row(const int& y) const {
    assert(0 <= y && y < entity.size());
    return entity[y];
  }
  inline const SimpleVector<T>  col(const int& x) const {
    assert(0 <= entity.size() && 0 <= x && x < ecols);
    SimpleVector<T> res(entity.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      res[i] = entity[i][x];
    return res;
  }
  inline       void             setCol(const int& x, const SimpleVector<T>& other) {
    assert(0 <= x && x < ecols && other.size() == entity.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i][x] = other[i];
    return;
  }
  // N.B. transpose : exhaust of the resource, so Eigen library handles better.
  inline       SimpleMatrix<T>  transpose() const {
    SimpleMatrix<T> res(ecols, entity.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < ecols; i ++) {
      SimpleVector<T>& resi(res.entity[i]);
      for(int j = 0; j < entity.size(); j ++)
        resi[j] = entity[j][i];
    }
    return res;
  }
  inline       SimpleMatrix<T>  subMatrix(const int& y, const int& x, const int& h, const int& w) const {
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
  inline       SimpleMatrix<T>& setMatrix(const int& y, const int& x, const SimpleMatrix<T>& d) {
    assert(0 <= y && y + d.rows() <= rows() && 0 <= x && x + d.cols() <= cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = y; i < y + d.rows(); i ++)
      for(int j = x; j < x + d.cols(); j ++)
        (*this)(i, j) = d(i - y, j - x);
    return *this;
  }
  inline       SimpleMatrix<T>& O(const T& r = T(int(0))) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < rows(); i ++)
      for(int j = 0; j < cols(); j ++)
        (*this)(i, j) = r;
    return *this;
  }
  inline       SimpleMatrix<T>& I(const T& r = T(int(1))) {
    const static T zero(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < rows(); i ++)
      for(int j = 0; j < cols(); j ++)
        (*this)(i, j) = (i == j ? r : zero);
    return *this;
  }
  inline       T                determinant(const bool& nonzero = false) const;
  inline       SimpleMatrix<T>  inverse() const {
    // XXX: extremely slow implementation.
    SimpleMatrix<T> result(entity.size(), ecols);
    result.I();
    for(int i = 0; i < result.cols(); i ++)
      result.setCol(i, solve(result.col(i)));
    return result;
  }
  inline       SimpleVector<T>  solve(SimpleVector<T> other) const;
  inline       SimpleVector<T>  projectionPt(const SimpleVector<T>& other) const;
  inline       SimpleMatrix<T>& fillP(const vector<int>& idx);
  inline       SimpleMatrix<T>  QR() const;
  inline       SimpleMatrix<T>  SVDleft1d() const;
  inline       pair<SimpleMatrix<T>, SimpleMatrix<T> > SVD1d() const;
  inline       SimpleMatrix<T> SVD() const;
  inline       pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SVD(const SimpleMatrix<T>& src) const;
  inline       SimpleVector<T>  zeroFix(const SimpleMatrix<T>& A, vector<pair<T, int> > fidx);
  inline       SimpleVector<T>  inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const;
  template <typename U> inline SimpleMatrix<U> real() const {
    assert(0 < entity.size() && 0 < ecols);
    SimpleMatrix<U> res(entity.size(), ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      for(int j = 0; j < ecols; j ++)
        res(i, j) = U(entity[i][j].real());
    return res;
  }
  template <typename U> inline SimpleMatrix<U> imag() const {
    assert(0 < entity.size() && 0 < ecols);
    SimpleMatrix<U> res(entity.size(), ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      for(int j = 0; j < ecols; j ++)
        res(i, j) = U(entity[i][j].imag());
    return res;
  }
  template <typename U> inline SimpleMatrix<U> cast() const {
    assert(0 < entity.size() && 0 < ecols);
    SimpleMatrix<U> res(entity.size(), ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      for(int j = 0; j < ecols; j ++)
        res(i, j) = U(entity[i][j]);
    return res;
  }
  inline const int rows() const {
    return entity.size();
  }
  inline const int cols() const {
    return ecols;
  }
  inline       void resize(const int& rows, const int& cols) {
    assert(0 <= rows && 0 <= cols);
    ecols = cols;
    entity.resize(rows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < entity.size(); i ++)
      entity[i].resize(ecols);
    return;
  }
  myfloat      epsilon() const {
#if defined(_FLOAT_BITS_)
    // N.B. conservative.
    static const auto eps(sqrt(myfloat(int(1)) >> myint(_FLOAT_BITS_ - 1)));
    // static const auto eps(myfloat(int(1)) >> myint(_FLOAT_BITS_ - 1));
#else
    // N.B. conservative.
    static const auto eps(sqrt(std::numeric_limits<myfloat>::epsilon()));
    // static const auto eps(std::numeric_limits<myfloat>::epsilon());
#endif
    return eps;
  }

  // friend ostream& operator << (ostream& os, const SimpleVector<T>& v);
  // friend istream& operator >> (istream& os, SimpleVector<T>& v);

private:
  // this isn't better idea for faster calculations.
  vector<SimpleVector<T> > entity;
  int ecols;
};

template <typename T> inline T SimpleMatrix<T>::determinant(const bool& nonzero) const {
  assert(0 <= entity.size() && 0 <= ecols && entity.size() == ecols);
  T det(1);
  auto work(*this);
  for(int i = 0; i < entity.size(); i ++) {
    int xchg = i;
    for(int j = i + 1; j < entity.size(); j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    swap(work.entity[i], work.entity[xchg]);
    const auto& ei(work.entity[i]);
    const auto& eii(ei[i]);
    if(! nonzero || ! i || pow(abs(det), T(int(1)) / T(int(i))) * epsilon() <= abs(eii))
      det *= eii;
    if(ei.dot(ei) * epsilon() < eii * eii) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = i + 1; j < entity.size(); j ++) {
        const auto ratio(work.entity[j][i] / eii);
        work.entity[j] -= ei * ratio;
      }
    }
  }
  return det;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::solve(SimpleVector<T> other) const {
  if(! (0 <= entity.size() && 0 <= ecols && entity.size() == ecols && entity.size() == other.size()) ) throw "SimpleMatrix<T>::Solve error";
  auto work(*this);
  for(int i = 0; i < entity.size(); i ++) {
    int xchg = i;
    for(int j = i + 1; j < entity.size(); j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    swap(work.entity[i], work.entity[xchg]);
    swap(other[i], other[xchg]);
    const auto& ei(work.entity[i]);
    const auto& eii(ei[i]);
    if(ei.dot(ei) * epsilon() < eii * eii) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = i + 1; j < entity.size(); j ++) {
        const auto ratio(work.entity[j][i] / eii);
        work.entity[j] -= ei       * ratio;
        other[j]       -= other[i] * ratio;
      }
    }
  }
  for(int i = entity.size() - 1; 0 <= i; i --) {
    if(work.entity[i][i] == T(int(0))) continue;
    const auto buf(other[i] / work.entity[i][i]);
    if(!isfinite(buf) || isnan(buf)) {
      // assert(!isfinite(work.entity[i][i] / other[i]) || isnan(work.entity[i][i] / other[i]));
      continue;
    }
    other[i] = buf;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = i - 1; 0 <= j; j --)
      other[j] -= other[i] * work.entity[j][i];
  }
  return other;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::projectionPt(const SimpleVector<T>& other) const {
  assert(0 < entity.size() && 0 < ecols && ecols == other.size());
  // also needs class or this->transpose() * (*this) == I assertion is needed.
  SimpleMatrix<T> work(entity.size(), ecols);
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
    res[i] = T(int(0));
    for(int j = 0; j < entity.size(); j ++)
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
    if(n2 <= epsilon()) continue;
    assert(0 <= idx[ii] && idx[ii] < this->rows());
    this->row(idx[ii ++]) = ek / sqrt(n2);
  }
  assert(idx.size() <= ii);
  return *this;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::QR() const {
  const auto norm2(norm2M(*this));
  if(! isfinite(norm2)) return *this;
  SimpleMatrix<T> Q(min(this->rows(), this->cols()), this->rows());
  Q.O();
  vector<int> residue;
  residue.reserve(Q.rows());
  for(int i = 0; i < Q.rows(); i ++) {
    const auto Atrowi(this->col(i));
    const auto work(Atrowi - Q.projectionPt(Atrowi));
    const auto n2(work.dot(work));
    if(n2 <= norm2 * epsilon())
      residue.emplace_back(i);
    else
      Q.row(i) = work / sqrt(n2);
  }
  return Q.fillP(residue);
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::SVDleft1d() const {
  // N.B. A = QR, (S - lambda I)x = 0 <=> R^t Q^t U = R^-1 Q^t U Lambda
  //        <=> R^t Q^t U Lambda' = R^-1 Q^t U Lambda'^(- 1)
  //      A := R^t, B := Q^t U, C := Lambda'
  //          (A + A^-t)*B*(C + C^-1) - A^-t B C - A B C^-1 = 2 A B C
  //        <=> (A + A^-t) * B * (C + C^-1) = (2I + 2A^-tA^-1) * ABC
  //        <=> B = A^-1 (2I + 2A^-(t+1))^-1 (A + A^-t) * B *
  //                (C + C^-1) * C^-1
  // N.B. since S is symmetric, singular value on SS^t = QRR^tQ^t is
  //      same square root as singular value on R.
  const auto S(*this * transpose());
  const auto Qt(S.QR());
  const auto A((Qt * S).transpose());
  const auto A1t(A * A.transpose());
        auto Left(A.inverse() * (SimpleMatrix<T>(A.rows(), A.cols()).I(T(int(2))) + A1t.inverse() * T(int(2))).inverse() * (A + A.transpose().inverse()));
        auto Right(SimpleMatrix<T>(Left.rows(), Left.cols()).O());
  for(int i = 0; i < Right.rows(); i ++)
    Right(i, i) = A(i, i) + T(int(1));
  Left  /= sqrt(norm2M(Left));
  Right /= sqrt(norm2M(Right));
  // N.B. now we have B = Left * B * Right.
  static const int p(ceil(sqrt(- log(epsilon()))));
  for(int i = 0; i < p; i ++) {
    Left  *= Left;
    Right *= Right;
  }
  return (Left * Right).QR() * Qt;
}

template <typename T> inline pair<SimpleMatrix<T>, SimpleMatrix<T> > SimpleMatrix<T>::SVD1d() const {
  if(this->rows() < this->cols()) {
    auto R(this->transpose().SVDleft1d().transpose());
    return make_pair(((*this) * R).QR(), move(R));
  }
  auto L(this->SVDleft1d());
  return make_pair(move(L), (L * (*this)).transpose().QR().transpose());
}

// XXX: O(n^4) over all, we need O(n^3) methods they make SVD1d as SVDnd.
template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::SVD() const {
  auto sym((*this) * this->transpose());
  assert(sym.rows() == sym.cols());
  auto res(sym);
  res.I();
  for(int i = 0; i <= sym.rows() + 1; i ++) {
    auto svd(sym.SVD1d());
    sym = (svd.first * sym * svd.second).transpose();
    if(i & 1)
      res = svd.second.transpose() * res;
    else
      res = svd.first * res;
  }
  return res;
}

template <typename T> inline pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SimpleMatrix<T>::SVD(const SimpleMatrix<T>& src) const {
  const auto norm2(max(norm2M(*this), norm2M(src)));
  if(! isfinite(norm2)) return *this;
  // refered from : https://en.wikipedia.org/wiki/Generalized_singular_value_decomposition .
  assert(this->cols() == src.cols());
  SimpleMatrix<T> C(this->rows() + src.rows(), this->cols());
  C.setMatrix(0, 0, *this);
  C.setMatrix(this->rows(), 0, src);
  const auto P(C.SVD());
  const auto Qt(C.transpose().SVD().transpose());
  const auto D(P.first * C * Qt.transpose());
  SimpleMatrix<T> P1(this->rows(), this->cols());
  SimpleMatrix<T> P2(src.rows(), this->cols());
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
  auto Wt(P1.transpose().SVD().transpose());
  auto U2(Wt * P2.transpose());
  vector<int> fill;
  fill.reserve(U2.rows());
  for(int i = 0; i < U2.rows(); i ++) {
    const auto n2(U2.row(i).dot(U2.row(i)));
    if(n2 <= epsilon())
      fill.emplace_back(i);
    else
      U2.row(i) /= sqrt(n2);
  }
  return make_pair(make_pair(move(U1), move(U2.fillP(fill))), (Wt * D).transpose().QR() * Qt);
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::zeroFix(const SimpleMatrix<T>& A, vector<pair<T, int> > fidx) {
  // N.B. we now have |[A -bb] [x t]| <= 1 condition.
  // N.B. there's no difference |[A - bb] [x t]|^2 <= 1 condition in this.
  //      but not with mixed condition.
  const auto R((*this) * A);
  SimpleVector<T> one(this->cols());
  one.I(T(int(1)));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < fidx.size(); i ++) {
    one[fidx[i].second] = - fidx[i].first;
    fidx[i].first = - T(int(1));
  }
  // we now have: Q [R [x t] ] <= {0, 1}^m cond.
  auto Pb(*this);
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
    if(n2 <= epsilon())
      continue;
    if(T(int(0)) < fidx[idx].first &&
       fidx[idx].first < sqrt(one.dot(one)) * epsilon()) {
      cerr << "linearInvariant: P matrix is orthogonal to 1 vector." << endl;
      *this = Pb;
      break;
    }
    Pb = *this;
    // N.B. O(mn) can be writed into O(lg m + lg n) in many core cond.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < this->cols(); j ++)
      this->setCol(j, this->col(j) - orth * this->col(j).dot(orth) / n2);
    if(T(int(0)) < fidx[idx].first) {
      const auto rfidxsz(fidx.size());
      fidx.resize(0);
      fidx.reserve(this->cols());
      const auto on(projectionPt(one));
      for(int j = 0; j < this->cols(); j ++)
        fidx.emplace_back(make_pair(abs(on[j]), i));
      sort(fidx.begin(), fidx.end());
      i -= rfidxsz - fidx.size();
    }
    i ++;
  }
  // N.B. now we have fix indexes that to be P R [x 1] * t == 0.
  return R.solve((*this) * one);
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const {
  assert(this->rows() == bl.size() && this->rows() == bu.size() &&
         0 < this->cols() && 0 < this->rows() && this->cols() < this->rows());
  // |(2 / bu) A x - 1 - bl / bu| <= |1 - bl / bu|
  // <=> with (-A, -bu, -bl), |bl| <= |bu|, |(2 / bu) A x - 2| <= 2(1 - bl / bu)
  auto bU(bu);
  auto bL(bl);
  SimpleMatrix<T> A(*this);
  vector<pair<T, int> > fidx;
  for(int i = 0; i < bU.size(); i ++) {
    if(abs(bu[i]) < abs(bl[i])) {
      bU[i] = - bl[i];
      bL[i] = - bu[i];
      A.row(i) = - this->row(i);
    } else if(bu[i] == bl[i])
      fidx.emplace_back(make_pair(- T(int(bu[i] == T(0) ? 0 : 1)), i));
    assert(bL[i] <= bU[i] && abs(bL[i]) <= abs(bU[i]));
    A.row(i) /= (T(2) * bU[i] - bL[i]) / T(2);
    assert(isfinite(A.row(i).dot(A.row(i))));
  }
  // N.B. in zeroFix, we get linear Invariant s.t. |Ax| <= 1 possible enough.
        auto res(A.QR().zeroFix(A, fidx));
  const auto z(*this * res * T(int(8)));
        T    t(int(1));
  for(int i = 0; i < z.size(); i ++)
    if(bu[i] * z[i] < T(int(0))) // N.B.: infeasible.
      continue;
    else if(z[i] != T(int(0)))
      t = bl[i] * z[i] < T(int(0)) ? min(t, bu[i] / z[i]) :
                  min(t, min(bu[i] / z[i], bl[i] / z[i]));
  return res *= t;
}

template <typename T> ostream& operator << (ostream& os, const SimpleMatrix<T>& v) {
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

template <typename T> istream& operator >> (istream& is, SimpleMatrix<T>& v) {
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
        v(i, j) = T(int(0));
  }
  return is;
}

template <typename T> static inline T norm2M(const SimpleMatrix<T>& m) {
  auto norm2(m.row(0).dot(m.row(0)));
  for(int i = 1; i < m.rows(); i ++)
    norm2 = max(norm2, m.row(i).dot(m.row(i)));
  return norm2;
}

template <typename T> static inline SimpleMatrix<T> log(const SimpleMatrix<T>& m) {
  static const int cut(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) * T(int(2)) );
  SimpleMatrix<T> res(m.rows(), m.cols());
  const auto c(sqrt(norm2M(m)) * T(2));
  const auto residue(SimpleMatrix<T>(m.rows(), m.cols()).I() - m / c);
        auto buf(residue);
  res.I(log(c));
  for(int i = 1; 0 < i && i < cut; i ++) {
    res -= buf / T(i);
    buf *= residue;
  }
  return res;
}

template <typename T> static inline SimpleMatrix<T> exp01(const SimpleMatrix<T>& m) {
  SimpleMatrix<T> res(m.rows(), m.cols());
  static const int cut(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) * T(int(2)) );
  auto buf(m);
  res.I();
  for(int i = 1; 0 < i && i < cut; i ++) {
    res += buf;
    buf *= m / T(i + 1);
  }
  return res;
}

template <typename T> static inline SimpleMatrix<T> exp(const SimpleMatrix<T>& m) {
  const auto p0(ceil(sqrt(norm2M(m))));
#if defined(_FLOAT_BITS_)
  // XXX:
  myuint p(p0.operator myint());
#else
  myuint p(p0);
#endif
  if(T(myint(int(1))) < abs(p0 - T(myint(p)))) throw "too large abs num in exp matrix";
  auto mm(exp01(m / T(myint(p))));
  auto res(m);
  for(res.I(); p; mm *= mm, p >>= 1)
    if(bool(p & myuint(myint(int(1)))))
      res *= mm;
  return res;
}

template <typename T> static inline SimpleMatrix<T> pow(const SimpleMatrix<T>& m, const T& p) {
  return exp(log(m) * p);
}

template <typename T> SimpleMatrix<complex<T> > dft(const int& size0) {
  const auto size(abs(size0));
  if(! size) {
    const static SimpleMatrix<complex<T> > m0;
    return m0;
  }
  SimpleMatrix<complex<T> > edft( size, size);
  SimpleMatrix<complex<T> > eidft(size, size);
  const auto file(string("./.cache/lieonn/dft-") + to_string(size) +
#if defined(_FLOAT_BITS_)
    string("-") + to_string(_FLOAT_BITS_)
#else
    string("-ld")
#endif
  );
  ifstream cache(file.c_str());
  if(cache.is_open()) {
    cache >> edft;
    cache >> eidft;
    cache.close();
  } else {
    static const auto Pi(T(int(4)) * atan2(T(int(1)), T(int(1))));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < edft.rows(); i ++) {
      for(int j = 0; j < edft.cols(); j ++) {
        const auto theta(- T(int(2)) * Pi * T(i) * T(j) / T(edft.rows()));
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
  const auto file(string("./.cache/lieonn/diff-") + to_string(size) +
#if defined(_FLOAT_BITS_)
    string("-") + to_string(_FLOAT_BITS_)
#else
    string("-ld")
#endif
  );
  ifstream cache(file.c_str());
  if(cache.is_open()) {
    cache >> dd;
    cache >> ii;
    cache.close();
  } else {
    // N.B. if we return recursive each size diff,
    //      taylor series should be broken.
    auto DD(dft<T>(size));
    auto II(dft<T>(size));
    static const auto Pi(T(int(4)) * atan2(T(int(1)), T(int(1))));
    // N.B. we should start this loop with i == 1 on integrate(diff) or inverse.
    //      we also should start with i == 0 on taylor series.
    //      we select latter one.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < DD.rows(); i ++)
      DD.row(i) *= complex<T>(T(int(0)), - T(int(2)) * Pi * T(i) / T(DD.rows()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 1; i < II.rows(); i ++)
      II.row(i) /= complex<T>(T(int(0)), - T(int(2)) * Pi * T(i) / T(DD.rows()));
    // N.B. if we apply DD onto 1 / (1 / f(x)) graph, it's reverse order.
    //      if we average them, it's the only 0 vector.
    // N.B. there exists also completely correct differential matrix,
    //      but it's also be only 0 vector.
    //      (because it's a imaginary part on originals.)
    // N.B. in continuous function, we don't divide dd by &pi;.
    //      (d/dx sum exp(2 Pi i x theta / N) exp(- 2 Pi i y theta / N) f(y)).
    // N.B. from some numerical test, sign on DD, II is reverse side.
    // N.B. In discrete function, we had be chosen ||dd|| := 1, and not now.
    //      (because the (left or right) differential itself cannot be
    //       larger than |x_{k+1}-x_k| < ||x|| in discrete,
    //       sum_0^1 - 2 Pi i (theta/n)^2/2 -&gt; Pi)
    // N.B. if we make plain differential with no error on cosine curve,
    //      it causes constant 0 vector.
    // N.B. if we don't take this real operator, we cannot get better accuracy
    //      result on taylor series.
    dd =   (dft<T>(- size) * DD).template real<T>();
    ii = - (dft<T>(- size) * II).template real<T>();
    ofstream ocache(file.c_str());
    ocache << dd;
    ocache << ii;
    ocache.close();
    cerr << "." << flush;
  }
  return size0 < 0 ? ii : dd;
}

template <typename T> SimpleMatrix<T> diffRecur(const int& size0) {
  const auto size(abs(size0));
  if(! size) {
    static const SimpleMatrix<T> m0;
    return m0;
  }
  SimpleMatrix<T> dd;
  SimpleMatrix<T> ii;
  const auto file(string("./.cache/lieonn/diffrecur-") + to_string(size) +
#if defined(_FLOAT_BITS_)
    string("-") + to_string(_FLOAT_BITS_)
#else
    string("-ld")
#endif
  );
  ifstream cache(file.c_str());
  if(cache.is_open()) {
    cache >> dd;
    cache >> ii;
    cache.close();
  } else {
    dd = SimpleMatrix<T>(size, size).O().setMatrix(0, 0,
         diffRecur<T>(size - 1));
    ii = SimpleMatrix<T>(size, size).O().setMatrix(0, 0,
         diffRecur<T>(- size + 1));
    cerr << "." << flush;
    if(3 < size) {
      dd += diff<T>(  size) +
        SimpleMatrix<T>(size, size).O().setMatrix(1, 1,
          diffRecur<T>(size - 1));
      ii += diff<T>(- size) +
        SimpleMatrix<T>(size, size).O().setMatrix(1, 1,
          diffRecur<T>(- size + 1));
      dd.row(0) /= T(int(2));
      dd.row(dd.rows() - 1) /= T(int(2));
      ii.row(0) /= T(int(2));
      ii.row(dd.rows() - 1) /= T(int(2));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 1; i < dd.rows() - 1; i ++)
        dd.row(i) /= T(int(3));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 1; i < ii.rows() - 1; i ++)
        ii.row(i) /= T(int(3));
      ofstream ocache(file.c_str());
      ocache << dd;
      ocache << ii;
      ocache.close();
    }
  }
  return size0 < 0 ? ii : dd;
}

template <typename T> static inline SimpleVector<T> taylor(const int& size, const T& step) {
  const int  step00(max(int(0), min(size - 1, int(floor(step)))));
  const auto residue0(step - T(step00));
  const auto step0(step00 == size - 1 || abs(residue0) <= T(int(1)) / T(int(2)) ? step00 : step00 + 1);
  const auto residue(step - T(step0));
  SimpleVector<T> res(size);
  res.ek(step0);
  if(residue == T(int(0))) return res;
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

template <typename T> static inline SimpleVector<T> linearInvariant(const SimpleMatrix<T>& in) {
  vector<pair<T, int> > sute;
  return in.QR().zeroFix(in, sute);
}

// N.B. please refer bitsofcotton/randtools.
template <typename T> static inline T makeProgramInvariantPartial(const T& in, const T& ratio, const bool& on01 = false) {
  auto res(on01 ? in :
    ((atan(- in) / atan(T(int(1))) / T(int(2))) + T(int(1))) / T(int(2)) );
  if(res == T(int(0)) ) res = T(int(1));
  assert(T(int(0)) < res && res <= T(int(1)));
  return res *= ratio;
}

template <typename T> static inline pair<SimpleVector<T>, T> makeProgramInvariant(const SimpleVector<T>& in, const T& index = - T(int(1)), const bool& on01 = false) {
  SimpleVector<T> res(in.size() + (T(int(0)) <= index ? 2 : 1));
  res.setVector(0, in);
  res[in.size()] = T(int(1));
  if(T(int(0)) <= index)
    res[in.size() + 1] = T(index);
  T ratio(0);
  for(int i = 0; i < res.size(); i ++)
    ratio += log(res[i] = makeProgramInvariantPartial<T>(res[i],
                            T(int(1)), on01));
  // N.B. x_1 ... x_n == 1.
  // <=> x_1 / (x_1 ... x_n)^(1/n) ... == 1.
  ratio = isfinite(ratio) ? exp(- ratio / T(res.size())) : T(int(1));
  return make_pair(res *= ratio, ratio);
}

template <typename T> static inline T revertProgramInvariant(const pair<T, T>& in, const bool& on01 = false) {
  const auto r0(in.first / in.second);
  const auto r1(T(int(0)) < r0 ? r0 - floor(r0) : ceil(- r0) + r0);
  const auto r(T(int(0)) == r1 ? T(int(1)) : r1);
  return on01 ? r :
      - tan(max(- T(int(1)) + sqrt(SimpleMatrix<T>().epsilon()),
            min(  T(int(1)) - sqrt(SimpleMatrix<T>().epsilon()),
            r * T(int(2)) - T(int(1)) ))
              * atan(T(int(1))) * T(int(2)) );
}

template <typename T> class idFeeder {
public:
  inline idFeeder(const int& size = 1) {
    res.resize(size);
    res.O();
    full = false;
    t = 0;
  }
  inline ~idFeeder() { ; }
  inline const SimpleVector<T>& next(const T& in) {
    if(t < res.size())
      res[t] = in;
    else {
      for(int i = 1; i < res.size(); i ++)
        res[i - 1] = move(res[i]);
      res[res.size() - 1] = in;
    }
    if(res.size() <= ++ t) full = true;
    return res;
  }
  SimpleVector<T> res;
  bool full;
private:
  int  t;
};

template <typename T> const T& sgn(const T& x) {
  static const T zero(0);
  static const T one(1);
  static const T mone(- T(int(1)));
  return x != zero ? (zero < x ? one : mone) : zero;
}

template <typename T> class SimpleSparseVector {
public:
  inline SimpleSparseVector() { ; }
  inline SimpleSparseVector(const int& sute) { assert(sute == 0); }
  inline SimpleSparseVector(const SimpleSparseVector<T>& other) { *this = other; }
  inline SimpleSparseVector(SimpleSparseVector<T>&& other) { *this = other; }
  inline ~SimpleSparseVector() { ; }
  
  inline SimpleSparseVector<T>  operator -  () const {
    auto res(*this);
    for(auto itr(res.entity.begin()); itr != res.entity.end(); ++ itr)
      itr->second = - itr->second;
    return res;
  }
  inline SimpleSparseVector<T>  operator +  (const SimpleSparseVector<T>& other) const {
    auto res(*this);
    return res += other;
  }
  inline SimpleSparseVector<T>& operator += (const SimpleSparseVector<T>& other) {
    for(auto itr(other.entity.begin()); itr != other.entity.end(); ++ itr) {
      if(itr->second == T(int(0))) continue;
      auto search(entity.lower_bound(itr->first));
      if(search == entity.end() || search->first != itr->first)
        (*this)[itr->first] = itr->second;
      else
        search->second += itr->second;
    }
    return *this;
  }
  inline SimpleSparseVector<T>  operator -  (const SimpleSparseVector<T>& other) const {
    auto res(*this);
    return res -= other;
  }
  inline SimpleSparseVector<T>& operator -= (const SimpleSparseVector<T>& other) {
    return *this += - other;
  }
  template <typename U> inline SimpleSparseVector<T>  operator *  (const U& other) const {
    auto res(*this);
    return res *= other;
  }
  template <typename U> inline SimpleSparseVector<T>& operator *= (const U& other) {
    for(auto itr(entity.begin()); itr != entity.end(); ++ itr)
      itr->second *= other;
    return *this;
  }
  template <typename U> inline SimpleSparseVector<T>  operator /  (const U& other) const {
    auto res(*this);
    return res /= other;
  }
  template <typename U> inline SimpleSparseVector<T>& operator /= (const U& other) {
    for(auto itr(entity.begin()); itr != entity.end(); ++ itr)
      itr->second /= other;
    return *this;
  }
  inline       SimpleSparseVector<T>& operator =  (const SimpleSparseVector<T>& other) {
    entity = other.entity;
    return *this;
  }
  inline       SimpleSparseVector<T>& operator =  (SimpleSparseVector<T>&& other) {
    entity = move(other.entity);
    return *this;
  }
  inline       bool                   operator != (const SimpleSparseVector<T>& other) const {
    for(auto itr(entity.begin()); itr != entity.end(); ++ itr)
      if(itr->second != other[itr->first])
        return true;
    for(auto itr(other.entity.begin()); itr != other.entity.end(); ++ itr)
      if(itr->second != const_cast<const SimpleSparseVector<T>&>(*this)[itr->first])
        return true;
    return false;
  }
  inline       bool                   operator == (const SimpleSparseVector<T>& other) const {
    return ! (*this != other);
  }
  inline       T  dot         (const SimpleSparseVector<T>& other) const {
    T res(0);
    for(auto itr(other.entity.begin()); itr < other.entity.end(); ++ itr) {
      auto search(entity.lower_bound(itr->first));
      if(search != entity.end() && search->first == itr->first)
        res += search->second * itr->second;
    }
    return res;
  }
  inline       T& operator [] (const int& idx) {
    assert(0 <= idx);
    const auto search(entity.lower_bound(idx));
    if(search != entity.end() && search->first == idx)
      return search->second;
    else
      entity[idx] = T(int(0));
    const auto search2(entity.lower_bound(idx));
    assert(search2 != entity.end() && search2->first == idx);
    return search2->second;
  }
  inline const T& operator [] (const int& idx) const {
    assert(0 <= idx);
    if(entity.size()) {
      const auto search(entity.lower_bound(idx));
      if(search != entity.end() && search->first == idx)
        return search->second;
    }
    const static T zero(0);
    return zero;
  }
  inline void     clear() {
    entity = map<int, T>();
    return;
  }
  inline       map<int, T>& iter() {
    return entity;
  }
  inline const map<int, T>& iter() const {
    return entity;
  }
private:
  map<int, T>  entity;
};

template <typename T> using SimpleSparseMatrix = SimpleSparseVector<SimpleSparseVector<T> >;
template <typename T> using SimpleSparseTensor = SimpleSparseVector<SimpleSparseVector<SimpleSparseVector<T> > >;


template <typename T> class CatG {
public:
  inline CatG() { ; }
  inline CatG(const int& size0, const vector<SimpleVector<T> >& in);
  inline ~CatG() { ; }
  inline T score(const SimpleVector<T>& in) {
    const auto size(cut.size() - 2);
    assert(0 < size);
    return makeProgramInvariant<T>(tayl(size, in.size()) * in).first.dot(cut) - origin;
  }
  const SimpleMatrix<T>& tayl(const int& size, const int& in) {
    static vector<SimpleMatrix<T> > t;
    if(in < t.size()) {
      if(t[in].rows() && t[in].cols())
        return t[in];
    } else
      t.resize(in + 1, SimpleMatrix<T>());
    t[in].resize(size, in);
    for(int i = 0; i < size; i ++)
      t[in].row(i) = taylor<T>(in, T(i) * T(in) / T(size));
    return t[in];
  }
  SimpleVector<T> cut;
  T   distance;
  T   origin;
};

template <typename T> inline CatG<T>::CatG(const int& size0, const vector<SimpleVector<T> >& in) {
  const auto size(abs(size0));
  SimpleMatrix<T> A(in.size(), size + 1);
  for(int i = 0; i < in.size(); i ++)
    tayl(size, in[i].size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < in.size(); i ++)
    A.row(i) = makeProgramInvariant(tayl(size, in[i].size()) * in[i]).first;
    // N.B. test for linear ones:
    // A.row(i).setVector(0, tayl(size, in[i].size()) * in[i]);
        auto Pt(A.QR());
        auto Ptb(Pt);
  const auto R(Pt * A);
  SimpleVector<T>   one(Pt.cols());
  SimpleVector<int> fix(one.size());
  one.I(T(int(1)));
  fix.I(false);
  for(int n_fixed = 0, idx = 0;
          n_fixed < Pt.rows() - 1;
          n_fixed ++, idx ++) {
    const auto on(Pt.projectionPt(one));
    if(on.dot(on) < one.dot(one) * Pt.epsilon()) {
      Pt = move(Ptb);
      break;
    }
    vector<pair<T, int> > fidx;
    vector<int> pidx;
    fidx.reserve(on.size());
    pidx.resize(on.size(), 0);
    if(0 < size0) {
      for(int i = 0; i < on.size(); i ++)
        fidx.emplace_back(make_pair(abs(on[i]), i));
    } else {
      for(int i = 0; i < on.size(); i ++) {
        T score(0);
        for(int j = 0; j < in.size(); j ++) {
          const auto lscore(abs(on[i] + on[j]));
          if(score == T(int(0)) || lscore < score) {
            score = lscore;
            pidx[i] = j;
          }
        }
        fidx.emplace_back(make_pair(score, i));
      }
    }
    sort(fidx.begin(), fidx.end());
    if(fidx.size() <= idx) break;
    const auto& iidx(fidx[idx].second);
    if(fix[iidx]) continue;
    const auto  orth(Pt.col(iidx));
    const auto  n2(orth.dot(orth));
    if(n2 <= Pt.epsilon()) continue;
    Ptb = Pt;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / n2);
    fix[iidx] = true;
    if(size0 < 0) fix[pidx[iidx]] = true;
  }
  auto ptone(Pt * one);
  if(ptone.dot(ptone) <= Pt.epsilon())
    ptone = Ptb * one;
  cut = R.solve(ptone);
  auto cutn(sqrt(cut.dot(cut)));
  if(cutn != T(int(0))) cut /= cutn;
  auto testv((A * cut).entity);
  sort(testv.begin(), testv.end());
  distance = T(int(0));
  for(int i = 1; i < testv.size(); i ++)
    if(distance <= testv[i] - testv[i - 1]) {
      distance =  testv[i] - testv[i - 1];
      origin   = (testv[i] + testv[i - 1]) / T(int(2));
    }
  return;
}

template <typename T> vector<pair<vector<SimpleVector<T> >, vector<int> > > crush(const vector<SimpleVector<T> >& v, const int& cs, const int& count = 0) {
  assert(0 <= count);
  vector<pair<vector<SimpleVector<T> >, vector<int> > > result;
  if(! v.size() || !v[0].size()) return result;
  int t(0);
  result.emplace_back(pair<vector<SimpleVector<T> >, vector<int> >());
  result[0].first.reserve(v.size());
  result[0].second.reserve(v.size());
  for(int i = 0; i < v.size(); i ++) {
    result[0].first.emplace_back(v[i]);
    result[0].second.emplace_back(i);
  }
  vector<pair<T, pair<int, bool> > > sidx;
  sidx.emplace_back(make_pair(T(int(0)), make_pair(0, false)));
  while(sidx.size() < (count ? count : int(sqrt(T(int(v.size())))))) {
    sort(sidx.begin(), sidx.end());
    int iidx(sidx.size() - 1);
    for( ; - 1 <= iidx; iidx --)
      if(iidx < 0 || (! sidx[iidx].second.second &&
        abs(cs) + 1 < result[sidx[iidx].second.first].first.size() ) ) break;
    if(iidx < 0) break;
    const auto& t(sidx[iidx].second.first);
    CatG<T> catg(cs, result[t].first);
    assert(catg.cut.size());
    vector<SimpleVector<T> > left;
    vector<SimpleVector<T> > right;
    vector<int> lidx;
    vector<int> ridx;
    left.reserve(result[t].first.size());
    right.reserve(result[t].first.size());
    lidx.reserve(result[t].first.size());
    ridx.reserve(result[t].first.size());
    for(int i = 0; i < result[t].first.size(); i ++) {
      const auto score(catg.score(result[t].first[i]));
      (score < T(int(0)) ? left : right).emplace_back(move(result[t].first[i]));
      (score < T(int(0)) ? lidx : ridx).emplace_back(move(result[t].second[i]));
    }
    if(left.size() && right.size()) {
      left.reserve(left.size());
      right.reserve(right.size());
      if(left.size() < right.size()) {
        swap(left, right);
        swap(lidx, ridx);
      }
      result[t].first  = move(left);
      result[t].second = move(lidx);
      sidx[iidx].first  = catg.distance;
      sidx[iidx].second = make_pair(t, false);
      result.emplace_back(make_pair(move(right), move(ridx)));
      sidx.emplace_back(make_pair(catg.distance, make_pair(sidx.size(), false)));
    } else {
      result[t].first  = move(left);
      result[t].second = move(lidx);
      result[t].first.reserve(result[t].first.size() + right.size());
      result[t].second.reserve(result[t].second.size() + ridx.size());
      for(int i = 0; i < right.size(); i ++) {
        result[t].first.emplace_back(move(right[i]));
        result[t].second.emplace_back(move(ridx[i]));
      }
      sidx[iidx].first = catg.distance;
      sidx[iidx].second.second = true;
    }
  }
  return result;
}

template <typename T> static inline vector<pair<vector<SimpleVector<T> >, vector<int> > > crush(const vector<SimpleVector<T> >& v) {
  return crush<T>(v, v[0].size());
}

template <typename T> vector<pair<vector<SimpleVector<T> >, vector<int> > > crushWithOrder(const vector<T>& v, const int& cs, const int& count) {
  vector<SimpleVector<T> > work;
  vector<int> edge;
  // N.B. it's O(v.size()^3 * cs^2).
  work.reserve((v.size() * v.size() * v.size() - 19 * v.size() + 30) / 6);
  edge.reserve(v.size() - 1);
  edge.emplace_back(0);
  for(int i = 3; i < v.size(); i ++) {
    SimpleVector<T> buf(i);
    for(int j = 0; j <= v.size() - i; j ++) {
      for(int k = j; k < j + i; k ++)
        buf[k - j] = v[k];
      if(count < 0) {
        vector<T> wbuf;
        wbuf.reserve(buf.size());
        for(int k = 0; k < buf.size(); k ++)
          wbuf.emplace_back(buf[k]);
        sort(wbuf.begin(), wbuf.end());
        for(int k = 0; k < buf.size(); k ++)
          buf[k] = wbuf[k];
      }
      work.emplace_back(buf);
    }
    edge.emplace_back(work.size());
  }
  auto whole_crush(crush<T>(work, cs, count == - 1 ? 0 : abs(count)));
  vector<pair<vector<SimpleVector<T> >, vector<int> > > res;
  res.reserve(whole_crush.size());
  for(int i = 0; i < whole_crush.size(); i ++) {
    vector<int> idx;
    const auto& sec(whole_crush[i].second);
    idx.reserve(sec.size());
    for(int j = 0; j < sec.size(); j ++)
      idx.emplace_back(sec[j] - *std::lower_bound(edge.begin(), edge.end(), sec[j]));
    res.emplace_back(make_pair(move(whole_crush[i].first), move(idx)));
  }
  return res;
}

template <typename T> static inline vector<pair<vector<SimpleVector<T> >, vector<int> > > crushWithOrder(const vector<T>& v, const int& cs) {
  return crushWithOrder<T>(v, cs, max(int(2), int(sqrt(T(v.size())))));
}

template <typename T> class P012L {
public:
  inline P012L(const int& step = 1, const int& var = 4) {
    assert(1 < var && 0 < step);
    varlen = var;
    this->step = step;
  }
  inline ~P012L() { ; }
  T next(const SimpleVector<T>& in);
private:
  int varlen;
  int step;
};

template <typename T> inline T P012L<T>::next(const SimpleVector<T>& d) {
  static const T zero(int(0));
         auto    M(zero);
  for(int i = 0; i < d.size(); i ++) {
    if(! isfinite(d[i])) return zero;
    M = max(M, abs(d[i]));
  }
  if(M <= zero) return zero;
  vector<SimpleVector<T> > cache;
  cache.reserve(d.size() - varlen + 2);
  for(int i = 0; i < d.size() - varlen - step + 2; i ++) {
    cache.emplace_back(d.subVector(i, varlen));
    cache[cache.size() - 1][varlen - 1] = d[i + varlen + step - 2];
  }
  const auto cat(crush<T>(cache, cache[0].size(), cache.size()));
  SimpleVector<T> work(varlen);
  for(int i = 1; i < work.size(); i ++)
    work[i - 1] = d[i - work.size() + d.size()];
  work[work.size() - 1] = zero;
  auto res(zero);
  auto sscore(zero);
  for(int i = 0; i < cat.size(); i ++) {
    if(! cat[i].first.size()) continue;
    if(! (cat[i].first.size() <= cat[i].first[0].size() + 1)) cerr << "!" << flush;
    SimpleVector<T> avg(cat[i].first[0].size() + 1);
    for(int j = 0; j < cat[i].first.size(); j ++)
      avg += makeProgramInvariant<T>(cat[i].first[j]).first;
    work[work.size() - 1] = T(int(0));
    const auto avg0(avg);
          auto last(sqrt(work.dot(work)));
    const auto navg(avg.dot(avg));
    if(! isfinite(navg) || navg == zero) continue;
    for(int ii = 0;
            ii < 2 * int(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) )
            && sqrt(work.dot(work) * SimpleMatrix<T>().epsilon()) <
                 abs(work[work.size() - 1] - last); ii ++) {
      last  = work[work.size() - 1];
      const auto vdp(makeProgramInvariant<T>(work));
      avg  = avg0 * sqrt(vdp.first.dot(vdp.first) / avg0.dot(avg0));
      work[work.size() - 1] =
        revertProgramInvariant<T>(make_pair(avg[varlen - 1] /
             T(int(avg.size())), vdp.second));
    }
    const auto vdp(makeProgramInvariant<T>(work));
    T score(0);
    for(int j = 0; j < work.size(); j ++)
      score += work[j] * revertProgramInvariant<T>(make_pair(avg[j], vdp.second));
    res += score * work[work.size() - 1];
    sscore += abs(score);
  }
  return sscore == zero ? sscore : res / sscore;
}


template <typename T> SimpleVector<T> pnext(const int& size, const int& step = 1, const int& r = 1) {
#if defined(_ADVANCE_PNEXT_BITS_)
  typedef DUInt<uint64_t, 64>          pnext_uint128_t;
  typedef DUInt<pnext_uint128_t, 128>  pnext_uint256_t;
  typedef DUInt<pnext_uint256_t, 256>  pnext_uint512_t;
  typedef Signed<pnext_uint512_t, 512> pnext_int512_t;
  typedef pnext_uint512_t pnext_uint;
  typedef pnext_int512_t  pnext_int;
  typedef SimpleFloat<pnext_uint, DUInt<pnext_uint, 512>, 512, pnext_int> pnext_float;
  SimpleVector<T> res;
  const auto file(string("./.cache/lieonn/pnextadv-") + to_string(size) +
    string("-") + to_string(step) + string("-") + to_string(r));
  ifstream cache(file.c_str());
  if(cache.is_open()) {
    cache >> res;
    cache.close();
  } else {
    cerr << "." << flush;
    auto work(taylor<pnext_float>(size * r, pnext_float(step * r < 0 ? step * r : (size + step) * r - 1)));
    for(int i = 1; i < r; i ++)
      work += taylor<pnext_float>(size * r, pnext_float(step * r < 0 ? step * r + i : (size + step) * r - 1 - i));
    const auto pn((dft<pnext_float>(- size * r).subMatrix(0, 0, size * r, size) * dft<pnext_float>(size)).template real<pnext_float>().transpose() * work);
    res.resize(pn.size());
    for(int i = 0; i < res.size(); i ++) {
      auto abspni(abs(pn[i]));
      res[i] = T(abspni.operator int());
      for(int j = 0; j < 512 / 16; j ++) {
        abspni -= floor(abspni);
        abspni *= pnext_float(int(65536));
        res[i] += T(abspni.operator int()) / pow(T(int(65536)), T(int(j + 1)) );
      }
      if(pn[i] < pnext_float(int(0))) res[i] = - res[i];
    }
    ofstream ocache(file.c_str());
    ocache << res;
    ocache.close();
  }
  return res;
#else
  auto work(taylor(size * r, T(step * r < 0 ? step * r : (size + step) * r - 1)));
  for(int i = 1; i < r; i ++)
    work += taylor(size * r, T(step * r < 0 ? step * r + i : (size + step) * r - 1 - i));
  return (dft<T>(- size * r).subMatrix(0, 0, size * r, size) * dft<T>(size)).template real<T>().transpose() * work;
#endif
}

template <typename T> SimpleVector<T> minsq(const int& size) {
  assert(1 < size);
  const T xsum(size * (size - 1) / 2);
  const T xdot(size * (size - 1) * (2 * size - 1) / 6);
  const auto denom(xdot * T(size) - xsum * xsum);
  SimpleVector<T> s(size);
  for(int i = 0; i < s.size(); i ++)
    s[i] = (T(i) * T(size) - xsum) / denom;
  return s;
}

template <typename T> const SimpleVector<T>& pnextcacher(const int& size, const int& step, const int& r) {
  assert(0 < size && 0 <= step && 0 < r);
  static vector<vector<vector<SimpleVector<T> > > > cp;
  if(cp.size() <= size)
    cp.resize(size + 1, vector<vector<SimpleVector<T> > >());
  if(cp[size].size() <= step)
    cp[size].resize(step + 1, vector<SimpleVector<T> >());
  if(cp[size][step].size() <= r)
    cp[size][step].resize(r + 1, SimpleVector<T>());
  if(cp[size][step][r].size()) return cp[size][step][r];
  return cp[size][step][r] = pnext<T>(size, step, r);
}

template <typename T> const SimpleVector<T>& mscache(const int& size) {
  assert(0 < size);
  static vector<SimpleVector<T> > ms;
  if(ms.size() <= size) ms.resize(size + 1, SimpleVector<T>());
  if(ms[size].size()) return ms[size];
  return ms[size] = minsq<T>(size);
}

template <typename T> const SimpleMatrix<complex<T> >& dftcache(const int& size) {
  assert(size != 0);
  static vector<SimpleMatrix<complex<T> > > cdft;
  static vector<SimpleMatrix<complex<T> > > cidft;
  if(0 < size) {
    if(cdft.size() <= size) cdft.resize(size + 1, SimpleMatrix<complex<T> >());
    if(cdft[size].rows() && cdft[size].cols()) return cdft[size];
    return cdft[size] = dft<T>(size);
  }
  if(cidft.size() <= abs(size)) cidft.resize(abs(size) + 1, SimpleMatrix<complex<T> >());
  if(cidft[abs(size)].rows() && cidft[abs(size)].cols()) return cidft[abs(size)];
  return cidft[abs(size)] = dft<T>(size);
}

template <typename T, int r = 4> class P0 {
public:
  inline P0(const int& step = 1) {
    this->step = step;
  }
  inline ~P0() { ; };
  inline T next(const SimpleVector<T>& in) {
    return pnextcacher<T>(in.size(), ((step - 1) % in.size()) + 1, r).dot(in);
  }
  int step;
};

template <typename T, typename P> class P0inv {
public:
  inline P0inv() { ; }
  inline P0inv(P&& p) { this->p = p; }
  inline ~P0inv() { ; }
  inline T next(const SimpleVector<T>& in) {
    static const T zero(int(0));
    static const T one(int(1));
    auto ff(in);
    for(int i = 0; i < in.size(); i ++) if(in[i] == zero) return in[in.size() - 1];
    else ff[i] = one / in[i];
    const auto pn(p.next(ff));
    if(pn == zero) return in[in.size() - 1];
    return one / pn;
  }
  P p;
};

template <typename T, typename P, typename feeder> class P0DFT {
public:
  inline P0DFT() { ; }
  inline P0DFT(P&& p, const int& size) {
    f = feeder(size);
    (this->p).resize(size, p);
    q = this->p;
  }
  inline ~P0DFT() { ; };
  inline T next(const T& in) {
    const auto& fn(f.next(in));
    if(! f.full) return T(int(0));
    auto ff(dftcache<T>(fn.size()) * fn.template cast<complex<T> >());
    assert(ff.size() == p.size() && p.size() == q.size());
    for(int i = 0; i < ff.size(); i ++)
      ff[i] = complex<T>(p[i].next(ff[i].real()), q[i].next(ff[i].imag()));
/*
      if(! (ff[i].real() == T(int(0)) && ff[i].imag() == T(int(0)) ) )
        ff[i] = abs(p[i].next(abs(ff[i]))) * exp(complex<T>(T(int(0)), q[i].next(arg(ff[i]))));
*/
    return dftcache<T>(- fn.size()).row(fn.size() - 1).dot(ff).real();
  }
  vector<P> p;
  vector<P> q;
  feeder f;
};

template <typename T, typename P> class northPole {
public:
  inline northPole() { ; }
  inline northPole(P&& p) { this->p = p; }
  inline ~northPole() { ; }
  inline T next(const SimpleVector<T>& in) {
    static const T zero(int(0));
    static const T one(int(1));
    static const T M(atan(one / sqrt(SimpleMatrix<T>().epsilon())));
    auto ff(in);
    for(int i = 0; i < in.size(); i ++)
      if(! isfinite(in[i]) || in[i] == zero) return in[in.size() - 1];
      else {
        ff[i] = atan(in[i]);
        // assert(- M < ff[i] && ff[i] < M);
        // N.B. we don't avoid right hand side, it's harmless.
        // ff[i] = atan(one / ff[i]);
        // assert(- M < ff[i] && ff[i] < M);
      }
    auto work(p.next(ff));
    // if(! isfinite(work) || work == zero) return in[in.size() - 1];
    if(! isfinite(work)) return in[in.size() - 1];
    // work = tan(max(- M, min(M, one / tan(max(- M, min(M, work))))));
    work = tan(max(- M, min(M, work)));
    if(isfinite(work)) return work;
    return in[in.size() - 1];
  }
  P p;
};

template <typename T, typename P, bool avg = false> class sumChain {
public:
  inline sumChain() { ; }
  inline sumChain(P&& p) { this->p = p; }
  inline ~sumChain() { ; }
  inline T next(const SimpleVector<T>& in) {
    auto ff(in);
    for(int i = 1; i < ff.size(); i ++)
      ff[i] += ff[i - 1];
    if(! avg) return p.next(ff) - ff[ff.size() - 1];
    const auto A(ff[ff.size() - 1] / T(ff.size()));
    for(int i = 0; i < ff.size(); i ++)
      ff[i] = in[i] - A;
    return p.next(ff) + A;
  }
  P p;
};

template <typename T, typename P> class logChain {
public:
  inline logChain() { ; }
  inline logChain(P&& p) { this->p = p; }
  inline ~logChain() { ; }
  inline T next(const SimpleVector<T>& in) {
    static const T zero(int(0));
    static const T one(int(1));
    auto ff(in);
    if(ff[0] == zero) return in[in.size() - 1];
    for(int i = 1; i < ff.size(); i ++)
      if((ff[i] += ff[i - 1]) == zero) return in[in.size() - 1];
    SimpleVector<T> gg(ff.size() - 1);
    gg.O();
    for(int i = 1; i < ff.size(); i ++)
      if(! isfinite(gg[i - 1] = ff[i] / ff[i - 1] - one)) return in[in.size() - 1];
    return p.next(gg) * ff[ff.size() - 1];
  }
  P p;
};

template <typename T> class P0maxRank0 {
public:
  inline P0maxRank0(const int& step = 1) {
    p = p0_0t(P0<T>(step));
    q = p0_i0t(p0_0t(P0<T>(step)));
  }
  inline ~P0maxRank0() { ; }
  inline T next(const SimpleVector<T>& in) {
    return (p.next(in) + q.next(in)) / T(int(2));
  }
  // N.B. on existing taylor series.
  //      if the sampling frequency is not enough, middle range of the original
  //      function frequency (enough large bands) will effect prediction fail.
  //      this is because we only observes highest and lowest frequency on
  //      sampling points, so omitted part exists.
  //      even if the parameter on P0 is large, situation unchange.
  //      so we should use sectional measurement for them.
  // N.B. the sectional measurament is done by following Ppad class.
  //      So this is only the raw prediction.
  typedef sumChain<T, P0<T>, true> p0_0t;
  typedef P0inv<T, p0_0t> p0_i0t;
  p0_0t p;
  p0_i0t q;
};

template <typename T> class P0maxRank {
public:
  inline P0maxRank(const int& step = 1) {
    p = p0_t(p0_2t(p0_1t(p0_0t(step))));
  }
  inline ~P0maxRank() { ; }
  inline T next(const SimpleVector<T>& in) {
    return p.next(in);
  }
/*
  // N.B. make information-rich not to associative/commutative.
  //      2 dimension semi-order causes (x, status) from input as sedenion.
  // N.B. we need only once P0DFT in general because associative condition
  //      is necessary for input ordering.
  typedef P0DFT<T, p0_1t, idFeeder<T> > p0_2t;
  // N.B. on any R to R into reasonable taylor.
  typedef northPole<T, p0_2t> p0_6t;
  typedef northPole<T, p0_6t> p0_7t;
  // N.B. we treat periodical part as non aligned complex arg part.
  typedef logChain<T, p0_7t>  p0_8t;
  typedef logChain<T, p0_8t>  p0_9t;
  // N.B. we make the prediction on (delta) summation.
  typedef sumChain<T, p0_9t>  p0_10t;
  // N.B. we take average as origin of input.
  typedef sumChain<T, p0_10t, true> p0_t;
  // N.B. this needs huge memory to run.
*/
  // N.B. plain complex form.
  typedef P0maxRank0<T> p0_0t;
  typedef northPole<T, p0_0t>  p0_1t;
/*
  typedef northPole<T, p0_1t> p0_2t;
  typedef logChain<T, p0_2t>  p0_3t;
  typedef logChain<T, p0_3t>  p0_4t;
  typedef sumChain<T, p0_4t>  p0_5t;
  typedef sumChain<T, p0_5t, true> p0_t;
*/
  // N.B. we only handle lebesgue measurable and R(finite)-valued functions.
  //      so worse structures are handled by P01.
  typedef sumChain<T, p0_1t>  p0_2t;
  typedef sumChain<T, p0_2t, true> p0_t;
  p0_t p;
};

// Get invariant structure that
// \[- &alpha, &alpha;\[ register computer with deterministic calculation.
// cf. bitsofcotton/randtools .
template <typename T, bool nonlinear = true> class P01 {
public:
  inline P01(const int& step = 1, const int& var = 4) {
    assert(0 < var && 0 < step);
    this->varlen = var;
    this->step = step;
  }
  inline ~P01() { ; }
  inline T next(const SimpleVector<T>& in) {
    static const T zero(0);
    static const T one(1);
    static const T two(2);
    // N.B. please use catgp to compete with over learning.
    // XXX: division accuracy glitch.
    const auto nin(sqrt(in.dot(in) * (one + SimpleMatrix<T>().epsilon())));
    if(! isfinite(nin) || nin == zero) return zero;
    SimpleMatrix<T> invariants(3, nonlinear ? varlen + 2 : varlen);
    invariants.O();
    for(int i0 = 0; i0 < invariants.rows(); i0 ++) {
      SimpleMatrix<T> toeplitz(in.size() - varlen - step + 2
                               - invariants.rows() + 1, invariants.cols());
      for(int i = i0; i < toeplitz.rows() + i0; i ++) {
        auto work(in.subVector(i, varlen));
        work[work.size() - 1] = in[i + varlen + step - 2];
        toeplitz.row(i - i0) = nonlinear ? makeProgramInvariant<T>(move(work),
          T(i + 1) / T(toeplitz.rows() + 1) ).first : move(work);
      }
      invariants.row(i0) = linearInvariant<T>(toeplitz);
    }
    SimpleVector<T> invariant(invariants.cols());
    invariant.O();
    for(int i = 0; i < invariants.cols(); i ++)
      invariant[i] = P0maxRank0<T>().next(invariants.col(i));
    if(invariant[varlen - 1] == zero) return zero;
    SimpleVector<T> work(varlen);
    for(int i = 1; i < work.size(); i ++)
      work[i - 1] = in[i - work.size() + in.size()];
    work[work.size() - 1] = zero;
    if(nonlinear) {
      auto last(sqrt(work.dot(work)));
      for(int ii = 0;
              ii < 2 * int(- log(SimpleMatrix<T>().epsilon()) / log(two) )
              && sqrt(work.dot(work) * SimpleMatrix<T>().epsilon()) <
                   abs(work[work.size() - 1] - last); ii ++) {
        last = work[work.size() - 1];
        const auto work2(makeProgramInvariant<T>(work, one));
        work[work.size() - 1] = revertProgramInvariant<T>(make_pair(
                 - (invariant.dot(work2.first) -
                        invariant[varlen - 1] * work2.first[varlen - 1]) /
                   invariant[varlen - 1], work2.second));
      }
      return work[work.size() - 1];
    }
    return - invariant.dot(work) / invariant[varlen - 1];
  }
private:
  int varlen;
  int step;
};

// N.B. we omit high frequency part (1/f(x) input) to be treated better in P.
template <typename T, typename P> class PBond {
public:
  inline PBond() { ; }
  inline PBond(const int& status, P&& p = P()) {
    assert(0 < status);
    this->p = p;
    f = idFeeder<T>(status);
    M = T(int(1));
  }
  inline ~PBond() { ; }
  inline T next(const T& in) {
    M = max(M, abs(in));
    auto g(f.next(in));
    if(! f.full) return T(int(0));
    // N.B. with 1-norm normalized input:
    T m(g[0] /= M);
    for(int i = 1; i < g.size(); i ++) m = min(m, g[i] /= M);
    // N.B. offset const.
    m -= T(int(1));
    // N.B. 0 < v, normalize with v's orthogonality:
    T mavg(log(g[0] - m));
    for(int i = 1; i < g.size(); i ++) mavg += log(g[i] - m);
    mavg /= T(int(g.size()));
    mavg  = exp(mavg);
    // N.B. we need nonlinear prediction, so * M before to predict.
    return max(- M, min(M, p.next(g / mavg * M) * mavg));
  }
  idFeeder<T> f;
  P p;
  T M;
};

// N.B. this class feeds new states into predictor as a new dimension of
//      linear sum.
template <typename T, typename P> class Pprogression {
public:
  inline Pprogression() { ; }
  inline Pprogression(const int& loop0, const int& istat) {
    assert(loop0);
    const auto loop(abs(loop0));
    p.reserve(loop);
    for(int i = 0; i < loop; i ++)
      p.emplace_back(PBond<T, P>(istat + i, P(i + 1)));
    h  = idFeeder<T>(loop);
    {
      vector<int> ph0;
      ph0.resize(loop, 0);
      ph.resize(loop, ph0);
      vector<T> eh0;
      eh0.resize(loop, T(int(0)));
      eh.resize(loop, eh0);
    }
    addp = 0 < loop0;
    bb.reserve(loop);
    for(int i = 0; i < loop; i ++)
      bb.emplace_back(idFeeder<T>(addp ? i + 1 : 1));
    t ^= t;
    this->istat = istat;
  }
  inline ~Pprogression() { ; }
  inline const T& progression(const SimpleVector<T>& h, const int& idx, const int& count) {
    assert(0 <= idx && 0 <= count);
    if(! count) return h[idx];
    if(ph[idx][count]) return eh[idx][count];
    ph[idx][count] = 1;
    return (eh[idx][count] = progression(h, idx, count - 1) - progression(h, idx - 1, count - 1));
  }
  inline T next(const T& in) {
    static const T zero(int(0));
    for(int i = 0; i < ph.size(); i ++)
      for(int j = 0; j < ph[i].size(); j ++)
        ph[i][j] = 0;
    const auto& hh(h.next(in));
    auto M(zero);
    if(! h.full) return M;
    for(int i = 0; i < p.size(); i ++)
      if(p.size() - 1 - i <= t)
        bb[i].next(p[i].next(progression(hh, hh.size() - 1, i)));
    if(p.size() <= t)
      for(int i = 0; i < p.size(); i ++) {
        M += bb[i].res[0];
        if(addp) for(int j = i - 1; 0 <= j; j --)
          M += progression(hh, hh.size() - 1, j);
      }
    t ++;
    return addp ? M /= T(int(p.size())) : M;
  }
  vector<PBond<T, P> > p;
  idFeeder<T> h;
  vector<vector<int> > ph;
  vector<vector<T> > eh;
  vector<idFeeder<T> > bb;
  int t;
  int istat;
  bool addp;
};

template <typename T, typename P> class PpersistentOnce {
public:
  inline PpersistentOnce() { ; }
  inline ~PpersistentOnce() { ; }
  inline const T& progression(const SimpleVector<T>& h, const int& idx, const int& count) {
    assert(0 <= idx && 0 <= count);
    if(! count) return h[idx];
    if(ph[idx][count]) return eh[idx][count];
    ph[idx][count] = 1;
    return (eh[idx][count] = progression(h, idx, count - 1) - progression(h, idx - 1, count - 1));
  }
  inline T next(const SimpleVector<T>& in, const int& istat) {
    {
      vector<T> eh0;
      vector<bool> ph0;
      eh0.resize(in.size(), T(int(0)));
      ph0.resize(in.size(), false);
      eh.resize(0);
      ph.resize(0);
      eh.resize(in.size(), eh0);
      ph.resize(in.size(), ph0);
    }
    // N.B. use full of the input to reduce counter measure.
    //      if our strategy is to feed sloppy predictor enough internal status,
    //      this can be reduced. either, if ours is to reduce enough,
    //      this can be increased. we balance them middle.
    const int nretry(in.size() / 2 - istat / 2);
    T res(int(0));
    for(int i = 0; i <= nretry; i ++) {
      idFeeder<T> buf(in.size() - nretry);
      for(int j = i; j < in.size() - nretry + i; j ++)
        buf.next(progression(in, j, i));
      assert(buf.full);
      res += P(nretry - i + 1).next(buf.res);
      for(int j = i - 1; 0 <= j; j --)
        res += progression(in, in.size() - 1, j);
    }
    return res /= T(nretry + 1);
  }
  vector<vector<T> >    eh;
  vector<vector<bool> > ph;
};

template <typename T, typename P, typename Q> class PAthenB {
public:
  inline PAthenB() { ; }
  inline PAthenB(P&& p, Q&& q) { this->p = p; this->q = q; M = T(int(0)); }
  inline ~PAthenB() { ; }
  inline T next(const T& in) {
    const auto M2(q.next(in * M));
    return M2 * (M = p.next(in));
  }
  P p;
  Q q;
  T M;
};

template <typename T> using P10 =
   PAthenB<T, Pprogression<T, P01<T> >,
              Pprogression<T, P0maxRank<T> > >;
template <typename T> using P210 =
   PAthenB<T, Pprogression<T, P012L<T> >, P10<T> >;

// N.B. invariant gathers some of the group on the input pattern.
template <typename T> SimpleMatrix<T> concat(const SimpleMatrix<T>& m0, const SimpleMatrix<T>& m1) {
  // det diag result = det diag m0 + det diag m1
  // [1 x x^reverse 1]
  // if we met rank shrink, assert exit then.
  // we can handle this with compiling operation adding period-depend values.
  assert(m0.rows() == m1.rows() && m0.cols() == m1.cols());
  SimpleMatrix<T> work0(m0);
  SimpleMatrix<T> work1(m1);
  auto res(m0);
  for(int i = 0; i < m0.rows(); i ++) {
    auto qw1(work1.transpose().QR());
    auto rw1(qw1 * work1.transpose());
    // XXX : assert exit here.
    work0 = (rw1.inverse() * qw1 * work0.transpose()).transpose();
    assert(work0.rows() == work0.cols());
    SimpleMatrix<T> lwork(work0.rows() * 2, work0.cols() * 2);
    const auto ii(SimpleMatrix<T>(work0.rows(), work0.cols()).I());
    lwork.setMatrix(0, 0, ii).setMatrix(0, work0.cols(), work0 - ii).setMatrix(work0.rows(), work0.cols(), ii);
    for(int j = 0; j < work0.rows(); j ++)
      for(int k = 0; k < work0.cols(); k ++)
        lwork(lwork.rows() - j, work0.cols() - k) = work0(j, k);
    for(int j = 0; j < lwork.rows() / 2; j ++)
      for(int k = 0; k < lwork.cols(); k ++)
        swap(lwork(j, k), lwork(lwork.rows() - j, k));
    for(int j = 0; j < lwork.rows(); j ++)
      for(int k = 0; k < lwork.cols() / 2; k ++)
        swap(lwork(j, k), lwork(j, lwork.cols() - k));
    // LDLt:
/*
    auto L();
    auto D();
    auto Linv();
    for(int j = 0; j < L.rows() / 2; j ++)
      for(int k = 0; k < L.cols(); k ++)
        swap(L(j, k), L(L.rows() - j, k));
    for(int j = 0; j < Linv.rows(); j ++)
      for(int k = 0; k < Linv.cols() / 2; k ++)
        swap(Linv(j, k), Linv(j, Linv.cols() - k));
*/
    // factor apply res, work0, work1:
  }
  return res;
}


template <typename T> SimpleMatrix<T> diff(const SimpleMatrix<T>& m, const int& idx) {
  SimpleMatrix<T> res(m.rows() - 1, m.cols());
  res.O();
  for(int i = 0; i < m.rows(); i ++) {
    auto lres(res);
    lres.O();
    for(int j = 0; j < i; j ++) lres.row(j) = m.row(j);
    for(int j = i + 1; j < m.rows(); j ++) lres.row(j - 1) = m.row(j);
    if(m(i, i) == T(int(0))) continue;
    else if(m(i, i) < T(int(0))) lres.row(0) = - lres.row(0);
    res = concat(res, lres /= pow(abs(m(i, i)), T(int(1)) / T(int(lres.rows()))));
  }
}

template <typename T> SimpleMatrix<T> integrate(const SimpleMatrix<T>& m, const int& idx, const int& stage = 0) {
  // N.B. S^x det diag Ax dx (= S u'v) =
  //  (S^x dx) * det diag Ax (= S(uv)') -
  //  S^x(x(det diag Ax)')dx (= S uv')
  //      S^x(x(det diag Ax)')dx   (= S u'v) =
  //  (S^x x dx) * (det diag Ax)'  (= S(uv)') -
  //  S^x(x^2/2 (det diag Ax)'')dx (= S uv')
  //    ...
  SimpleMatrix<T> factorial(m.rows(), m.cols());
  factorial.O();
  for(int i = 0; i < factorial.rows(); i ++)
    factorial(i, idx) = T(int(1)) / T(int(i + 1));
  if(stage == m.rows() - 1) return factorial * m(m.rows() - 1, idx);
  return concat(m, factorial.setMatrix(stage + 1, 0, integrate(diff(m, idx), idx, stage + 1)), true);
}

// N.B. we need huge computing power depends on m at least O(m.rows()^4).
template <typename T> SimpleVector<T> reduce(const SimpleMatrix<T> m) {
  SimpleMatrix<T> work(m);
  for(int i = 0; i < m.rows() - 1; i ++)
    work = integrate(work, i);
  for(int i = 0; i < m.rows() - 1; i ++)
    work = diff(work, i);
  return work.row(0);
}


template <typename T> class Decompose {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline Decompose(const int& size = 0) {
    assert(0 <= size);
    this->size = size;
  }
  inline ~Decompose() { ; }
  const vector<Mat>& A() {
    static vector<vector<Mat> > mA;
    if(mA.size() <= size) mA.resize(size + 1, vector<Mat>());
    auto& a(mA[size]);
    if(a.size() < size) {
      a.reserve(size);
      for(int i = 0; i < size; i ++) {
        SimpleMatrix<T> AA(size, size);
        for(int j = 0; j < AA.rows(); j ++) {
          const auto jj(T(j) * T(i + 2) / T(size + 1));
          AA.row(j) = taylor<T>(AA.cols(), (jj - floor(jj)) * T(size - 1));
        }
        a.emplace_back(move(AA));
      }
      for(int i = 1; i < a.size(); i ++)
        swap(a[a.size() - i], a[a.size() - i - 1]);
    }
    return a;
  }
  inline Vec  mimic(const Vec& dst, const Vec& src, const T& intensity = T(1)) {
    const auto size2(dst.size() / size);
    const auto size3(src.size() / size);
          auto res(dst);
    for(int i = 0; i < size2; i ++) {
      const auto dd(prepare(dst, i));
      apply(res, synth(mother(prepare(src, i * size3 / size2)),
                       freq(mother(dd), dd)) * intensity +
                 dd * (T(1) - intensity), dd, i);
    }
    return res;
  }
  inline Vec  emphasis(const Vec& dst, const T& intensity = T(1)) {
    const auto size2(dst.size() / size);
          auto res(dst);
    for(int i = 0; i < size2; i ++) {
      const auto dd(prepare(dst, i));
            auto lfreq(dd);
      for(int j = 0; j < lfreq.size(); j ++)
        lfreq[j] = T(j) / T(lfreq.size());
      apply(res, synth(mother(dd), lfreq) * intensity +
                 dd * (T(1) - intensity), dd, i);
    }
    return res;
  }
  inline Vec  mother(const Vec& in) {
    vector<Mat> A0;
    if(A0.size() <= size) A0.resize(size + 1, Mat());
    auto& a0(A0[size]);
    if(a0.rows() != size || a0.cols() != size) {
      const auto& a(A());
      a0 = a[0];
      for(int i = 1; i < a.size(); i ++)
        a0 += a[i];
    }
    return a0.solve(in);
  }
  inline Vec  freq(const Vec& mother, const Vec& in) {
    assert(size == mother.size() && size == in.size());
    Mat work(size, size);
    for(int i = 0; i < size; i ++)
      work.setCol(i, A()[i] * mother);
    return work.solve(in);
  }
  inline Vec  synth(const Vec& mother, const Vec& in) {
    assert(size == mother.size() && size == in.size());
    Vec res(size);
    for(int i = 0; i < size; i ++)
      res[i] = T(0);
    for(int i = 0; i < size; i ++)
      res += A()[i] * mother * in[i];
    return res;
  }
  inline Mat subImage(const Mat& img, const int& x, const int& y, const int& r) const {
    Mat res(size, size);
    for(int i = 0; i < res.rows(); i ++)
      for(int j = 0; j < res.cols(); j ++) {
        const auto rr(T(j + 1) / T(res.cols()) * T(r));
        const auto th(T(i) / T(res.rows()) * T(2) * T(4) * atan2(T(1), T(1)));
        res(i, j) = img(flip(x + int(rr * cos(th)), img.rows()),
                        flip(y + int(rr * sin(th)), img.cols()));
      }
    return res;
  }
  Vec  enlarge(const Vec& in, const int& r = 2);
  Mat  represent(const Mat& img, const int& depth = 3);
private:
  inline Vec  prepare(const Vec& in, const int& idx = 0) const {
    const auto cnt(in.size() / size);
    assert(0 < cnt);
    Vec res(size);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < size; i ++)
      res[i] = i * cnt + idx < in.size() ? in[i * cnt + idx] : T(0);
    return res;
  }
  inline void apply(Vec& v, const Vec& dst, const Vec& src, const int& idx = 0) const {
    assert(size && dst.size() == size && src.size() == size);
    const auto cnt(v.size() / size);
    assert(0 < cnt);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < size; i ++)
      if(i * cnt + idx < v.size())
        v[i * cnt + idx] += dst[i] - src[i];
    return;
  }
  inline int flip(const int& x, const int& s) const {
    const int xx(abs(x % (s * 2)));
    return s <= xx ? s * 2 - xx - 1 : xx;
  }
  int  size;
};

template <typename T> typename Decompose<T>::Vec Decompose<T>::enlarge(const Vec& in, const int& r) {
  assert(0 < r && size == in.size());
  static vector<vector<Mat> > p;
  if(p.size() <= size)
    p.resize(size + 1, vector<Mat>());
  if(p[size].size() <= r)
    p[size].resize(r + 1, Mat());
  auto& pp(p[size][r]);
  if(pp.rows() < size * r) {
    pp.resize(size * r, size);
    for(int i = 0; i < pp.rows(); i ++)
      pp.row(i) = taylor<T>(size, T(i) * T(size - 1) / T(pp.rows() - 1));
  }
  const auto m(mother(in));
  const auto f2(freq(m, in));
  const auto bm(pp * m);
        auto ff(bm);
  for(int i = 0; i < ff.size(); i ++)
    ff[i] = i ? f2[(i % (f2.size() - 1)) + 1] : f2[i];
  Decompose<T> ee(size * r);
  auto result(ee.synth(bm, ff));
  return result *= sqrt(in.dot(in) / result.dot(result) * T(r));
}

template <typename T> typename Decompose<T>::Mat Decompose<T>::represent(const Mat& img, const int& depth) {
  Mat res0(1, size);
  Mat w00(img.rows() - size * 2, size);
  const auto int4(diff<T>(- size) * diff<T>(- size) * diff<T>(- size) * diff<T>(- size));
  const auto int4t(int4.transpose());
  for(int i = size; i < img.rows() - size; i ++) {
    Mat w0(img.cols() - size * 2, size);
    for(int j = size; j < img.cols() - size; j ++) {
      vector<Vec> w1;
      for(int r = size;
              r < min(min(i, j),
                    min(img.rows() - i - 1, img.cols() - j - 1));
              r += max(int(1), min(img.rows() / size, img.cols() / size))) {
        // integrate 4th times because we svd 4 times.
        // svd takes bothside transform, we suppose them as differential op.
        const auto part(int4 * subImage(img, i, j, r) * int4t);
        const auto left(part.SVD() * part);
              Vec  work(left.rows());
        for(int k = 0; k < work.size(); k ++)
          work[k] = sqrt(left.row(k).dot(left.row(k))) + T(1);
        // N.B. normalized singular values on the image with circle region.
        //      If this is flat, the data we have is flat.
        //      If this is edged, the data we have has some data.
        work = mother(work);
        // N.B. enlarging specific bias.
        //      recursive on them.
        w1.emplace_back(move(work /= sqrt(work.dot(work))));
      }
      if(! w1.size())
        for(int k = 0; k < w0.cols(); k ++)
          w0(j - size, k) = T(1) / sqrt(T(w0.cols()));
      else if(w1.size() == 1)
        // N.B. line intensity.
        w0.row(j - size) = move(w1[0]);
      else {
        Mat w1m(w1.size(), size);
        for(int i = 0; i < w1m.rows(); i ++)
          w1m.row(i) = move(w1[i]);
        w1m = w1m.transpose();
        const auto left(w1m.SVD() * w1m);
        for(int k = 0; k < left.rows(); k ++)
          w0(j - size, k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
        w0.row(j - size)  = mother(w0.row(j - size));
        w0.row(j - size) /= sqrt(w0.row(j - size).dot(w0.row(j - size)));
      }
    }
    // N.B. do same on x axis: w0 = w0.transpose();
    for(int j = 0; j < w0.rows(); j ++)
      for(int k = 0; k < w0.cols(); k ++)
        assert(isfinite(w0(j, k)) && ! isnan(w0(j, k)));
    const auto left(w0.SVD() * w0);
    for(int k = 0; k < left.rows(); k ++)
      w00(i - size, k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
    w00.row(i - size)  = mother(w00.row(i - size));
    w00.row(i - size) /= sqrt(w00.row(i - size).dot(w00.row(i - size)));
  }
  // N.B. do same on whole image:
  w00 = w00.transpose();
  const auto left(w00.SVD() * w00);
  for(int k = 0; k < left.rows(); k ++)
    res0(0, k) = sqrt(left.row(k).dot(left.row(k))) + T(1);
  res0.row(0)  = mother(res0.row(0));
  res0.row(0) /= sqrt(res0.row(0).dot(res0.row(0)));
  // N.B. recursive on them.
  if(0 < depth && size * 4 <= min(img.rows(), img.cols()) / 2) {
    Mat dimg[5];
    for(int i = 0; i < 5; i ++)
      dimg[i] = Mat(img.rows() / 2, img.cols() / 2);
    for(int i = 0; i < dimg[0].rows(); i ++)
      for(int j = 0; j < dimg[0].cols(); j ++) {
        dimg[0](i, j) = img(i, j);
        dimg[1](i, j) = img(i - dimg[0].rows() + img.rows(), j);
        dimg[2](i, j) = img(i, j - dimg[0].cols() + img.cols());
        dimg[3](i, j) = img(i - dimg[0].rows() + img.rows(),
                            j - dimg[0].cols() + img.cols());
        dimg[4](i, j) = img(i + (img.rows() - dimg[0].rows()) / 2,
                            j + (img.cols() - dimg[0].cols()) / 2);
      }
    Mat dres[5];
    for(int i = 0; i < 5; i ++)
      dres[i] = represent(dimg[i], depth - 1);
    Mat res(1 + dres[0].rows() * 5, size);
    res.row(0) = res0.row(0);
    for(int i = 0; i < 5; i ++)
      for(int j = 0; j < dres[i].rows(); j ++)
        res.row(1 + i * dres[i].rows() + j) = dres[i].row(j);
    return res;
  }
  return res0;
}


static inline bool whiteline(const string& s) {
  for(auto ss(s.begin()); ss < s.end(); ++ ss)
    if(! std::isspace(* ss) && *ss != '\n')
      return false;
  return true;
}

template <typename T> static inline bool loadstub(istream& input, const int& nmax, const int& ncolor, vector<SimpleMatrix<T> >& datas) {
  int i = 0, j = 0, k = 0;
  char buf;
  int  work = 0;
  bool mode = false;
  while(input.get(buf) && j < datas[0].rows()) {
    if('0' <= buf && buf <= '9') {
      work *= 10;
      work += buf - '0';
      mode  = true;
      continue;
    } else if(mode) {
      mode = false;
      datas[k](j, i) = T(work) / (T(nmax) + T(int(1)));
      work = 0;
      if(++ k >= ncolor) {
        if(++ i >= datas[0].cols()) {
          if(++ j >= datas[0].rows())
            return true;
          i = 0;
        }
        k = 0;
      }
    }
  }
  return true;
}

template <typename T> bool loadp2or3(vector<SimpleMatrix<T> >& data, istream& input) {
  string line;
  string line2;
  string line3;
  try {
    data.resize(3, SimpleMatrix<T>());
    getline(input, line);
    while((whiteline(line) || line[0] == '#') && getline(input, line) && !input.eof() && !input.bad()) ;
    getline(input, line2);
    while((whiteline(line2) || line2[0] == '#') && getline(input, line2) && !input.eof() && !input.bad()) ;
    getline(input, line3);
    while((whiteline(line3) || line3[0] == '#') && getline(input, line3) && !input.eof() && !input.bad()) ;
    istringstream iline2(line2);
    int w, h;
    iline2 >> w;
    iline2 >> h;
    if(line.size() < 2 || w <= 0 || h <= 0) {
      cerr << "unknown size." << endl;
      return false;
    }
    istringstream iline3(line3);
    int nmax;
    iline3 >> nmax;
    if(line[0] == 'P') {
      if(line[1] == '2') {
        data.resize(1);
        data[0] = SimpleMatrix<T>(h,w ).O();
        loadstub<T>(input, nmax, 1, data);
      } else if(line[1] == '3') {
        for(int i = 0; i < 3; i ++)
          data[i] = SimpleMatrix<T>(h, w).O();
        loadstub<T>(input, nmax, 3, data);
      } else {
        cerr << "unknown file type." << endl;
        return false;
      }
    } else {
      cerr << "unknown file type." << endl;
      return false;
    }
  } catch (...) {
    cerr << "Exception while reading." << endl;
    return false;
  }
  return true;
}

template <typename T> bool loadp2or3(vector<SimpleMatrix<T> >& data, const char* filename) {
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    if(! loadp2or3<T>(data, static_cast<istream&>(input))) {
      input.close();
      return false;
    }
    input.close();
  } else {
    cerr << "Unable to open file for read: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> bool savep2or3(const char* filename, const vector<SimpleMatrix<T> >& data, const int& depth = 65535) {
  ofstream output;
  output.open(filename);
  if(output.is_open()) {
    try {
      if(data.size() == 1)
        output << "P2" << "\n";
      else
        output << "P3" << "\n";
      output << data[0].cols() << " " << data[0].rows() << "\n";
      output << depth << "\n";
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if(data.size() == 1)
            output << min(int(depth), int(data[0](i, j) * (T(depth) + T(int(1)))) ) << "\n";
          else
            for(int k = 0; k < 3; k ++)
              output << min(int(depth), int(data[k](i, j) * (T(depth) + T(int(1)))) ) << "\n";
    } catch (...) {
      cerr << "An error has occured while writing file." << endl;
    }
    output.close();
  } else {
    cerr << "Unable to open file for write: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> static inline vector<vector<SimpleMatrix<T> > > normalize(const vector<vector<SimpleMatrix<T> > >& data, const T& upper = T(1)) {
  T MM(0), mm(0);
  bool fixed(false);
  for(int kk = 0; kk < data.size(); kk ++)
    for(int k = 0; k < data[kk].size(); k ++)
      for(int i = 0; i < data[kk][k].rows(); i ++)
        for(int j = 0; j < data[kk][k].cols(); j ++)
          if(! fixed || (isfinite(data[kk][k](i, j)) &&
               ! isinf(data[kk][k](i, j)) && ! isnan(data[kk][k](i, j)))) {
            if(! fixed)
              MM = mm = data[kk][k](i, j);
            else {
              MM = max(MM, data[kk][k](i, j));
              mm = min(mm, data[kk][k](i, j));
            }
            fixed = true;
          }
  if(MM == mm || ! fixed)
    return data;
  auto result(data);
  for(int kk = 0; kk < data.size(); kk ++)
    for(int k = 0; k < data[kk].size(); k ++)
      for(int i = 0; i < data[kk][k].rows(); i ++)
        for(int j = 0; j < data[kk][k].cols(); j ++) {
          if(isfinite(result[kk][k](i, j)) && ! isinf(data[kk][k](i, j)) && ! isnan(result[kk][k](i, j)))
            result[kk][k](i, j) -= mm;
          else
            result[kk][k](i, j)  = T(0);
          assert(T(0) <= result[kk][k](i, j) && result[kk][k](i, j) <= MM - mm);
          result[kk][k](i, j) *= upper / (MM - mm);
        }
  return result;
}

template <typename T> static inline vector<SimpleMatrix<T> > normalize(const vector<SimpleMatrix<T> >& data, const T& upper = T(1)) {
  vector<vector<SimpleMatrix<T> > > w;
  w.emplace_back(data);
  return normalize<T>(w, upper)[0];
}

template <typename T> static inline SimpleVector<T> normalize(const SimpleVector<T>& in, const T& upper = T(1)) {
  vector<vector<SimpleMatrix<T> > > w;
  w.resize(1);
  w[0].resize(1);
  w[0][0].resize(1, in.size());
  w[0][0].row(0) = in;
  return normalize<T>(w, upper)[0][0].row(0);
}

template <typename T> static inline vector<SimpleMatrix<T> > autoLevel(const vector<SimpleMatrix<T> >& data, const int& count = 0) {
  vector<T> res;
  res.reserve(data[0].rows() * data[0].cols() * data.size());
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        res.emplace_back(data[k](i, j));
  sort(res.begin(), res.end());
  auto result(data);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        result[k](i, j) = max(min(data[k](i, j), res[res.size() - count - 1]), res[count]);
  return result;
}

template <typename T> static inline SimpleVector<T> autoLevel(const SimpleVector<T>& data, const int& count = 0) {
  vector<SimpleMatrix<T> > b;
  b.resize(1);
  b[0].resize(1, data.size());
  b[0].row(0) = data;
  return autoLevel<T>(b, count)[0].row(0);
}

template <typename T> static inline vector<SimpleMatrix<T> > autoGamma(const vector<SimpleMatrix<T> >& data, const T& ratio = T(int(1)) / T(int(2)) ) {
  T r(int(0));
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++) {
        assert(T(int(0)) <= data[k](i, j) && data[k](i, j) <= T(int(1)) );
        r += log(data[k](i, j) + T(int(1)) / T(int(65536)) );
      }
  r /= T(int(data.size() * data[0].rows() * data[0].cols()));
  r  = log(ratio) / r;
  auto result(data);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        result[k](i, j) = max(T(int(0)), min(T(int(1)), pow(data[k](i, j) + T(int(1)) / T(int(65536)), r) ));
  return result;
}

template <typename T> static inline SimpleVector<T> autoGamma(const SimpleVector<T>& data, const T& r = T(int(1)) / T(int(2)) ) {
  vector<SimpleMatrix<T> > b;
  b.resize(1);
  b[0].resize(1, data.size());
  b[0].row(0) = data;
  return autoGamma<T>(b, r)[0].row(0);
}

template <typename T> static inline T getImgPt(const T& y, const T& h) {
  auto yy(y % (2 * h));
  if(yy < 0)
    yy = - yy;
  if(yy >= h)
    yy = h - (yy - h);
  return yy % h;
}

template <typename T> SimpleMatrix<T> rotate(const SimpleMatrix<T>& d, const T& 
theta) {
  assert(abs(theta) < atan(T(1)));
  const auto c(cos(theta));
  const auto s(sin(theta));
  const auto h0(abs(int(c * T(d.rows()) - s * T(d.cols()))));
  const auto h1(h0 + abs(int(s * T(d.cols()))) * 2);
  const auto w0(abs(int(s * T(d.rows()) + c * T(d.cols()))));
  const auto w1(w0 + abs(int(s * T(d.rows()))) * 2);
  SimpleMatrix<T> res(h0 < d.rows() ? h1 : h0,
                      w0 < d.cols() ? w1 : w0);
  const T offy(h0 < d.rows() ? abs(int(s * T(d.cols()))) : 0);
  const T offx(w0 < d.cols() ? abs(int(s * T(d.rows()))) : 0);
  res.O();
  const auto diag(int(sqrt(res.rows() * res.rows() +
                           res.cols() * res.cols())) + 1);
  for(int j = - diag; j < diag; j ++)
    for(int k = - diag; k < diag; k ++) {
      const int yy(c * T(j) - s * T(k) + offy);
      const int xx(s * T(j) + c * T(k) + offx);
      if(0 <= yy && yy < res.rows() &&
         0 <= xx && xx < res.cols()) {
        const auto dyy(getImgPt<int>(j, d.rows()));
        const auto dxx(getImgPt<int>(k, d.cols()));
        {
          res(yy, xx) = res(min(yy + 1, int(res.rows()) - 1), xx) =
            res(yy, min(xx + 1, int(res.cols()) - 1)) =
            res(min(yy + 1, int(res.rows()) - 1),
                min(xx + 1, int(res.cols()) - 1)) =
              d(dyy, dxx);
        }
      }
    }
  return res;
}

template <typename T> static inline SimpleMatrix<T> center(const SimpleMatrix<T>& dr, const SimpleMatrix<T>& d) {
  SimpleMatrix<T> res(d.rows(), d.cols());
  for(int i = 0; i < res.rows(); i ++)
    for(int j = 0; j < res.cols(); j ++)
      res(i, j) = dr(max(int(0), min(i + (dr.rows() - d.rows()) / 2, dr.rows() - 1)),
                     max(int(0), min(j + (dr.cols() - d.cols()) / 2, dr.cols() - 1)));
  return res;
}

// N.B. the raw P01 predictor is useless because of their sloppiness.
template <typename T> pair<vector<SimpleVector<T> >, vector<SimpleVector<T> > > predv(const vector<SimpleVector<T> >& in) {
  // N.B. we need to initialize p0 vector.
  SimpleVector<T> init(3);
  for(int i = 0; i < init.size(); i ++)
    init[i] = T(int(i));
  cerr << "Coherent: P0: " << P0maxRank0<T>().next(init) << endl;
  SimpleVector<T> seconds(in.size());
  seconds.O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < in.size(); i ++)  {
    seconds[i] = makeProgramInvariant<T>(in[i], - T(int(1)), true).second;
  }
  // N.B. we need Ppersistent<..., P01...> with maximum range.
  const int istat(5 * 5 - 4 + 2);
  vector<SimpleVector<T> > p;
  vector<SimpleVector<T> > q;
  p.resize(1, SimpleVector<T>(in[0].size()).O());
  q.resize(1, SimpleVector<T>(in[0].size()).O());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int j = 0; j < in[0].size(); j ++) {
    cerr << j << " / " << in[0].size() << endl;
    idFeeder<T> buf(in.size());
    for(int i = 0; i < in.size(); i ++)
      buf.next(makeProgramInvariantPartial<T>(in[i][j], seconds[i], true));
    p[0][j] = PpersistentOnce<T, P01<T, true> >().next(buf.res, istat);
    q[0][j] = PpersistentOnce<T, P01<T, true> >().next(buf.res.reverse(), istat);
  }
  // N.B. We believe predictor for each pixel context strong and
  //      only the 0 <= x in R^n each orthogonality condition.
  //      We should revertProgramInvariant in general, however, we bet
  //      each pixel context strong.
  const auto nseconds(sqrt(seconds.dot(seconds)));
  const auto rp(PpersistentOnce<T, P01<T, true> >().next(seconds / nseconds, istat) * nseconds);
  const auto rq(PpersistentOnce<T, P01<T, true> >().next(seconds.reverse() / nseconds, istat) * nseconds);
  const auto p0(makeProgramInvariant<T>(normalize<T>(p[0]), - T(int(1)), true).first);
  const auto q0(makeProgramInvariant<T>(normalize<T>(q[0]), - T(int(1)), true).first);
  for(int j = 0; j < p[0].size(); j ++) {
    p[0][j] = revertProgramInvariant<T>(make_pair(p0[j], rp), true);
    q[0][j] = revertProgramInvariant<T>(make_pair(q0[j], rq), true);
  }
  return make_pair(move(p), move(q));
}

template <typename T> pair<vector<vector<SimpleVector<T> > >, vector<vector<SimpleVector<T> > > > predVec(const vector<vector<SimpleVector<T> > >& in0) {
  assert(in0.size() && in0[0].size() && in0[0][0].size());
  vector<SimpleVector<T> > in;
  in.resize(in0.size());
  for(int i = 0; i < in0.size(); i ++) {
    assert(in0[i].size() == in0[0].size() &&
           in0[i][0].size() == in0[0][0].size());
    in[i].resize(in0[i].size() * in0[i][0].size());
    for(int j = 0; j < in0[i].size(); j ++) {
      assert(in0[i][0].size() == in0[i][j].size());
      in[i].setVector(j * in0[i][0].size(), in0[i][j]);
    }
  }
  const auto p(predv<T>(in));
  pair<vector<vector<SimpleVector<T> > >, vector<vector<SimpleVector<T> > > > res;
  res.first.resize(p.first.size());
  res.second.resize(p.second.size());
  for(int i = 0; i < res.first.size(); i ++) {
    res.first[i].resize(in0[0].size());
    res.second[i].resize(in0[0].size());
    for(int j = 0; j < in0[0].size(); j ++) {
      res.first[i][j] =
        p.first[i].subVector(in0[0][0].size() * j, in0[0][0].size());
      res.second[i][j] =
        p.second[i].subVector(in0[0][0].size() * j, in0[0][0].size());
    }
  }
  return res;
}

template <typename T> pair<vector<vector<SimpleMatrix<T> > >, vector<vector<SimpleMatrix<T> > > > predMat(const vector<vector<SimpleMatrix<T> > >& in0) {
  assert(in0.size() && in0[0].size() && in0[0][0].rows() && in0[0][0].cols());
  vector<SimpleVector<T> > in;
  in.resize(in0.size());
  for(int i = 0; i < in0.size(); i ++) {
    assert(in0[i].size() == in0[0].size());
    in[i].resize(in0[i].size() * in0[i][0].rows() * in0[i][0].cols());
    for(int j = 0; j < in0[i].size(); j ++) {
      assert(in0[i][j].rows() == in0[0][0].rows() &&
             in0[i][j].cols() == in0[0][0].cols());
      for(int k = 0; k < in0[i][j].rows(); k ++)
        in[i].setVector(j * in0[i][0].rows() * in0[i][0].cols() +
                        k * in0[i][0].cols(), in0[i][j].row(k));
    }
  }
  const auto p(predv<T>(in));
  pair<vector<vector<SimpleMatrix<T> > >, vector<vector<SimpleMatrix<T> > > > res;
  res.first.resize(p.first.size());
  res.second.resize(p.second.size());
  for(int i = 0; i < res.first.size(); i ++) {
    res.first[i].resize(in0[0].size());
    res.second[i].resize(in0[0].size());
    for(int j = 0; j < res.first[i].size(); j ++) {
      res.first[i][j].resize(in0[0][0].rows(), in0[0][0].cols());
      res.second[i][j].resize(in0[0][0].rows(), in0[0][0].cols());
      for(int k = 0; k < in0[0][0].rows(); k ++) {
        res.first[i][j].row(k) =
          p.first[i].subVector(
            j * in0[0][0].rows() * in0[0][0].cols() + k * in0[0][0].cols(),
            in0[0][0].cols());
        res.second[i][j].row(k) =
          p.second[i].subVector(
            j * in0[0][0].rows() * in0[0][0].cols() + k * in0[0][0].cols(),
            in0[0][0].cols());
      }
    }
  }
  return res;
}

template <typename T> pair<vector<SimpleSparseTensor<T> >, vector<SimpleSparseTensor<T> > > predSTen(const vector<SimpleSparseTensor<T> >& in0, const vector<int>& idx) {
  assert(idx.size() && in0.size());
  // N.B. we don't do input scaling.
  // N.B. the data we target is especially string stream corpus.
  //      they are incontinuous one, so complementing with continuous stream
  //      shouldn't improve outputs.
  vector<SimpleVector<T> > in;
  vector<pair<int, pair<int, int> > > attend;
  in.resize(in0.size());
  attend.reserve(idx.size() * idx.size() * idx.size());
  for(int i = 0; i < idx.size(); i ++)
    for(int j = 0; j < idx.size(); j ++)
      for(int k = 0; k < idx.size(); k ++) {
        for(int ii = 0; ii < in0.size(); ii ++)
          if(in0[ii][idx[i]][idx[j]][idx[k]] != T(int(0)))
            goto next;
        continue;
       next:
        attend.emplace_back(make_pair(i, make_pair(j, k)));
      }
  sort(attend.begin(), attend.end());
  for(int i = 0; i < in0.size(); i ++) {
    in[i].resize(attend.size());
    for(int j = 0, cnt = 0; j < idx.size(); j ++)
      for(int k = 0; k < idx.size(); k ++)
        for(int m = 0; m < idx.size(); m ++)
          if(binary_search(attend.begin(), attend.end(),
              make_pair(j, make_pair(k, m))))
            in[i][cnt ++] =
              (in0[i][idx[j]][idx[k]][idx[m]] + T(int(1))) / T(int(2));
  }
  auto p(predv<T>(in));
  in.resize(0);
  pair<vector<SimpleSparseTensor<T> >, vector<SimpleSparseTensor<T> > > res;
  res.first.resize(p.first.size());
  res.second.resize(p.second.size());
  for(int i = 0; i < res.first.size(); i ++)
    for(int j = 0, cnt = 0; j < idx.size(); j ++)
      for(int k = 0; k < idx.size(); k ++)
        for(int m = 0; m < idx.size(); m ++)
          if(binary_search(attend.begin(), attend.end(),
               make_pair(j, make_pair(k, m)))) {
            res.first[i][idx[j]][idx[k]][idx[m]] = p.first[i][cnt] * T(int(2)) - T(int(1));
            res.second[i][idx[j]][idx[k]][idx[m]] = p.second[i][cnt ++] * T(int(2)) - T(int(1));
          }
  return res;
}

template <typename T> static inline vector<SimpleMatrix<T> > rgb2xyz(const vector<SimpleMatrix<T> >& rgb) {
  // CIE 1931 XYZ from wikipedia.org
  SimpleMatrix<T> mRGB2XYZ(3, 3);
  mRGB2XYZ(0, 0) = T(49000);
  mRGB2XYZ(0, 1) = T(31000);
  mRGB2XYZ(0, 2) = T(20000);
  mRGB2XYZ(1, 0) = T(17697);
  mRGB2XYZ(1, 1) = T(81240);
  mRGB2XYZ(1, 2) = T( 1063);
  mRGB2XYZ(2, 0) = T(0);
  mRGB2XYZ(2, 1) = T( 1000);
  mRGB2XYZ(2, 2) = T(99000);
  mRGB2XYZ /= T(17697);
  assert(rgb.size() == 3);
  assert(rgb[0].rows() == rgb[1].rows() && rgb[1].rows() == rgb[2].rows());
  assert(rgb[0].cols() == rgb[1].cols() && rgb[1].cols() == rgb[2].cols());
  auto xyz(rgb);
  xyz[0] = rgb[0] * mRGB2XYZ(0, 0) + rgb[1] * mRGB2XYZ(0, 1) + rgb[2] * mRGB2XYZ(0, 2);
  xyz[1] = rgb[0] * mRGB2XYZ(1, 0) + rgb[1] * mRGB2XYZ(1, 1) + rgb[2] * mRGB2XYZ(1, 2);
  xyz[2] = rgb[0] * mRGB2XYZ(2, 0) + rgb[1] * mRGB2XYZ(2, 1) + rgb[2] * mRGB2XYZ(2, 2);
  assert(xyz.size() == 3);
  assert(xyz[0].rows() == xyz[1].rows() && xyz[1].rows() == xyz[2].rows());
  assert(xyz[0].cols() == xyz[1].cols() && xyz[1].cols() == xyz[2].cols());
  return xyz;
}

template <typename T> static inline vector<SimpleMatrix<T> > xyz2rgb(const vector<SimpleMatrix<T> >& xyz) {
  // CIE 1931 XYZ from wikipedia.org
  SimpleMatrix<T> mRGB2XYZ(3, 3);
  mRGB2XYZ(0, 0) = T(49000);
  mRGB2XYZ(0, 1) = T(31000);
  mRGB2XYZ(0, 2) = T(20000);
  mRGB2XYZ(1, 0) = T(17697);
  mRGB2XYZ(1, 1) = T(81240);
  mRGB2XYZ(1, 2) = T( 1063);
  mRGB2XYZ(2, 0) = T(0);
  mRGB2XYZ(2, 1) = T( 1000);
  mRGB2XYZ(2, 2) = T(99000);
  mRGB2XYZ /= T(17697);
  const auto mXYZ2RGB(mRGB2XYZ.inverse());
  assert(xyz.size() == 3);
  assert(xyz[0].rows() == xyz[1].rows() && xyz[1].rows() == xyz[2].rows());
  assert(xyz[0].cols() == xyz[1].cols() && xyz[1].cols() == xyz[2].cols());
  auto rgb(xyz);
  rgb[0] = xyz[0] * mXYZ2RGB(0, 0) + xyz[1] * mXYZ2RGB(0, 1) + xyz[2] * mXYZ2RGB(0, 2);
  rgb[1] = xyz[0] * mXYZ2RGB(1, 0) + xyz[1] * mXYZ2RGB(1, 1) + xyz[2] * mXYZ2RGB(1, 2);
  rgb[2] = xyz[0] * mXYZ2RGB(2, 0) + xyz[1] * mXYZ2RGB(2, 1) + xyz[2] * mXYZ2RGB(2, 2);
  assert(rgb.size() == 3);
  assert(rgb[0].rows() == rgb[1].rows() && rgb[1].rows() == rgb[2].rows());
  assert(rgb[0].cols() == rgb[1].cols() && rgb[1].cols() == rgb[2].cols());
  return rgb;
}

static const vector<int>& pnTinySingle(const int& upper = 1) {
  static vector<int> pn;
  if(! pn.size()) pn.emplace_back(2);
  pn.reserve(upper);
  for(int i = pn.size(); i < upper; i ++) {
    for(int j = pn[pn.size() - 1] + 1; 0 <= j; j ++) {
      for(int k = 0; k < pn.size(); k ++)
        if(! (j % pn[k])) goto next_pn;
      pn.emplace_back(j);
      break;
     next_pn:
      ;
    }
  }
  return pn;
}

#define _SIMPLELIN_
#endif


/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/* BSD 3-Clause License:
 * Copyright (c) 2013 - 2021, kazunobu watatsu.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(_SIMPLELIN_)

using std::move;
using std::swap;
using std::vector;
using std::sort;
using std::pair;
using std::make_pair;
using std::min;
using std::max;
using std::map;
using std::stringstream;
using std::endl;

using std::cerr;
using std::flush;

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
  assert(size > 0);
  this->entity = new T[size];
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
  SimpleVector<T> res(*this);
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
  SimpleVector<T> res(*this);
  return res -= other;
}

template <typename T> inline SimpleVector<T>& SimpleVector<T>::operator -= (const SimpleVector<T>& other) {
  return *this += - other;
}

template <typename T> inline SimpleVector<T> SimpleVector<T>::operator * (const T& other) const {
  SimpleVector<T> res(*this);
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
  SimpleVector<T> res(*this);
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
  for(int i = 0; i < esize; i ++)
    result.entity[i] = U(entity[i].real());
  return result;
}

template <typename T> template <typename U> inline SimpleVector<U> SimpleVector<T>::imag() const {
  SimpleVector<U> result(esize);
  for(int i = 0; i < esize; i ++)
    result.entity[i] = U(entity[i].imag());
  return result;
}

template <typename T> template <typename U> inline SimpleVector<U> SimpleVector<T>::cast() const {
  SimpleVector<U> result(esize);
  for(int i = 0; i < esize; i ++)
    result.entity[i] = U(entity[i]);
  return result;
}

template <typename T> inline const int& SimpleVector<T>::size() const {
  return esize;
}

template <typename T> inline void SimpleVector<T>::resize(const int& size) {
  assert(size > 0);
  if(size != esize) {
    esize = size;
    if(entity)
      delete[] entity;
    entity = new T[esize];
  }
  return;
}

template <typename T> std::ostream& operator << (std::ostream& os, const SimpleVector<T>& v) {
  SimpleVector<std::string> buf(v.size());
  int M(0);
  for(int i = 0; i < v.size(); i ++) {
    stringstream ss;
    ss << v[i];
    buf[i] = ss.str();
    M = max(int(buf[i].size()), M);
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
    if(c == ' ' || c == '\t' || c == '[' || c == ',' || c == ']') continue;
    is >> v[i ++];
  }
  if(i < v.size()) {
    cerr << "XXX SimpleVector<T>::operator >>" << flush;
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
  inline       T                determinant() const;
  inline       SimpleMatrix<T>  inverse() const;
  inline       SimpleVector<T>  solve(SimpleVector<T> other) const;
  inline       SimpleVector<T>  projectionPt(const SimpleVector<T>& other) const;
  inline       SimpleMatrix<T>& fillP(const vector<int>& idx);
  inline       SimpleMatrix<T>  QR() const;
  inline       SimpleMatrix<T>  SVD() const;
  inline       pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SVD(const SimpleMatrix<T>& src) const;
  inline       vector<int>      innerFix(const SimpleMatrix<T>& A, vector<pair<T, int> >& fidx);
  inline       SimpleVector<T>  inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const;
  template <typename U> inline SimpleMatrix<U> real() const;
  template <typename U> inline SimpleMatrix<U> imag() const;
  template <typename U> inline SimpleMatrix<U> cast() const;
  inline const int& rows() const;
  inline const int& cols() const;
  inline       void resize(const int& rows, const int& cols);
  num_t        epsilon;

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
  epsilon = num_t(1) >> int64_t(mybits - 1);
#else
  epsilon = std::numeric_limits<num_t>::epsilon();
#endif
  return;
}

template <typename T> inline SimpleMatrix<T>::SimpleMatrix(const int& rows, const int& cols) {
  assert(rows > 0 && cols > 0);
  entity = new SimpleVector<T>[rows];
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < rows; i ++)
    entity[i].resize(cols);
  erows = rows;
  ecols = cols;
#if defined(_FLOAT_BITS_)
  epsilon = num_t(1) >> int64_t(mybits - 1);
#else
  epsilon = std::numeric_limits<num_t>::epsilon();
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
  SimpleMatrix<T> res(*this);
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
  SimpleMatrix<T> res(*this);
  return res -= other;
}

template <typename T> inline SimpleMatrix<T>& SimpleMatrix<T>::operator -= (const SimpleMatrix<T>& other) {
  return *this += - other;
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::operator * (const T& other) const {
  SimpleMatrix<T> res(*this);
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
  SimpleMatrix<T> derived(other.transpose());
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
  SimpleMatrix<T> res(*this);
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
    if(entity)
      delete[] entity;
    if(other.erows)
      entity = new SimpleVector<T>[other.erows];
    else
      entity = NULL;
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

template <typename T> inline T SimpleMatrix<T>::determinant() const {
  assert(0 <= erows && 0 <= ecols && erows == ecols);
  T det(1);
  SimpleMatrix<T> work(*this);
  for(int i = 0; i < erows; i ++) {
    int xchg = i;
    for(int j = i + 1; j < erows; j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    SimpleVector<T> buf(work.entity[i]);
    work.entity[i]    = work.entity[xchg];
    work.entity[xchg] = buf;
    const SimpleVector<T>& ei(work.entity[i]);
    const T&               eii(ei[i]);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = i + 1; j < erows; j ++) {
      const T ratio(work.entity[j][i] / eii);
      work.entity[j] -= ei * ratio;
    }
    det *= ei[i];
  }
  return isfinite(det) ? det : T(0);
}

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::inverse() const {
  // XXX: extremely slow implementation.
  SimpleMatrix<T> result(erows, ecols);
  for(int i = 0; i < result.rows(); i ++)
    for(int j = 0; j < result.cols(); j ++)
      result(i, j) = i == j ? T(1) : T(0);
  for(int i = 0; i < result.cols(); i ++)
    result.setCol(i, solve(result.col(i)));
  return result;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::solve(SimpleVector<T> other) const {
  assert(0 <= erows && 0 <= ecols && erows == ecols && entity && erows == other.size());
  SimpleMatrix<T> work(*this);
  for(int i = 0; i < erows; i ++) {
    int xchg = i;
    for(int j = i + 1; j < erows; j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    swap(work.entity[i], work.entity[xchg]);
    swap(other[i], other[xchg]);
    const SimpleVector<T>& ei(work.entity[i]);
    const T&               eii(ei[i]);
    if(epsilon < eii * eii) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = i + 1; j < erows; j ++) {
        const T ratio(work.entity[j][i] / eii);
        work.entity[j] -= ei       * ratio;
        other[j]       -= other[i] * ratio;
      }
    }
  }
  for(int i = erows - 1; 0 <= i; i --) {
    if(work.entity[i][i] == T(0)) continue;
    const T buf(other[i] / work.entity[i][i]);
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
  for(int j = 0; j < this->cols() && ii < idx.size(); j ++) {
    SimpleVector<T> ek(this->cols());
    for(int k = 0; k < this->cols(); k ++)
      ek[k] = T(j == k ? 1 : 0);
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
  for(int i = 0; i < Q.rows(); i ++)
    for(int j = 0; j < Q.cols(); j ++)
      Q(i, j) = T(0);
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
  if(this->cols() < this->rows()) {
    auto res((* this) * this->transpose().SVD());
    vector<int> residue;
    residue.reserve(res.cols());
    for(int i = 0; i < res.cols(); i ++) {
      const auto r2(res.col(i).dot(res.col(i)));
      if(epsilon < r2)
        res.setCol(i, res.col(i) / sqrt(r2));
      else
        residue.emplace_back(i);
    }
    if(residue.size())
      return res.transpose().fillP(residue).transpose();
    return res;
  }
  for(int i = 0; i < this->rows(); i ++)
    for(int j = 0; j < this->cols(); j ++)
      assert(isfinite((*this)(i, j)));
  auto s((*this) * this->transpose());
  for(int i = 0; i < s.rows(); i ++)
    for(int j = 0; j < s.cols(); j ++)
      assert(isfinite(s(i, j)));
  SimpleMatrix<T> left(s.rows(), s.rows());
  for(int i = 0; i < left.rows(); i ++)
    for(int j = 0; j < left.cols(); j ++)
      left(i, j) = T(i == j ? 1 : 0);
  for(int ii = 0; ii < s.rows(); ii ++) {
    // find eigen max on working matrix:
    SimpleMatrix<T> work(s.rows() - ii, s.rows() - ii);
    for(int i = 0; i < work.rows(); i ++)
      for(int j = 0; j < work.cols(); j ++)
        work(i, j) = s(ii + i, ii + j);
    SimpleVector<T>   d(work.rows());
    SimpleVector<int> idx(work.rows());
    // LDLt:
    for(int i = 0; i < work.rows(); i ++) {
      int xchg = i;
      for(int j = i + 1; j < work.rows(); j ++)
        if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
          xchg = j;
      swap(work.row(i), work.row(idx[i] = xchg));
      const auto& ei(work.row(i));
      const auto& eii(ei[i]);
      d[i] = eii;
      if(epsilon < eii * eii) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int j = i + 1; j < work.rows(); j ++) {
          const auto ratio(work(j, i) / eii);
          assert(isfinite(ratio));
          work.row(j) -= ei * ratio;
        }
      }
    }
    for(int i = 0; i < idx.size(); i ++) {
      swap(work.row(idx[idx.size() - 1 - i]), work.row(idx.size() - 1 - i));
      swap(d[idx[idx.size() - 1 - i]], d[idx.size() - 1 - i]);
    }
    auto eMidx(0);
    for(int i = 1; i < d.size(); i ++)
      if(abs(d[eMidx]) < abs(d[i]))
        eMidx = i;
    SimpleVector<T> ev(work.rows());
    for(int i = 0; i < ev.size(); i ++)
      ev[i] = i == eMidx ? T(1) : T(0);
    ev    = work.solve(ev);
    ev   /= sqrt(ev.dot(ev));
    work *= T(0);
    work.row(0) = ev;
    vector<int> fp;
    fp.reserve(work.rows() - 1);
    for(int i = 1; i < work.rows(); i ++)
      fp.emplace_back(i);
    work.fillP(fp);
    auto mul(left);
    for(int i = 0; i < mul.rows(); i ++)
      for(int j = 0; j < mul.cols(); j ++)
        mul(i, j) = i < ii || j < ii ? T(i == j ? 1 : 0) :
          work(i - ii, j - ii);
    left = left * mul;
    s    = mul.transpose() * s * mul;
  }
  for(int i = 0; i < left.rows(); i ++) {
    int jj(i);
    for(int j = i + 1; j < left.cols(); j ++)
      if(abs(left(jj, i)) < abs(left(j, i)))
        jj = j;
    swap(left.row(i), left.row(jj));
  }
  return left;
}

// N.B. for full rank.
template <typename T> inline pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SimpleMatrix<T>::SVD(const SimpleMatrix<T>& src) const {
  // refered from : https://en.wikipedia.org/wiki/Generalized_singular_value_decomposition .
  assert(this->cols() == src.cols());
  SimpleMatrix<T> C(this->rows() + src.rows(), this->cols());
  for(int i = 0; i < this->rows(); i ++)
    C.row(i) = this->row(i);
  for(int i = 0; i < src.rows(); i ++)
    C.row(i + this->rows()) = src.row(i);
  const auto P(C.SVD());
  SimpleVector<T> d(this->cols());
        auto Qt(P.transpose() * C);
  for(int i = 0; i < d.size(); i ++)
    Qt.row(i) /= (d[i] = sqrt(Qt.row(i).dot(Qt.row(i))));
  const auto D(P.transpose() * C * Qt.transpose());
  SimpleMatrix<T> P1(this->rows(), d.size());
  SimpleMatrix<T> P2(src.rows(), d.size());
  for(int i = 0; i < P1.rows(); i ++)
    P1.row(i) = P.row(i);
  for(int i = 0; i < P2.rows(); i ++)
    P2.row(i) = P.row(i + P1.rows());
  auto U1(P1.SVD());
  auto Wt(U1.transpose() * P1);
  for(int i = 0; i < Wt.rows(); i ++)
    Wt.row(i) /= sqrt(Wt.row(i).dot(Wt.row(i)));
  auto U2(P2 * Wt.transpose());
  for(int i = 0; i < U2.cols(); i ++) {
    const auto u2i(U2.col(i));
    U2.setCol(u2i / sqrt(u2i.dot(u2i)));
  }
  return make_pair(make_pair(move(U1), move(U2)), (Wt * D).transpose().QR().transpose() * Qt);
}

template <typename T> template <typename U> inline SimpleMatrix<U> SimpleMatrix<T>::real() const {
  assert(0 < erows && 0 < ecols);
  SimpleMatrix<U> res(erows, ecols);
  for(int i = 0; i < erows; i ++)
    for(int j = 0; j < ecols; j ++)
      res(i, j) = U(entity[i][j].real());
  return res;
}

template <typename T> template <typename U> inline SimpleMatrix<U> SimpleMatrix<T>::imag() const {
  assert(0 < erows && 0 < ecols);
  SimpleMatrix<U> res(erows, ecols);
  for(int i = 0; i < erows; i ++)
    for(int j = 0; j < ecols; j ++)
      res(i, j) = U(entity[i][j].imag());
  return res;
}

template <typename T> template <typename U> inline SimpleMatrix<U> SimpleMatrix<T>::cast() const {
  assert(0 < erows && 0 < ecols);
  SimpleMatrix<U> res(erows, ecols);
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
  assert(rows > 0 && cols > 0);
  if(rows != erows) {
    erows = rows;
    if(entity)
      delete[] entity;
    entity = new SimpleVector<T>[erows];
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

template <typename T> inline vector<int> SimpleMatrix<T>::innerFix(const SimpleMatrix<T>& A, vector<pair<T, int> >& fidx) {
  // N.B. we now have |[A -bb] [x t]| <= 1 condition.
  // N.B. there's no difference |[A - bb] [x t]|^2 <= 1 condition in this.
  //      but not with mixed condition.
  const auto R((*this) * A);
  SimpleVector<T> one(this->cols());
  for(int i = 0; i < one.size(); i ++)
    one[i] = T(1);
  // we now have: Q [R [x t] ] <= {0, 1}^m cond.
  const auto on(projectionPt(one));
  fidx.reserve(fidx.size() + on.size());
  for(int i = 0; i < on.size(); i ++)
    fidx.emplace_back(make_pair(abs(on[i]), i));
  sort(fidx.begin(), fidx.end());
  vector<int> fix;
  fix.reserve(this->rows());
  // sort by: |<Q^t(1), q_k>|, we subject to minimize each, to do this,
  //   maximize minimum q_k orthogonality.
  for(int idx = 0; fix.size() < this->rows() - 1 && idx < fidx.size(); idx ++) {
    const auto& iidx(fidx[idx].second);
    const auto  orth(this->col(iidx));
    const auto  n2(orth.dot(orth));
    if(n2 <= epsilon)
      continue;
    fix.emplace_back(fidx[idx].second);
    // N.B. O(mn) can be writed into O(lg m + lg n) in many core cond.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < this->cols(); j ++)
      this->setCol(j, this->col(j) - orth * this->col(j).dot(orth) / n2);
  }
  // N.B. now we have fix indexes that to be A_fix [x 1] * t == 0.
  assert(fix.size() == this->rows() - 1);
  return fix;
}

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const {
  static const auto ee(pow(epsilon, T(1) / T(8)));
  assert(this->rows() == bl.size() && this->rows() == bu.size() &&
         0 < this->cols() && 0 < this->rows());
  // bu - bb == A, bl - bb == - A <=> bu - bl == 2 A. 
  const auto bb((bu + bl) / T(2));
  const auto upper(bu - bb);
  SimpleMatrix<T> A(this->rows() * 2 - 1 + this->cols() * 2, this->cols() + 1);
  vector<pair<T, int> > fidx;
  fidx.reserve(this->rows());
  for(int i = 0; i < this->rows(); i ++) {
    for(int j = 0; j < this->cols(); j ++)
      A(i, j) = (*this)(i, j);
    A(i, this->cols()) = - bb[i];
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
    const auto A2(A.row(i).dot(A.row(i)));
    assert(isfinite(A2) && ! isnan(A2));
    if(this->rows() - 1 <= i) break;
    A.row(i + this->rows()) = - A.row(i);
  }
  // this is tricky but one of them is fixed if it is needed.
  for(int i = this->rows() * 2 - 1; i < A.rows(); i ++) {
    const auto ii(i - (this->rows() * 2 - 1));
    for(int j = 0; j < this->cols(); j ++)
      A(i, j) = T(j == ii / 2 ? (ii & 1 ? - 1 : 1) : 0);
    A(i, this->cols()) = T(ii & 1 ? 1 : - 1);
    A.row(i) /= ee;
  }
  // N.B. we now have |[A -bb] [x t]| <= 1 condition.
  // N.B. there's no difference |[A - bb] [x t]|^2 <= 1 condition in this.
  //      but not with mixed condition.
        auto Pt(A.QR());
  const auto fix(Pt.innerFix(A, fidx));
  // N.B. now we have fix indexes that to be A_fix [x 1] * t == 0.
  assert(fix.size() == this->cols());
  SimpleMatrix<T> Afix(fix.size(), fix.size());
  SimpleVector<T> f(Afix.rows());
  for(int i = 0; i < fix.size(); i ++) {
    for(int j = 0; j < Afix.cols(); j ++)
      Afix(i, j) = A(fix[i], j);
    f[i] = - A(fix[i], Afix.cols());
  }
  return Afix.solve(f);
}

template <typename T> std::ostream& operator << (std::ostream& os, const SimpleMatrix<T>& v) {
  SimpleMatrix<std::string> buf(v.rows(), v.cols());
  int M(0);
  for(int i = 0; i < v.rows(); i ++)
    for(int j = 0; j < v.cols(); j ++) {
      stringstream ss;
      ss << v(i, j);
      buf(i, j) = ss.str();
      M = max(int(buf(i, j).size()), M);
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
    if(c == ' ' || c == '\t' || c == '(') continue;
    else {
      is.unget();
      break;
    }
  }
  int r, c;
  is >> r;
  if(r <= 0) return is;
  while(! is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == ',') continue;
    else {
      is.unget();
      break;
    }
  }
  is >> c;
  while(! is.eof() && ! is.bad()) {
    const auto c(is.get());
    if(c == ' ' || c == '\t' || c == ')') continue;
    else {
      is.unget();
      break;
    }
  }
  if(c <= 0) return is;
  v.resize(r, c);
  int i(0);
  int j(0);
  for( ; i < v.rows() && ! is.eof() && ! is.bad(); i ++) {
    while(! is.eof() && ! is.bad()) {
      const auto c(is.get());
      if(c == ' ' || c == '\t' || c == '[') continue;
      else {
        is.unget();
        break;
      }
    }
    for( ; j < v.cols() && ! is.eof() && ! is.bad(); j ++) {
      is >> v(i, j);
      if(v.cols() - 1 <= j) break;
      while(! is.eof() && ! is.bad()) {
        const auto c(is.get());
        if(c == ' ' || c == '\t' || c == ',') continue;
        else {
          is.unget();
          break;
        }
      }
    }
    while(! is.eof() && ! is.bad()) {
      const auto c(is.get());
      if(c == ' ' || c == '\t' || c == ']') continue;
      else {
        is.unget();
        break;
      }
    }
  }
  if(i < v.rows() || j < v.cols()) {
    cerr << "XXX SimpleMatrix<T>::operator >>" << flush;
    for( ; i < v.rows(); i ++)
      for( ; j < v.cols(); j ++)
        v(i, j) = T(0);
  }
  return is;
}


// some non class functions:
template <typename T> const SimpleMatrix<complex<T> >& dft(const int& size0) {
  const auto size(abs(size0));
  if(! size) {
    const static SimpleMatrix<complex<T> > m0;
    return m0;
  }
  static vector<SimpleMatrix<complex<T> > > adft;
  static vector<SimpleMatrix<complex<T> > > aidft;
  if(adft.size() <= size)
    adft.resize(size + 1, SimpleMatrix<complex<T> >());
  if(aidft.size() <= size)
    aidft.resize(size + 1, SimpleMatrix<complex<T> >());
  auto& edft(  adft[size]);
  auto& eidft(aidft[size]);
  if(edft.rows() != size || edft.cols() != size) {
    edft.resize( size, size);
    eidft.resize(size, size);
    static const auto Pi(T(4) * atan2(T(1), T(1)));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < edft.rows(); i ++) {
      for(int j = 0; j < edft.cols(); j ++) {
        const auto theta(- T(2) * Pi * T(i) * T(j) / T(edft.rows()));
        edft( i, j) = complex<T>(cos(  theta), sin(  theta));
        eidft(i, j) = complex<T>(cos(- theta), sin(- theta)) / complex<T>(T(size));
      }
    }
  }
  return size0 < 0 ? eidft : edft;
}

template <typename T> const SimpleMatrix<T>& diff(const int& size0) {
  const auto size(abs(size0));
  if(! size) {
    static const SimpleMatrix<T> m0;
    return m0;
  }
  static vector<SimpleMatrix<T> > D;
  static vector<SimpleMatrix<T> > I;
  if(D.size() <= size)
    D.resize(size + 1, SimpleMatrix<T>());
  if(I.size() <= size)
    I.resize(size + 1, SimpleMatrix<T>());
  auto& dd(D[size]);
  auto& ii(I[size]);
  if(dd.rows() != size || dd.cols() != size) {
    auto DD(dft<T>(size));
    auto II(dft<T>(size));
    static const auto Pi(T(4) * atan2(T(1), T(1)));
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
    for(int i = 0; i < DD.rows(); i ++)
      DD.row(i) *= - complex<T>(T(0), T(2) * Pi * T(i) / T(DD.rows()));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
    for(int i = 1; i < II.rows(); i ++)
      II.row(i) /= - complex<T>(T(0), T(2) * Pi * T(i) / T(DD.rows()) / Pi);
    dd = (dft<T>(-size) * DD).template real<T>() / Pi;
    ii = (dft<T>(-size) * II).template real<T>();
  }
  return size0 < 0 ? ii : dd;
}

template <typename T> inline SimpleVector<T> taylor(const int& size, const T& step) {
  const int  step00(max(int(0), min(size - 1, int(floor(step)))));
  const auto residue0(step - T(step00));
  const auto step0(step00 == size - 1 || abs(residue0) <= T(1) / T(2) ? step00 : step00 + 1);
  const auto residue(step - T(step0));
  SimpleVector<T> res(size);
  for(int i = 0; i < res.size(); i ++)
    res[i] = T(i == step0 ? 1 : 0);
  if(residue == T(0)) return res;
  const auto& D(diff<T>(size));
        auto  dt(D.col(step0) * residue);
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
    dt   = D * dt * residue / T(i);
  }
  return res;
}

template <typename T> SimpleVector<T> linearInvariant(const vector<SimpleVector<T> >& in) {
  SimpleMatrix<T> A(in.size() * 2 - 1, in[0].size());
  assert(in[0].size() <= in.size());
  SimpleVector<T> fvec(A.cols());
  SimpleVector<T> one(A.rows());
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < fvec.size(); i ++)
    fvec[i] = T(0);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < A.rows(); i ++)
    one[i] = T(1);
  for(int i = 0; i < in.size(); i ++) {
    assert(in[i].size() == in[0].size());
    for(int j = 0; j < in[0].size(); j ++)
      A(i, j) = in[i][j];
    assert(isfinite(A.row(i).dot(A.row(i))));
    if(i + in.size() < A.rows())
      A.row(i + in.size()) = - A.row(i);
  }
        auto Pt(A.QR());
  const auto R(Pt * A);
  vector<pair<T, int> > sute;
  Pt.innerFix(A, sute);
  fvec = R.solve(Pt * one);
  return isfinite(fvec.dot(fvec)) ? fvec : fvec * T(0);
}

// N.B. please refer bitsofcotton/randtools.
template <typename T> SimpleVector<T> makeProgramInvariant(const SimpleVector<T>& in, const int& complexity, const T& index = T(1)) {
  SimpleVector<T> res(in.size() + (T(0) <= index ? 2 : 1));
  for(int i = 0; i < in.size(); i ++) {
    assert(- T(1) <= in[i] && in[i] <= T(1) && isfinite(in[i]));
    res[i] = tan((in[i] + T(1)) / T(2) * atan2(T(1), T(1)));
  }
  res[in.size()] = T(1);
  if(T(0) <= index)
    res[in.size() + 1] = index;
  T pd(0);
  for(int i = 0; i < res.size(); i ++)
    pd += log(abs(res[i]));
  if(0 <= complexity)
    res *= exp((complexity ? T(complexity) * pd : - pd) / T(res.size()));
  return res;
}


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
  // already initialized by compiler, initializer doesn't need this.
//  entity = map<int, T>();
  return;
}

template <typename T> inline SimpleSparseVector<T>::SimpleSparseVector(const int& sute) {
  assert(sute == 0);
  // already initialized by compiler, initializer doesn't need this.
//  entity = map<int, T>();
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
  SimpleSparseVector<T> res(*this);
  for(auto itr(res.entity.begin()); itr != res.entity.end(); ++ itr)
    itr->second = - itr->second;
  return res;
}

template <typename T> inline SimpleSparseVector<T> SimpleSparseVector<T>::operator + (const SimpleSparseVector<T>& other) const {
  SimpleSparseVector<T> res(*this);
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
  SimpleSparseVector<T> res(*this);
  return res -= other;
}

template <typename T> inline SimpleSparseVector<T>& SimpleSparseVector<T>::operator -= (const SimpleSparseVector<T>& other) {
  return *this += - other;
}

template <typename T> template <typename U> inline SimpleSparseVector<T> SimpleSparseVector<T>::operator * (const U& other) const {
  SimpleSparseVector<T> res(*this);
  return res *= other;
}

template <typename T> template <typename U> inline SimpleSparseVector<T>& SimpleSparseVector<T>::operator *= (const U& other) {
  for(auto itr(entity.begin()); itr != entity.end(); ++ itr)
    itr->second *= other;
  return *this;
}

template <typename T> template <typename U> inline SimpleSparseVector<T> SimpleSparseVector<T>::operator / (const U& other) const {
  SimpleSparseVector<T> res(*this);
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


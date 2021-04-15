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
  inline       SimpleMatrix<T>  LSVD() const;
  inline       pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > GSVD(const SimpleMatrix<T>& src) const;
  inline       SimpleVector<T>  inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const;
  template <typename U> inline SimpleMatrix<U> real() const;
  template <typename U> inline SimpleMatrix<U> imag() const;
  template <typename U> inline SimpleMatrix<U> cast() const;
  inline const int& rows() const;
  inline const int& cols() const;
  inline       void resize(const int& rows, const int& cols);
private:
  // this isn't better idea for faster calculations.
  SimpleVector<T>* entity;
  int              erows;
  int              ecols;
  num_t            epsilon;
};

template <typename T> inline SimpleMatrix<T>::SimpleMatrix() {
  erows  = 0;
  ecols  = 0;
  entity = NULL;
#if defined(_FLOAT_BITS_)
  epsilon = T(1) >> int64_t(mybits - 1);
#else
  epsilon = std::numeric_limits<T>::epsilon();
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
  epsilon = T(1) >> int64_t(mybits - 1);
#else
  epsilon = std::numeric_limits<T>::epsilon();
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

template <typename T> inline SimpleMatrix<T> SimpleMatrix<T>::LSVD() const {
  if(this->cols() < this->rows()) {
    auto res((* this) * this->transpose().LSVD());
    for(int i = 0; i < res.cols(); i ++) {
      const auto r2(res.col(i).dot(res.col(i)));
      if(epsilon < r2)
        res.setCol(i, res.col(i) / sqrt(r2));
      else {
        res = res.transpose();
        for(int j = i; j < res.rows(); j ++)
          res.row(j) *= T(0);
        SimpleVector<T> ek(res.cols());
        for(int j = 0; j < res.rows() && i < res.rows(); j ++) {
          for(int k = 0; k < ek.size(); k ++)
            ek[k] = j == k ? T(1) : T(0);
          const auto ekp(ek - res.projectionPt(ek));
          const auto e2(ekp.dot(ekp));
          if(e2 <= epsilon) continue;
          res.row(i ++) = ekp / sqrt(e2);
        }
        res = res.transpose();
        break;
      }
    }
    return res;
  }
  auto s((*this) * this->transpose());
  SimpleMatrix<T> left(s.rows(), s.rows());
  for(int i = 0; i < left.rows(); i ++)
    for(int j = 0; j < left.cols(); j ++)
      left(i, j) = i == j ? T(1) : T(0);
  for(int i = 0; i < this->rows(); i ++)
    for(int j = 0; j < this->cols(); j ++)
      assert(isfinite((*this)(i, j)));
  for(int i = 0; i < s.rows(); i ++)
    for(int j = 0; j < s.cols(); j ++)
      assert(isfinite(s(i, j)));
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
          const T ratio(work(j, i) / eii);
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
    SimpleVector<T> ek(ev.size());
    for(int i = 0, i2 = 1; i < work.rows() && i2 < work.rows(); i ++) {
      for(int j = 0; j < ek.size(); j ++)
        ek[j] = i == j ? T(1) : T(0);
      ek -= work.projectionPt(ek);
      if(ek.dot(ek) <= epsilon)
        continue;
      work.row(i2 ++) = ek / sqrt(ek.dot(ek));
    }
    SimpleMatrix<T> mul(s.rows(), s.cols());
    for(int i = 0; i < mul.rows(); i ++)
      for(int j = 0; j < mul.cols(); j ++)
        mul(i, j) = i < ii || j < ii ? (i == j ? T(1) : T(0)) :
          work(i - ii, j - ii);
    left = left * mul.transpose();
    s    = mul * s * mul.transpose();
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
template <typename T> inline pair<pair<SimpleMatrix<T>, SimpleMatrix<T> >, SimpleMatrix<T> > SimpleMatrix<T>::GSVD(const SimpleMatrix<T>& src) const {
  // refered from : https://en.wikipedia.org/wiki/Generalized_singular_value_decomposition .
  assert(this->cols() == src.cols());
  SimpleMatrix<T> C(this->rows() + src.rows(), this->cols());
  for(int i = 0; i < this->rows(); i ++)
    C.row(i) = this->row(i);
  for(int i = 0; i < src.rows(); i ++)
    C.row(i + this->rows()) = src.row(i);
  const auto P(C.LSVD());
  SimpleVector<T> d(this->cols());
  const auto PtC(P.transpose() * C);
        auto Qt(PtC);
  for(int i = 0; i < d.size(); i ++)
    Qt.row(i) /= (d[i] = sqrt(PtC.row(i).dot(PtC.row(i))));
  SimpleMatrix<T> P1(this->rows(), d.size());
  for(int i = 0; i < P1.rows(); i ++)
    P1.row(i) = P.row(i);
  auto U1(P1.LSVD());
  auto Wt(U1.transpose() * P1);
  for(int i = 0; i < Wt.rows(); i ++)
    Wt.row(i) /= sqrt(Wt.row(i).dot(Wt.row(i)));
  SimpleMatrix<T> P2(src.rows(), d.size());
  for(int i = 0; i < P2.rows(); i ++)
    P2.row(i) = P.row(i + P1.rows());
  const auto P21W(P2 * Wt.transpose());
  const auto U2(P21W.LSVD());
  auto& WtD(Wt);
  for(int i = 0; i < WtD.cols(); i ++)
    WtD.setCol(i, WtD.col(i) * d[i]);
  assert(WtD.rows() == WtD.cols());
  SimpleMatrix<T> QtonWtDt(WtD.rows(), WtD.cols());
  for(int i = 0; i < QtonWtDt.rows(); i ++)
    for(int j = 0; j < QtonWtDt.cols(); j ++)
      QtonWtDt(i, j) = T(0);
  QtonWtDt.row(0) = WtD.row(0) / sqrt(WtD.row(0).dot(WtD.row(0)));
  for(int i = 1; i < QtonWtDt.rows(); i ++) {
    QtonWtDt.row(i)  = WtD.row(i) - QtonWtDt.projectionPt(WtD.row(i));
    QtonWtDt.row(i) /= sqrt(QtonWtDt.row(i).dot(QtonWtDt.row(i)));
  }
  return make_pair(make_pair(move(U1), move(U2)), QtonWtDt * Qt);
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

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::inner(const SimpleVector<T>& bl, const SimpleVector<T>& bu) const {
  static const auto ee(pow(epsilon, T(1) / T(8)));
  assert(this->rows() == bl.size() && this->rows() == bu.size() &&
         0 < this->cols() && 0 < this->rows());
  // bu - bb == A, bl - bb == - A <=> bu - bl == 2 A. 
  const auto bb((bu + bl) / T(2));
  const auto upper(bu - bb);
  SimpleMatrix<T> A(this->rows() * 2 - 1 + this->cols() * 2, this->cols() + 1);
  SimpleVector<T> one(A.rows());
  vector<pair<T, int> > fidx;
  fidx.reserve(this->rows());
  for(int i = 0; i < this->rows(); i ++) {
    for(int j = 0; j < this->cols(); j ++)
      A(i, j) = (*this)(i, j);
    A(i, this->cols()) = - bb[i];
    one[i]          = T(1);
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
    one[i + this->rows()] = T(1);
  }
  // this is tricky but one of them is fixed if it is needed.
  for(int i = this->rows() * 2 - 1; i < A.rows(); i ++) {
    const auto ii(i - (this->rows() * 2 - 1));
    for(int j = 0; j < this->cols(); j ++)
      A(i, j) = T(j == ii / 2 ? (ii & 1 ? - 1 : 1) : 0);
    A(i, this->cols()) = T(ii & 1 ? 1 : - 1);
    A.row(i) /= ee;
    one[i] = T(1);
  }
  // N.B. we now have |[A -bb] [x t]| <= 1 condition.
  // N.B. there's no difference |[A - bb] [x t]|^2 <= 1 condition in this.
  //      but not with mixed condition.
  SimpleMatrix<T> Pt(A.cols(), A.rows());
  for(int i = 0; i < Pt.rows(); i ++)
    for(int j = 0; j < Pt.cols(); j ++)
      Pt(i, j) = T(0);
  vector<int> residue;
  for(int i = 0; i < Pt.rows(); i ++) {
    const auto Atrowi(A.col(i));
    const auto work(Atrowi - Pt.projectionPt(Atrowi));
    const auto n2(work.dot(work));
    if(n2 <= epsilon) {
      residue.emplace_back(i);
      continue;
    }
    Pt.row(i) = work / sqrt(n2);
  }
  int ii(0);
  for(int j = 0; j < Pt.cols() && ii < residue.size(); j ++) {
    SimpleVector<T> ek(Pt.cols());
    for(int k = 0; k < Pt.cols(); k ++)
      ek[k] = T(j == k ? 1 : 0);
    ek -= Pt.projectionPt(ek);
    const auto n2(ek.dot(ek));
    if(n2 <= epsilon) continue;
    Pt.row(residue[ii ++]) = ek / sqrt(n2);
  }
  assert(residue.size() <= ii);
  const auto R(Pt * A);
  // we now have: Q [R [x t] ] <= {0, 1}^m cond.
  const auto on(Pt.projectionPt(one));
  fidx.reserve(fidx.size() + on.size());
  for(int i = 0; i < on.size(); i ++)
    if(isfinite(on[i]) && ! isnan(on[i]))
      fidx.emplace_back(make_pair(abs(on[i]), i));
  if(fidx.size())
    sort(fidx.begin(), fidx.end());
  vector<int> fix;
  fix.reserve(this->cols());
  // sort by: |<Q^t(1), q_k>|, we subject to minimize each, to do this,
  //   maximize minimum q_k orthogonality.
  for(int idx = 0; fix.size() < this->cols() && idx < fidx.size(); idx ++) {
    const auto& iidx(fidx[idx].second);
    const auto  orth(Pt.col(iidx));
    const auto  n2(orth.dot(orth));
    if(n2 <= epsilon)
      continue;
    fix.emplace_back(fidx[idx].second);
    // N.B. O(mn) can be writed into O(lg m + lg n) in many core cond.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / n2);
  }
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

#define _SIMPLELIN_
#endif


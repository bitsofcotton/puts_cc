/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/* BSD 3-Clause License:
 * Copyright (c) 2013 - 2018, kazunobu watatsu.
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
private:
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
  inline       SimpleVector<T>  solve(SimpleVector<T> other) const;
  inline       SimpleVector<T>  projectionPt(const SimpleVector<T>& other) const;
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
};

template <typename T> inline SimpleMatrix<T>::SimpleMatrix() {
  erows  = 0;
  ecols  = 0;
  entity = NULL;
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
  return;
}

template <typename T> inline SimpleMatrix<T>::SimpleMatrix(SimpleMatrix<T>&& other) {
  erows  = move(other.erows);
  ecols  = move(other.ecols);
  entity = move(other.entity);
  other.entity = NULL;
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

template <typename T> inline SimpleVector<T> SimpleMatrix<T>::solve(SimpleVector<T> other) const {
  assert(0 <= erows && 0 <= ecols && erows == ecols && entity && erows == other.size());
  SimpleMatrix<T> work(*this);
  for(int i = 0; i < erows; i ++) {
    int xchg = i;
    for(int j = i + 1; j < erows; j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    SimpleVector<T> buf(work.entity[i]);
    T               buf2(other[i]);
    work.entity[i]    = work.entity[xchg];
    other[i]          = other[xchg];
    work.entity[xchg] = buf;
    other[xchg]       = buf2;
    const SimpleVector<T>& ei(work.entity[i]);
    const T&               eii(ei[i]);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = i + 1; j < erows; j ++) {
      const T ratio(work.entity[j][i] / eii);
      work.entity[j] -= ei       * ratio;
      other[j]       -= other[i] * ratio;
    }
  }
  for(int i = erows - 1; 0 <= i; i --) {
    const T buf(other[i] / work.entity[i][i]);
    if(!isfinite(buf) || isnan(buf)) {
      // throw "Non full rank";
      // assert(!isfinite(work.entity[i][i] / other[i]) || isnan(work.entity[i][i] / other[i]));
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

#define _SIMPLELIN_
#endif

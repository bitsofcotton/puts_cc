/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/* BSD 3-Clause License:
 * Copyright (c) 2018, bitsofcotton.
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

#if !defined(_SPARSE_TENSOR_)

#include <cstdio>
#include <vector>
#include <utility>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;
using std::sort;
using std::pair;
using std::make_pair;
using std::distance;
using std::upper_bound;

template <typename T> bool orderPairVec(const T& a, const T& b) {
  return a.first < b.first;
}

template <typename T> class SimpleSparseVector {
public:
  SimpleSparseVector();
  SimpleSparseVector(const int& sute);
  SimpleSparseVector(const SimpleSparseVector<T>& other);
  SimpleSparseVector(SimpleSparseVector<T>&& other);
  ~SimpleSparseVector();
  
        SimpleSparseVector<T>  operator -  () const;
        SimpleSparseVector<T>  operator +  (const SimpleSparseVector<T>& other) const;
  const SimpleSparseVector<T>& operator += (const SimpleSparseVector<T>& other);
        SimpleSparseVector<T>  operator -  (const SimpleSparseVector<T>& other) const;
  const SimpleSparseVector<T>& operator -= (const SimpleSparseVector<T>& other);
  template <typename U>       SimpleSparseVector<T>  operator *  (const U& other) const;
  template <typename U> const SimpleSparseVector<T>& operator *= (const U& other);
  template <typename U>       SimpleSparseVector<T>  operator /  (const U& other) const;
  template <typename U> const SimpleSparseVector<T>& operator /= (const U& other);
        SimpleSparseVector<T>& operator =  (const SimpleSparseVector<T>& other);
        bool                   operator != (const SimpleSparseVector<T>& other) const;
        bool                   operator == (const SimpleSparseVector<T>& other) const;
        T  dot         (const SimpleSparseVector<T>& other) const;
        T& operator [] (const int& idx);
  const T  operator [] (const int& idx) const;
private:
  vector<pair<int, T> > entity;
};

template <typename T> SimpleSparseVector<T>::SimpleSparseVector() {
  entity = vector<pair<int, T> >();
  return;
}

template <typename T> SimpleSparseVector<T>::SimpleSparseVector(const int& sute) {
  assert(sute == 0);
  entity = vector<pair<int, T> >();
  return;
}

template <typename T> SimpleSparseVector<T>::SimpleSparseVector(const SimpleSparseVector<T>& other) {
  entity = vector<pair<int, T> >(other.entity);
  return;
}

template <typename T> SimpleSparseVector<T>::SimpleSparseVector(SimpleSparseVector<T>&& other) {
  entity = other.entity;
  return;
}

template <typename T> SimpleSparseVector<T>::~SimpleSparseVector() {
  return;
}

template <typename T> SimpleSparseVector<T> SimpleSparseVector<T>::operator - () const {
  SimpleSparseVector<T> res(*this);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < res.entity.size(); i ++)
    res.entity[i].second = - res.entity[i].second;
  return res;
}

template <typename T> SimpleSparseVector<T> SimpleSparseVector<T>::operator + (const SimpleSparseVector<T>& other) const {
  SimpleSparseVector<T> res(*this);
  return res += other;
}

template <typename T> const SimpleSparseVector<T>& SimpleSparseVector<T>::operator += (const SimpleSparseVector<T>& other) {
  for(int i = 0; i < other.entity.size(); i ++) {
    const int j(distance(entity.begin(), upper_bound(entity.begin(), entity.end(), other.entity[i], orderPairVec<pair<int, T> >)) - 1);
    if(0 <= j && j < entity.size() && entity[j].first == other.entity[i].first)
      entity[j].second += other.entity[i].second;
    else if(other.entity[i].second != T(0))
      (*this)[other.entity[i].first] = other.entity[i].second;
  }
  return *this;
}

template <typename T> SimpleSparseVector<T> SimpleSparseVector<T>::operator - (const SimpleSparseVector<T>& other) const {
  SimpleSparseVector<T> res(*this);
  return res -= other;
}

template <typename T> const SimpleSparseVector<T>& SimpleSparseVector<T>::operator -= (const SimpleSparseVector<T>& other) {
  return *this += - other;
}

template <typename T> template <typename U> SimpleSparseVector<T> SimpleSparseVector<T>::operator * (const U& other) const {
  SimpleSparseVector<T> res(*this);
  return res *= other;
}

template <typename T> template <typename U> const SimpleSparseVector<T>& SimpleSparseVector<T>::operator *= (const U& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < entity.size(); i ++)
    entity[i].second *= other;
  return *this;
}

template <typename T> template <typename U> SimpleSparseVector<T> SimpleSparseVector<T>::operator / (const U& other) const {
  SimpleSparseVector<T> res(*this);
  return res /= other;
}

template <typename T> SimpleSparseVector<T>& SimpleSparseVector<T>::operator = (const SimpleSparseVector<T>& other) {
  entity = other.entity;
  return *this;
}

template <typename T> bool SimpleSparseVector<T>::operator != (const SimpleSparseVector<T>& other) const {
  // XXX : this is imcomplete.
  return entity != other.entity;
}

template <typename T> bool SimpleSparseVector<T>::operator == (const SimpleSparseVector<T>& other) const {
  return ! (*this != other);
}

template <typename T> template <typename U> const SimpleSparseVector<T>& SimpleSparseVector<T>::operator /= (const U& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < entity.size(); i ++)
    entity[i].second /= other;
  return *this;
}

template <typename T> T SimpleSparseVector<T>::dot(const SimpleSparseVector<T>& other) const {
  T res(0);
  SimpleSparseVector<T> work;
  for(int i = 0; i < other.entity.size(); i ++) {
    const int j(distance(entity.begin(), upper_bound(entity.begin(), entity.end(), other.entity[i], orderPairVec<pair<int, T> >)) - 1);
    if(0 <= j && j < entity.size() && entity[j].first == other.entity[i].first)
      work[entity[j].first] = entity[j].second * other.entity[i].second;
  }
  for(int i = 0; i < work.entity.size(); i ++)
    res += work.entity[i].second;
  return res;
}

template <typename T> T& SimpleSparseVector<T>::operator [] (const int& idx) {
  assert(0 <= idx);
  const auto lidx(make_pair(int(idx), T(0)));
  const int  ii(distance(entity.begin(), upper_bound(entity.begin(), entity.end(), lidx, orderPairVec<pair<int, T> >)) - 1);
  if(0 <= ii && ii < entity.size() && entity[ii].first == idx)
    return entity[ii].second;
  else {
    entity.push_back(make_pair(int(idx), T(0)));
    sort(entity.begin(), entity.end(), orderPairVec<pair<int, T> >);
  }
  const int  i(distance(entity.begin(), upper_bound(entity.begin(), entity.end(), lidx, orderPairVec<pair<int, T> >)) - 1);
  return entity[i].second;
}

template <typename T> const T SimpleSparseVector<T>::operator [] (const int& idx) const {
  assert(0 <= idx);
  const auto lidx(make_pair(int(idx), T(0)));
  const int  i(distance(entity.begin(), upper_bound(entity.begin(), entity.end(), lidx, orderPairVec<pair<int, T> >)) - 1);
  if(0 <= i && i < entity.size() && entity[i].first == idx)
    return entity[i].second;
  return T(0);
}

template <typename T> using SimpleSparseMatrix = SimpleSparseVector<SimpleSparseVector<T> >;
template <typename T> using SimpleSparseTensor = SimpleSparseVector<SimpleSparseVector<SimpleSparseVector<T> > >;

#define _SPARSE_TENSOR_
#endif


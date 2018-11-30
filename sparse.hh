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

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::map;
using std::move;

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

#define _SPARSE_TENSOR_
#endif


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
#if !defined(_CORPUS_)

#include <Eigen/Core>
#include <Eigen/SVD>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iterator>
#include <iostream>

using std::cerr;
using std::endl;
using std::flush;
using std::string;
using std::vector;
using std::sort;
using std::distance;
using std::equal_range;
using std::upper_bound;
using std::unique;
using std::find;
using std::pair;
using std::make_pair;
using std::to_string;
using std::isfinite;
using std::sqrt;
using std::min;
using std::max;
using std::abs;
using std::exp;
using std::log;

template <typename T> bool equalStrClip(const T& a, const T& b) {
  int cmp(0), jidx(0);
  for( ; !cmp && jidx < min(a.size(), b.size()); jidx ++)
    cmp = a[jidx] - b[jidx];
  return !cmp && min(a.size(), b.size()) <= jidx;
}

template <typename T> bool lessEqualStrClip(const T& a, const T& b) {
  return a < b ||  equalStrClip<T>(a, b);
}

template <typename T> bool lessNotEqualStrClip(const T& a, const T& b) {
  return a < b && !equalStrClip<T>(a, b);
}

template <typename T, typename U> class corpus {
public:
  typedef Eigen::Matrix<T,   Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T,   Eigen::Dynamic, 1>              Vec;
  typedef Eigen::Matrix<Vec, Eigen::Dynamic, Eigen::Dynamic> Tensor;
  
  corpus();
  ~corpus();
  
  void init(const vector<U>& words0, const int& nthresh, const int& Nthresh, const int& Mwords = 150);
  corpus<T, U>&    operator = (const corpus<T, U>& other);
  const void       compute(const U& input, const vector<U>& delimiter = vector<U>());
  const U          getAttributed(const vector<U>& highlight) const;
  const vector<U>& getWords() const;
  const Tensor&    getCorpus() const;
private:
  U                    orig;
  vector<U>            words0;
  vector<U>            words;
  vector<vector<int> > ptrs0;
  vector<vector<int> > ptrs;
  vector<int>          pdelim;
  Tensor               corpust;
  int                  nthresh;
  int                  Nthresh;
  int                  Mwords;
  T                    Midx;
   
  void getWordPtrs(const U& input, const vector<U>& delimiter);
  void corpusEach();
};

template <typename T, typename U> corpus<T,U>::corpus() {
  init(vector<U>(), 1, 1, 1);
}

template <typename T, typename U> corpus<T,U>::~corpus() {
  // auto called destructors for string.
  ;
}

template <typename T, typename U> void corpus<T,U>::init(const vector<U>& words0, const int& nthresh, const int& Nthresh, const int& Mwords) {
  this->orig    = U();
  this->words0  = words0;
  sort(this->words0.begin(), this->words0.end());
  this->words0.erase(unique(this->words0.begin(), this->words0.end()), this->words0.end());
  this->words   = vector<U>();
  this->ptrs0   = vector<vector<int> >();
  for(int i = 0; i < words0.size(); i ++)
    ptrs0.push_back(vector<int>());
  this->ptrs    = vector<vector<int> >();
  this->pdelim  = vector<int>();
  this->corpust = Tensor();
  this->nthresh = nthresh;
  this->Nthresh = Nthresh;
  this->Midx    = 1;
  this->Mwords  = Mwords;
  return;
}

template <typename T, typename U> corpus<T, U>& corpus<T, U>::operator = (const corpus<T, U>& other) {
  words0  = other.words0;
  words   = other.words;
  ptrs0   = other.ptrs0;
  ptrs    = other.ptrs;
  corpust = other.corpust;
  nthresh = other.nthresh;
  Nthresh = other.Nthresh;
  return *this;
}

template <typename T, typename U> const void corpus<T,U>::compute(const U& input, const vector<U>& delimiter) {
  cerr << "Getting word pointers " << flush;
  getWordPtrs(input, delimiter);
  if(Mwords < words.size()) {
    cerr << "exceeds Mwords." << endl;
    return;
  }
  cerr << "Corpus" << flush;
  corpusEach();
  return;
}

template <typename T, typename U> const U corpus<T,U>::getAttributed(const vector<U>& highlight) const {
  vector<int> idxs;
  for(int i = 0; i < highlight.size(); i ++)
    if(highlight[i] != U("^") && highlight[i] != U("$"))
      for(int j = 0; j < words.size(); j ++)
        if(highlight[i] == words[j]) {
          idxs.push_back(j);
          break;
        }
  sort(idxs.begin(), idxs.end());
  idxs.erase(unique(idxs.begin(), idxs.end()), idxs.end());
  U result;
  int    i;
  for(i = 0; i < orig.size(); ) {
    int idx(- 1);
    for(int j = 0; j < idxs.size(); j ++)
      if(binary_search(ptrs[idxs[j]].begin(), ptrs[idxs[j]].end(), i - 1)) {
        idx = idxs[j];
        break;
      }
    if(0 <= idx) {
      result += U("<font class=\"match\">");
      result += words[idx];
      result += U("</font>");
      i      += words[idx].size();
    } else
      result += orig[i ++];
  }
  return result;
}

template <typename T, typename U> const vector<U>& corpus<T,U>::getWords() const {
  return words;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic>& corpus<T,U>::getCorpus() const {
  return corpust;
}

template <typename T, typename U> void corpus<T,U>::getWordPtrs(const U& input, const vector<U>& delimiter) {
  sort(words0.begin(), words0.end());
  pdelim = vector<int>();
  pdelim.push_back(0);
  U work;
  vector<int> matchwidx;
  vector<int> matchidxs;
  int dM(0);
  for(int i = 0; i < delimiter.size(); i ++)
    dM = max(dM, int(delimiter[i].size()));
  vector<U> workd;
  for(int i = 0; i < dM; i ++) {
    workd.push_back(U(""));
    for(int j = i; j < dM; j ++)
      workd[i] += U(" ");
  }
  orig = U("");
  for(int i = 0; i < input.size(); i ++)
    orig += input[i];
  int  i(0), i0(0);
  for( ; i < input.size(); i ++) {
    work += input[i];
    for(int ii = 0; ii < workd.size(); ii ++) {
      workd[ii]  = workd[ii].substr(1, workd[ii].size() - 1);
      workd[ii] += input[i];
      for(int j = 0; j < delimiter.size(); j ++)
        if(workd[ii] == delimiter[j] && pdelim[pdelim.size() - 1] < i)
          pdelim.push_back(i);
    }
    auto lo(words0.begin() + distance(words0.begin(), upper_bound(words0.begin(), words0.end(), work, lessEqualStrClip<U>)));
    auto up(words0.begin() + distance(words0.begin(), upper_bound(words0.begin(), words0.end(), work, lessNotEqualStrClip<U>)));
    bool match(false);
    for(auto itr(lo); itr < up; ++ itr)
      if(equalStrClip<U>(work, *itr)) {
        if(work.size() == itr->size()) {
          matchwidx.push_back(distance(words0.begin(), itr));
          matchidxs.push_back(i0);
        } else if(work.size() < itr->size())
          match = true;
      }
    if(match && input[i + 1])
      continue;
    if(matchwidx.size() > 0) {
      int j = matchwidx.size() - 1;
      ptrs0[matchwidx[j]].push_back(matchidxs[j]);
      Midx = matchidxs[j];
      matchwidx = vector<int>();
      matchidxs = vector<int>();
    }
    i     -= work.size() - 1;
    i0     = i;
    work   = U();
  }
  {
    if(matchwidx.size() > 0) {
      int j = matchwidx.size() - 1;
      ptrs0[matchwidx[j]].push_back(matchidxs[j]);
      Midx = matchidxs[j];
      matchwidx = vector<int>();
      matchidxs = vector<int>();
    }
  }
  words = vector<U>();
  ptrs  = vector<vector<int> >();
  vector<int> head, tail;
  words.push_back(U("^"));
  head.push_back(0);
  ptrs.push_back(head);
  for(auto itr = words0.begin(); itr != words0.end(); ++ itr) {
    const int idx = distance(words0.begin(), itr);
    if(itr->size() && ptrs0[idx].size()) {
      words.push_back(*itr);
      ptrs.push_back(ptrs0[idx]);
    }
  }
  words.push_back(U("$"));
  tail.push_back(Midx + 1);
  ptrs.push_back(tail);
  pdelim.push_back(Midx + 2);
  cerr << words.size() - 2 << " words used, ptr: " << i << endl;
  return;
}

template <typename T, typename U> void corpus<T,U>::corpusEach() {
  corpust = Tensor(words.size(), words.size());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < words.size(); i ++) {
    cerr << "." << flush;
    for(int j = 0; j < words.size(); j ++) {
      int kk(0);
      corpust(i, j) = Vec(words.size());
      for(int k = 0; k < corpust(i, j).size(); k ++)
        corpust(i, j)[k] = 0;
      if(!ptrs[i].size() || !ptrs[j].size())
        continue;
      for(int k = 0; k < words.size(); k ++) {
        if(!ptrs[k].size())
          continue;
        int ctru = 0;
        int ctrv = 0;
        for(auto itr = ptrs[k].begin(); itr != ptrs[k].end(); ++ itr) {
          while(ctru < ptrs[i].size() && ptrs[i][ctru] < *itr) ctru ++;
          ctru --;
          if(ctru < 0 || ptrs[i].size() <= ctru) {
            if(ctru < 0)
              ctru = 0;
            continue;
          }
          while(ctrv < ptrs[j].size() && ptrs[j][ctrv] < *itr) ctrv ++;
          if(ptrs[j].size() <= ctrv || ptrs[j][ctrv] < *itr)
            break;
          for( ; kk < pdelim.size() - 1; kk ++)
            if(pdelim[kk] <= *itr && *itr <= pdelim[kk + 1])
              break;
          if(! (pdelim[kk] <= ptrs[i][ctru] &&
                              ptrs[i][ctru] <= pdelim[kk + 1] &&
                pdelim[kk] <= ptrs[j][ctrv] &&
                              ptrs[j][ctrv] <= pdelim[kk + 1]) )
            continue;
          // XXX configure me:
          const T buf0(log(T(abs(*itr + .5 - ptrs[i][ctru])) * T(2) * exp(T(1))));
          const T buf1(log(T(abs(*itr + .5 - ptrs[j][ctrv])) * T(2) * exp(T(1))));
          // const T buf0(abs(*itr + .5 - ptrs[i][ctru]));
          // const T buf1(abs(*itr + .5 - ptrs[j][ctrv]));
          const T work(T(1) / (buf0 * buf0 + buf1 * buf1));
          if(isfinite(work))
            corpust(i, j)[k] += sqrt(work) / Midx;
        }
      }
    }
  }
  return;
}


template <typename T, typename U> class corpushl {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>   Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                Vec;
  typedef Eigen::Matrix<Vec, Eigen::Dynamic, Eigen::Dynamic> Tensor;
  
  corpushl();
  corpushl(const corpus<T, U>&   obj);
  corpushl(const corpushl<T, U>& obj);
  ~corpushl();
  
  const corpushl<T, U>  cast(const vector<U>& words) const;
  const corpushl<T, U>& operator += (const corpushl<T, U>& other);
  const corpushl<T, U>& operator -= (const corpushl<T, U>& other);
  const corpushl<T, U>& operator *= (const T& t);
  const bool            operator <  (const corpushl<T, U>& other) const;
  const corpushl<T, U>  operator +  (const corpushl<T, U>& other) const;
  const corpushl<T, U>  operator -  () const;
  const corpushl<T, U>  operator -  (const corpushl<T, U>& other) const;
  const corpushl<T, U>  operator *  (const T& t)                  const;
  const corpushl<T, U>& operator =  (const corpushl<T, U>& other);
  const bool            operator == (const corpushl<T, U>& other) const;
  const bool            operator != (const corpushl<T, U>& other) const;
  const corpushl<T, U>  withDetail(const U& word, const corpushl<T, U>& other);
  const T               cdot(const corpushl<T, U>& other) const;
  const corpushl<T, U>  match2relPseudo(const corpushl<T, U>& other) const;
  const T               prej(const corpushl<T, U>& prejs) const;
  const T               prej2(const vector<corpushl<T, U> >& prej0, const vector<corpushl<T, U> >& prej1, const T& thresh) const;
  const corpushl<T, U>  lackOfInsist() const;
  const corpushl<T, U>  invertInsist() const;
  const T               culturalConflicts(const corpushl<T, U>& base) const;
  const corpushl<T, U>  conflictPart();
  const vector<U>&      getWords() const;
  const Tensor&         getCorpus() const;
  const U               serialize() const;
  const corpushl<T, U>  abbrev(const U& word, const corpushl<T, U>& mean);
  const vector<U>       reverseLink(const corpushl<T, U>& orig) const;
  const U               reverseLink(const corpus<T, U>& orig) const;
  const pair<T, T>      compareStructure(const corpushl<T, U>& src, const T& thresh = T(1e-4), const T& thresh2 = T(.125)) const;
  const corpushl<T, U>  reDig(const T& ratio);
  const corpushl<T, U>  simpleThresh(const T& ratio) const;

private:
  const U               serializeSub(const vector<int>& idxs) const;
  const Vec             singularValues() const;
  vector<U>             words;
  Tensor                corpust;
   
  vector<U>      gatherWords(const vector<U>& in0, const vector<U>& in1, vector<int>& ridx0, vector<int>& ridx1) const;
  const Tensor   prepareDetail(const corpushl<T, U>& other, const vector<U>& workwords, const int& eidx, const vector<int>& ridx0, const vector<int>& ridx1, const vector<int>& ridx2);
};

template <typename T, typename U> corpushl<T,U>::corpushl() {
  words   = vector<U>();
  corpust = Tensor();
}

template <typename T, typename U> corpushl<T,U>::~corpushl() {
  // auto called destructors for string.
  ;
}

template <typename T, typename U> corpushl<T,U>::corpushl(const corpus<T, U>& obj) {
  words   = vector<U>(obj.getWords());
  corpust = Tensor(obj.getCorpus());
}

template <typename T, typename U> corpushl<T,U>::corpushl(const corpushl<T, U>& obj) {
  *this = obj;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::cast(const vector<U>& words) const {
  cerr << "cast" << endl;
  corpushl<T, U> result;
  vector<U>      sword(words);
  vector<int>    idxs;
  sort(sword.begin(), sword.end());
  for(int i = 0; i < sword.size(); i ++) {
    auto p(find(words.begin(), words.end(), sword[i]));
    if(*p == sword[i]) {
      idxs.push_back(distance(words.begin(), p));
      result.words.push_back(*p);
    }
  }
  result.corpust = Tensor(idxs.size(), idxs.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < idxs.size(); i ++)
    for(int j = 0; j < idxs.size(); j ++) {
      result.corpust(i, j) = Vec(idxs.size());
      for(int k = 0; k < idxs.size(); k ++)
        result.corpust(i, j)[k] = corpust(idxs[i], idxs[j])[idxs[k]];
    }
  return result;
}

template <typename T, typename U> const corpushl<T, U>& corpushl<T, U>::operator = (const corpushl<T, U>& other) {
  words   = vector<U>(other.words);
  corpust = Tensor(other.corpust);
  return *this;
}

template <typename T, typename U> const bool corpushl<T, U>::operator == (const corpushl<T, U>& other) const {
  return ! (*this != other);
}

template <typename T, typename U> const bool corpushl<T, U>::operator != (const corpushl<T, U>& other) const {
  return words != other.words || corpust != other.corpust;
}

template <typename T, typename U> const corpushl<T, U>& corpushl<T, U>::operator += (const corpushl<T, U>& other) {
  return *this = *this + other;
}

template <typename T, typename U> const corpushl<T, U>& corpushl<T, U>::operator -= (const corpushl<T, U>& other) {
  return *this = *this - other;
}

template <typename T, typename U> const corpushl<T, U>& corpushl<T, U>::operator *= (const T& t) {
  cerr << "*=" << flush;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      corpust(i, j) *= t;
  return *this;
}

template <typename T, typename U> const bool corpushl<T, U>::operator < (const corpushl<T, U>& other) const {
  return corpust.rows() < other.corpust.rows();
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator + (const corpushl<T, U>& other) const {
  corpushl<T, U> result;
  vector<int>    ridx0, ridx1;
  result.words   = gatherWords(words, other.words, ridx0, ridx1);
  result.corpust = Tensor(result.words.size(), result.words.size());
  cerr << "+" << flush;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.corpust.rows(); i ++)
    for(int j = 0; j < result.corpust.cols(); j ++) {
      result.corpust(i, j) = Vec(result.words.size());
      for(int k = 0; k < result.corpust(i, j).size(); k ++)
        result.corpust(i, j)[k] = T(0);
      if(0 <= ridx0[i] && 0 <= ridx0[j])
        for(int k = 0; k < result.corpust(i, j).size(); k ++)
          if(0 <= ridx0[k])
            result.corpust(i, j)[k] += corpust(ridx0[i], ridx0[j])[ridx0[k]];
      if(0 <= ridx1[i] && 0 <= ridx1[j])
        for(int k = 0; k < result.corpust(i, j).size(); k ++)
          if(0 <= ridx1[k])
            result.corpust(i, j)[k] += other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
    }
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator - () const {
  cerr << "-" << flush;
  corpushl<T, U> result(*this);
  result.corpust = - result.corpust;
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator - (const corpushl<T, U>& other) const {
  return (*this) + (- other);
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator * (const T& t) const {
  corpushl<T, U> work(*this);
  return work *= t;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::withDetail(const U& word, const corpushl<T, U>& other) {
  if(words.size() <= 0 || other.words.size() <= 0)
    return *this;
  auto itr(find(words.begin(), words.end(), word));
  int  fidx(distance(words.begin(), itr));
  if(!(0 <= fidx && fidx < words.size() && *itr == word))
    return *this;
  cerr << "withDetail : " << *itr << ", " << word << ": " << itr->size() << " / " << word.size() << endl;
  corpushl<T, U> result;
  vector<int> ridx0, ridx1, ridx2;
  vector<U>   workwords(gatherWords(words, other.words, ridx0, ridx1));
  int  eidx(- 1);
  bool flag1(false);
  for(int i = 0, ii = 0; i < workwords.size(); i ++) {
    if(workwords[i] == word) {
      ridx2.push_back(- 1);
      eidx = i;
      continue;
    }
    if(ridx1[i] >= 0)
      flag1 = true;
    ridx2.push_back(ii ++);
  }
  if(!flag1 || eidx < 0 || ridx0[eidx] < 0)
    return result;
  result.corpust = Tensor(workwords.size() - 1, workwords.size() - 1);
  for(int i = 0; i < workwords.size(); i ++)
    if(ridx0[i] >= 0 && ridx2[i] >= 0)
      result.words.push_back(words[ridx0[i]]);
    else if(ridx1[i] >= 0 && ridx2[i] >= 0)
      result.words.push_back(other.words[ridx1[i]]);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.corpust.rows(); i ++)
    for(int j = 0; j < result.corpust.cols(); j ++) {
      result.corpust(i, j) = Vec(workwords.size() - 1);
      for(int k = 0; k < result.corpust(i, j).size(); k ++)
        result.corpust(i, j)[k] = T(0);
    }
#if defined(_OPENMP)
#pragma omp paralell for schedule(static, 1)
#endif
  for(int i = 0; i < workwords.size(); i ++) if(0 <= ridx2[i])
    for(int j = 0; j < workwords.size(); j ++) if(0 <= ridx2[j])
      for(int k = 0; k < result.corpust(ridx2[i], ridx2[j]).size(); k ++) if(ridx2[k] >= 0)
        if(ridx0[i] >= 0 && ridx0[j] >= 0 && ridx0[k] >= 0)
          result.corpust(ridx2[i], ridx2[j])[ridx2[k]] = corpust(ridx0[i], ridx0[j])[ridx0[k]];
  result.corpust += prepareDetail(other, workwords, eidx, ridx0, ridx1, ridx2);
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::match2relPseudo(const corpushl<T, U>& other) const {
  cerr << "m2r" << flush;
  corpushl<T, U> result(*this);
  vector<int>    ridx0, ridx1;
  vector<U>      words(gatherWords(result.words, other.words, ridx0, ridx1));
  Tensor mul(result.corpust);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < mul.rows(); i ++)
    for(int j = 0; j < mul.cols(); j ++)
      for(int k = 0; k < mul(i, j).size(); k ++)
        mul(i, j)[k] = T(0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < words.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0)
    for(int j = 0; j < words.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < words.size(); k ++) if(ridx0[k] >= 0) {
        if(ridx0[k] < mul(ridx0[i], ridx0[j]).size())
          mul(ridx0[i], ridx0[j])[ridx0[k]] = 1.;
        if(ridx0[i] < mul(ridx0[j], ridx0[k]).size())
          mul(ridx0[j], ridx0[k])[ridx0[i]] = 1.;
        if(ridx0[j] < mul(ridx0[k], ridx0[i]).size())
          mul(ridx0[k], ridx0[i])[ridx0[j]] = 1.;
      }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.corpust.rows(); i ++) if(0 <= ridx0[i])
    for(int j = 0; j < result.corpust.cols(); j ++) if(0 <= ridx0[j])
      for(int k = 0; k < result.corpust(i, j).size(); k ++) if(0 <= ridx0[k] &&
          ridx0[k] < result.corpust(ridx0[i], ridx0[j]).size())
        result.corpust(ridx0[i], ridx0[j])[ridx0[k]] *=
          mul(ridx0[i], ridx0[j])[ridx0[k]];
  return result;
}

template <typename T, typename U> const T corpushl<T, U>::cdot(const corpushl<T, U>& other) const {
  T res(0);
  vector<int> ridx0, ridx1;
  vector<U>   drop(gatherWords(words, other.words, ridx0, ridx1));
  for(int i = 0; i < drop.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0)
    for(int j = 0; j < drop.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < drop.size(); k ++) if(ridx0[k] >= 0 && ridx1[k] >= 0 &&
          ridx0[k] < corpust(ridx0[i], ridx0[j]).size() &&
          ridx1[k] < other.corpust(ridx1[i], ridx1[j]).size())
        res += corpust(ridx0[i], ridx0[j])[ridx0[k]] * other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
  return res;
}

template <typename T, typename U> const T corpushl<T, U>::prej(const corpushl<T, U>& prejs) const {
  // XXX this is broken method. DO NOT USE THIS ONLY.: fixme...
  cerr << "XXX : confirm me corpushl::prej function." << endl;
  // XXX confirm me: need some other counting methods?
  // XXX fixme: amount of the word that is not said in the context is
  //            also important.
  const auto n2this(cdot(*this));
  if(n2this == T(0))
    return - T(4);
  const auto n2p(prejs.cdot(prejs));
  if(n2p == T(0))
    return T(0);
  const auto ex0(*this - match2relPseudo(prejs));
  const auto n2e0(ex0.cdot(ex0));
  if(n2e0 == T(0))
    return   T(4);
  const auto ex1(prejs - prejs.match2relPseudo(*this));
  const auto n2e1(ex1.cdot(ex1));
  if(n2e1 == T(0))
    return - T(4);
  // XXX is sign correct?
  return cdot(prejs) / sqrt(n2this * n2p) - ex0.cdot(ex1) / sqrt(n2e0 * n2e1);
}

template <typename T, typename U> const T corpushl<T, U>::prej2(const vector<corpushl<T, U> >& prej0, const vector<corpushl<T, U> >& prej1, const T& thresh) const {
  cerr << "XXX confirm me: corpushl::prej2" << endl;
  // XXX confirm me: is this correct counting method?
  // XXX fixme: also with prej func.
  corpushl<T, U> p0(*this), p1(*this);
  for(int i = 0; i < prej0.size(); i ++)
    p0 = p0.abbrev(string("P") + to_string(i), prej0[i]);
  for(int i = 0; i < prej1.size(); i ++)
    p1 = p1.abbrev(string("Q") + to_string(i), prej1[i]);
  p0 = p0.simpleThresh(thresh);
  p1 = p1.simpleThresh(thresh);
  return T(p0.words.size() - prej0.size()) / T(p1.words.size() - prej1.size());
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::lackOfInsist() const {
  cerr << "XXX confirm me: corpushl::lackOfInsist" << endl;
  assert(0);
  // XXX confirm me: this method can only count what's on the table.
  //                 so the things that isn't on the table will not be counted.
  corpushl<T, U> result;
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::invertInsist() const {
  cerr << "XXX confirm me: corpushl::invertInsist" << endl;
  assert(0);
  // XXX confirm me: this method cannot calculate in logically correct
  //                 because of it's method.
  cerr << "STUB INVERT INSIST." << endl;
  corpushl<T, U> result;
  return result;
}

template <typename T, typename U> const T corpushl<T, U>::culturalConflicts(const corpushl<T, U>& base) const {
  cerr << "XXX confirm me: corpushl::culturalConflicts" << endl;
  assert(0);
  return T(0);
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::conflictPart() {
  cerr << "XXX confirm me: corpushl::conflictPart" << endl;
  assert(0);
  // search conflict parts.
  // dictionary base of the word 'NOT' is needed.
  corpushl<T, U> result;
  return result;
}

template <typename T, typename U> const U corpushl<T, U>::serialize() const {
  cerr << " serialize" << flush;
  if(corpust.rows() < 1 || corpust.cols() < 1)
    return U("*NULL*");
  auto plus(*this), minus(*this);
  for(int i = 0; i < plus.corpust.rows(); i ++)
    for(int j = 0; j < plus.corpust.cols(); j ++)
      for(int k = 0; k < plus.corpust(i, j).size(); k ++)
        if(plus.corpust(i, j)[k] < T(0)) {
          plus.corpust(i, j)[k]  = T(0);
          minus.corpust(i, j)[k] = - minus.corpust(i, j)[k];
        } else
          minus.corpust(i, j)[k] = T(0);
  vector<int> entire;
  for(int i = 0; i < corpust.rows(); i ++)
    entire.push_back(i);
  return plus.serializeSub(entire)  + U(".<br/>\n-") +
         minus.serializeSub(entire) + U(".");
}

template <typename T, typename U> const U corpushl<T, U>::serializeSub(const vector<int>& idxs) const {
  cerr << "." << flush;
  if(idxs.size() <= 1) {
    if(idxs.size())
      return words[idxs[0]];
    return U();
  }
  vector<pair<int, int> > cscore;
  // N.B. i0 - i1 - i2 is stored in corpust(i0, i2)[i1].
  for(int i = 0; i < idxs.size(); i ++) {
    int lscore(0);
    for(int j = 0; j < idxs.size(); j ++)
      if(0 <= idxs[j] && idxs[j] < corpust.rows()) {
        for(int k = 0; k < idxs.size(); k ++)
          if(j != k &&
             0 <= idxs[k] && idxs[k] < corpust.cols() &&
             0 <= idxs[i] && idxs[i] < corpust(idxs[j], idxs[k]).size() &&
             corpust(idxs[j], idxs[k])[idxs[i]] != T(0))
            lscore --;
      }
    cscore.push_back(make_pair(lscore, idxs[i]));
  }
  sort(cscore.begin(), cscore.end());
  for(int si = 0; si < cscore.size(); si ++) {
    vector<int> middle, left, right;
    middle.push_back(cscore[si].second);
    if(!cscore[si].first) {
      string result;
      for(int i = 0; i < cscore.size(); i ++)
        if(cscore[i].first != T(0))
          result += words[cscore[i].second];
      return result;
    }
    vector<pair<int, int> > score;
    for(int i = 0; i < idxs.size(); i ++)
      if(idxs[i] != middle[0] && 0 <= idxs[i] && idxs[i] < corpust.rows()) {
        int lscore(0);
        for(int j = 0; j < idxs.size(); j ++)
          if(idxs[j] != middle[0] && 0 <= idxs[j] && idxs[j] < corpust.cols()) {
            if(corpust(idxs[i], idxs[j])[middle[0]] != T(0))
              lscore --;
            if(corpust(idxs[j], idxs[i])[middle[0]] != T(0))
              lscore ++;
          }
        score.push_back(make_pair(lscore, idxs[i]));
      }
    sort(score.begin(), score.end());
    int i(0);
    for( ; i < score.size() && score[i].first < 0; i ++)
      left.push_back(score[i].second);
    for( ; i < score.size() && !score[i].first; i ++)
      middle.push_back(score[i].second);
    for( ; i < score.size(); i ++)
      right.push_back(score[i].second);
    if(left.size() || right.size())
      return serializeSub(left) + serializeSub(middle) + serializeSub(right);
  }
  return U("SYMMETRIC");
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::abbrev(const U& word, const corpushl<T, U>& mean) {
  cerr << "a" << flush;
  corpushl<T, U> work(match2relPseudo(mean));
  // XXX fixme: need to add abbreved word for work.
  const T tn(cdot(work)), td(work.cdot(work));
  if(isfinite(tn / td))
    return *this - (work * (tn / td));
  // can't abbrev.
  cerr << "XXX: can't abbrev : " << td << endl;
  // assert(td == T(0));
  return *this;
}

template <typename T, typename U> const vector<U> corpushl<T, U>::reverseLink(const corpushl<T, U>& orig) const {
  vector<U>      res;
  vector<int>    ridx0, ridx1;
  corpushl<T, U> work(match2relPseudo(corpushl<T, U>(orig)));
  vector<U>      rwords(gatherWords(words, work.words, ridx0, ridx1));
  return rwords;
}

template <typename T, typename U> const U corpushl<T, U>::reverseLink(const corpus<T, U>& orig) const {
  vector<int>    ridx0, ridx1;
  corpushl<T, U> work(match2relPseudo(orig));
  return orig.getAttributed(gatherWords(words, work.words, ridx0, ridx1));
}

template <typename T, typename U> const pair<T, T> corpushl<T, U>::compareStructure(const corpushl<T, U>& src, const T& thresh, const T& thresh2) const {
  // get H-SVD singular values for each of them and sort:
  const Vec s0(singularValues()), s1(src.singularValues());
  
  // get compared.
  pair<T, T> result;
  result.first = result.second = T(0);
  Mat S0(s0.size(), s0.size()), S1(s1.size(), s1.size());
  for(int i = 0; i < S0.rows(); i ++)
    for(int j = 0; j < S0.cols(); j ++) {
      S0(i, j) = s0[i] / s0[j];
      if(!isfinite(S0(i, j)) || T(1) / thresh < abs(S0(i, j)))
        S0(i, j) = T(0);
    }
  for(int i = 0; i < S1.rows(); i ++)
    for(int j = 0; j < S1.cols(); j ++) {
      S1(i, j) = s1[i] / s1[j];
      if(!isfinite(S1(i, j)) || T(1) / thresh < abs(S1(i, j)))
        S1(i, j) = T(0);
    }
  Eigen::JacobiSVD<Mat> svd0(S0, 0);
  Eigen::JacobiSVD<Mat> svd1(S1, 0);
  const Vec ss0(svd0.singularValues());
  const Vec ss1(svd1.singularValues());
  int i(0), j(0);
  for( ; i < ss0.size() && j < ss1.size(); )
    if(abs(ss0[i] - ss1[j]) / max(abs(ss0[i]), abs(ss1[i])) < thresh2) {
      result.first  += ss0[i] * ss0[i] + ss1[j] * ss1[j];
      i ++, j ++;
    } else {
      result.second += ss0[i] * ss0[i] + ss1[j] * ss1[j];
      if(ss0[i] > ss1[j])
        i ++;
      else
        j ++;
    }
  for( ; i < ss0.size(); i ++)
    result.second += ss0[i] * ss0[i];
  for( ; j < ss1.size(); j ++)
    result.second += ss1[j] * ss1[j];
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::reDig(const T& ratio) {
#if defined(_OPENMP)
#pragma omp paralell for schedule(static, 1)
#endif
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      for(int k = 0; k < corpust(i, j).size(); k ++)
        if(corpust(i, j)[k] != T(0))
          corpust(i, j)[k] = (corpust(i, j)[k] < T(0) ? - T(1) : T(1)) * exp(log(abs(corpust(i, j)[k])) * ratio);
  return *this;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::simpleThresh(const T& ratio) const {
  T absmax(0);
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      for(int k = 0; k < corpust(i, j).size(); k ++)
        if(absmax < abs(corpust(i, j)[k]))
          absmax = abs(corpust(i, j)[k]);
  Tensor work(corpust);
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      for(int k = 0; k < corpust(i, j).size(); k ++)
        if(abs(corpust(i, j)[k]) < absmax * ratio)
          work(i, j)[k] = T(0);
  vector<int> okidx;
  for(int i = 0; i < corpust.rows(); i ++) {
    bool flag(true);
    for(int j = 0; j < corpust.rows(); j ++) {
      if(work(i, j).dot(work(i, j)) != T(0) ||
         work(j, i).dot(work(j, i)) != T(0)) {
        flag = false;
        break;
      }
      for(int k = 0; k < corpust.cols(); k ++)
        if(i < work(j, k).size() && work(j, k)[i] != T(0)) {
          flag = false;
          break;
        }
      if(!flag)
        break;
    }
    if(!flag)
      okidx.push_back(i);
  }
  corpushl<T, U> result;
  result.words   = vector<U>();
  result.corpust = Tensor(okidx.size(), okidx.size());
  for(int i = 0; i < okidx.size(); i ++) {
    result.words.push_back(words[okidx[i]]);
    for(int j = 0; j < okidx.size(); j ++) {
      result.corpust(i, j) = Vec(okidx.size());
      for(int k = 0; k < okidx.size(); k ++)
        result.corpust(i, j)[k] = corpust(okidx[i], okidx[j])[okidx[k]];
    }
  }
  return result;
}

template <typename T, typename U> vector<U> corpushl<T, U>::gatherWords(const vector<U>& in0, const vector<U>& in1, vector<int>& ridx0, vector<int>& ridx1) const {
  cerr << "g" << flush;
  vector<U> result;
  ridx0 = vector<int>();
  ridx1 = vector<int>();
  if(!in0.size()) {
    if(!in1.size())
      return result;
    for(int i = 0; i < in1.size(); i ++) {
      ridx0.push_back(- 1);
      ridx1.push_back(i);
      result.push_back(in1[i]);
    }
    return result;
  }
  if(!in1.size()) {
    for(int i = 0; i < in0.size(); i ++) {
      ridx0.push_back(i);
      ridx1.push_back(- 1);
      result.push_back(in0[i]);
    }
    return result;
  }
  vector<U> sin0(in0), sin1(in1);
  sort(sin0.begin(), sin0.end());
  sort(sin1.begin(), sin1.end());
  int rbufsize(in0.size() + in1.size() + 1);
  for(int i = 0; i < rbufsize; i ++) {
    ridx0.push_back(- 1);
    ridx1.push_back(- 1);
  }
  int i(0), j(0);
  for( ; i < sin0.size(); i ++) {
    for(; i < sin0.size() && j < sin1.size(); j ++) {
      if(sin0[i] == sin1[j]) {
        result.push_back(U(sin0[i]));
        const int d0(distance(in0.begin(), equal_range(in0.begin(), in0.end(), sin0[i]).first));
        const int d1(distance(in1.begin(), equal_range(in1.begin(), in1.end(), sin1[j]).first));
        if(0 <= d0 && d0 < in0.size())
          ridx0[result.size() - 1] = d0;
        if(0 <= d1 && d1 < in1.size())
          ridx1[result.size() - 1] = d1;
        i ++;
        continue;
      } else if(sin0[i] > sin1[j]) {
        result.push_back(sin1[j]);
        const int d1(distance(in1.begin(), equal_range(in1.begin(), in1.end(), sin1[j]).first));
        if(0 <= d1 && d1 < in1.size())
          ridx1[result.size() - 1] = d1;
        continue;
      }
      break;
    }
    if(sin0.size() <= i) break;
    result.push_back(sin0[i]);
    const int d0(distance(in0.begin(), equal_range(in0.begin(), in0.end(), sin0[i]).first));
    if(0 <= d0 && d0 < in0.size())
      ridx0[result.size() - 1] = distance(in0.begin(), equal_range(in0.begin(), in0.end(), sin0[i]).first);
  }
  return result;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic> corpushl<T, U>::prepareDetail(const corpushl<T, U>& other, const vector<U>& workwords, const int& eidx, const vector<int>& ridx0, const vector<int>& ridx1, const vector<int>& ridx2) {
  cerr << "p" << flush;
  // initialize result tensor with 0.
  Tensor res(workwords.size() - 1, workwords.size() - 1);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < workwords.size(); i ++) if(ridx2[i] >= 0)
    for(int j = 0; j < workwords.size(); j ++) if(ridx2[j] >= 0) {
      res(ridx2[i], ridx2[j]) = Vec(workwords.size() - 1);
      for(int k = 0; k < res(ridx2[i], ridx2[j]).size(); k ++) if(ridx2[k] >= 0)
        res(ridx2[i], ridx2[j])[ridx2[k]] = T(0);
    }
  // Sum-up detailed word into result pool without definition row, col.
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < workwords.size(); i ++) if(ridx0[i] >= 0 && ridx2[i] >= 0)
    for(int j = 0; j < workwords.size(); j ++) if(ridx0[j] >= 0 && ridx2[j] >= 0)
      for(int k = 0; k < workwords.size(); k ++) if(ridx0[k] >= 0) {
        const T ratio(corpust(ridx0[i], ridx0[j])[ridx0[k]]);
        if(ridx2[k] >= 0 && ridx1[k] >= 0 && ridx1[i] >= 0 && ridx1[j] >= 0)
          res(ridx2[i], ridx2[j])[ridx2[k]] +=
            ratio * other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
      }

  // Sum-up words defined in definition row, col without crossing point.
  {
    const int i(eidx);
    for(int j = 0; j < workwords.size(); j ++) if(ridx2[j] >= 0) {
      for(int k = 0; k < workwords.size(); k ++) if(ridx2[k] >= 0) {
        // if all side exists on both corpust, weight and add.
        if(ridx0[i] >= 0 && ridx0[j] >= 0 && ridx0[k] >= 0 &&
           ridx1[j] >= 0 && ridx1[k] >= 0) {
          const T& r0(corpust(ridx0[i], ridx0[j])[ridx0[k]]);
          const T& r1(corpust(ridx0[j], ridx0[i])[ridx0[k]]);
          for(int kk = 0; kk < workwords.size(); kk ++)
            if(ridx1[kk] >= 0 && ridx2[kk] >= 0) {
              res(ridx2[kk], ridx2[j])[ridx2[k]] +=
                r0 * other.corpust(ridx1[kk], ridx1[j])[ ridx1[k]];
              res(ridx2[j], ridx2[kk])[ridx2[k]] +=
                r1 * other.corpust(ridx1[j], ridx1[kk])[ ridx1[k]];
            }
        // if only other have words defined, just add.
        } else if(ridx1[j] >= 0 && ridx1[k] >= 0 && ridx2[i] >= 0) {
          res(ridx2[i], ridx2[j])[ridx2[k]] += other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
          res(ridx2[j], ridx2[i])[ridx2[k]] += other.corpust(ridx1[j], ridx1[i])[ridx1[k]];
        }
      }
    }
  }
  return res;
}

template <typename T, typename U> const Eigen::Matrix<T, Eigen::Dynamic, 1> corpushl<T, U>::singularValues() const {
  Mat planes(corpust.rows(), corpust.cols());
  for(int i = 0; i < corpust.rows(); i ++) {
    Mat buf(corpust.rows(), corpust.cols());
    for(int j = 0; j < corpust.cols(); j ++) {
      for(int k = 0; k < corpust(i, j).size(); k ++)
        if(isfinite(corpust(i, j)[k]))
          buf(k, j) = corpust(i, j)[k];
        else {
          cerr << "nan" << flush;
          buf(k, j) = T(0);
        }
      for(int k = corpust(i, j).size(); k < buf.rows(); k ++)
        buf(k, j) = T(0);
    }
    Eigen::JacobiSVD<Mat> svd(buf, 0);
    planes.col(i) = svd.singularValues();
  }
  Eigen::JacobiSVD<Mat> svd(planes, 0);
  return svd.singularValues();
}

template <typename T, typename U> const vector<U>& corpushl<T, U>::getWords() const {
  return words;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic>& corpushl<T, U>::getCorpus() const {
  return corpust;
}



template <typename T, typename U> void getAbbreved(vector<corpushl<T, U> >& cstat, const vector<U>& words, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& delimiter, const int& szwindow) {
  assert(detailtitle.size() == detail.size());
  cerr << "getAbbreved : preparing datas";
  vector<vector<corpushl<T, U> > > details;
  for(int i = 0; i < detail.size(); i ++) {
    cerr << "." << flush;
    details.push_back(vector<corpushl<T, U> >());
    for(int j = 0; j < detail[i].size() / szwindow + 1; j ++) {
      corpus<T, U> lstat;
      lstat.init(words, 0, 120);
      lstat.compute(detail[i].substr(j * szwindow, szwindow).c_str(), delimiter);
      details[details.size() - 1].push_back(corpushl<T, U>(lstat));
    }
  }
  
  cerr << " getting abbrevs";
  for(int i = 0; i < details.size(); i ++) {
    cerr << "." << flush;
    for(int j = 0; j < cstat.size(); j ++)
      for(int k = 0; k < details[i].size(); k ++)
        cstat[j] = cstat[j].abbrev(detailtitle[i], details[i][k]);
  }
  cerr << endl;
  return;
}

template <typename T, typename U> void getDetailed(vector<corpus<T, U> >& cstat0, vector<corpushl<T, U> >& cstat, const U& input, const vector<U>& words, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& delimiter, const int& szwindow) {
  assert(detailtitle.size() == detail.size());
  cerr << "getDetailed : preparing datas";
  vector<vector<corpushl<T, U> > > details;
  for(int i = 0; i < detail.size(); i ++) {
    cerr << "." << flush;
    details.push_back(vector<corpushl<T, U> >());
    for(int j = 0; j < detail[i].size() / szwindow + 1; j ++) {
      corpus<T, U> lstat;
      lstat.init(words, 0, 120);
      lstat.compute(detail[i].substr(j * szwindow, szwindow).c_str(), delimiter);
      details[details.size() - 1].push_back(corpushl<T, U>(lstat));
    }
  }
  cstat0 = vector<corpus<T, U> >();
  cstat  = vector<corpushl<T, U> >();
  for(int i = 0; i < input.size() / szwindow + 1; i ++) {
    corpus<T, U> lstat;
    lstat.init(words, 0, 120);
    lstat.compute(input.substr(i * szwindow, szwindow).c_str(), delimiter);
    cstat0.push_back(lstat);
    cstat.push_back(corpushl<T, U>(lstat));
  }

  cerr << " getting details";
  for(int i = 0; i < details.size(); i ++) {
    cerr << "." << flush;
    for(int j = 0; j < cstat0.size(); j ++)
      for(int k = 0; k < details[i].size(); k ++)
        cstat[j] = cstat[j].withDetail(detailtitle[i], details[i][k]);
  }
  cerr << endl;
  return;
}

template <typename T, typename U> U preparedTOC(const U& input, const vector<U>& words, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& topictitle, const vector<U>& topics, const vector<U>& delimiter, const int& szwindow, const T& thresh, const T& redig = T(1)) {
  cerr << "preparedToc: parsing input" << endl;
  assert(detailtitle.size() == detail.size());
  assert(topictitle.size()  == topics.size());
  
  vector<corpus<T, U> >   cstat0;
  vector<corpushl<T, U> > cstat;
  getDetailed<T, U>(cstat0, cstat, input, words, detailtitle, detail, delimiter, szwindow);
  for(int i = 0; i < cstat0.size(); i ++)
    cstat[i].reDig(redig);
  
  cerr << "preparedToc: analysing input text" << flush;
  U result;
  for(int i = 0; i < topics.size(); i ++) {
     cerr << "." << flush;
     vector<corpus<T, U> >   tstat0;
     vector<corpushl<T, U> > tstat;
     getDetailed<T, U>(tstat0, tstat, topics[i], words, detailtitle, detail, delimiter, szwindow);
     vector<pair<T, pair<int, int> > > scores;
     for(int j = 0; j < tstat.size(); j ++)
       for(int k = 0; k < cstat.size(); k ++) {
         const auto score(- cstat[k].prej(tstat[j]));
         // XXX: quick and dirty.
         if(abs(score) < T(4))
           scores.push_back(make_pair(score, make_pair(j, k)));
       }
     sort(scores.begin(), scores.end());
     result += topictitle[i] + U(" : (") + to_string(scores[0].first);
     T   sum(0);
     int cnt(0);
     for(int j = 0; j < scores.size(); j ++) {
       sum += scores[j].first;
       cnt ++;
     }
     result += U(", ") + to_string(sum / cnt);
     result += U(", ") + to_string(scores[scores.size() - 1].first);
     result += U(")<br/><span class=\"small\">\n");;
     for(int j = 0; j < scores.size(); j ++) {
        if(- scores[j].first < thresh)
          break;
        const auto& work(cstat[scores[j].second.second] + tstat[scores[j].second.first]);
        result += to_string(scores[j].first) + U(" : ");
        result += work.serialize() + U("<br/>\n");
        result += work.reverseLink(cstat0[scores[j].second.second]) + U("<br/>\n");
        result += work.reverseLink(tstat0[scores[j].second.first])  + U("<br/><br/>\n");
     }
     result += U("</span><br/>\n");
  }
  return result;
}

template <typename T, typename U> U optimizeTOC(const U& input, const vector<U>& words, const vector<U>& detail, const vector<U>& detailtitle, const vector<U>& delimiter, const int& szwindow, const int& depth, const int& Mgather = 8, const T& redig = T(1)) {
  // prepare dictionaries:
  cerr << "optimizeToc: parsing input" << endl;
  vector<corpus<T, U> >   cstat0;
  vector<corpushl<T, U> > cstat;
  getDetailed<T, U>(cstat0, cstat, input, words, detailtitle, detail, delimiter, szwindow);
  
  for(int i = 0; i < cstat0.size(); i ++)
    cstat[i].reDig(redig);

  cerr << "optimizeToc: analysing input text." << flush;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> cstats(cstat.size(), cstat.size());
  for(int i = 0; i < cstat0.size(); i ++) {
    cerr << "." << flush;
    for(int j = 0; j < i; j ++)
      cstats(i, j) = cstats(j, i);
    cstats(i, i) = T(0);
    for(int j = i + 1; j < cstat.size(); j ++) {
      const auto score(- cstat[i].prej(cstat[j]));
      // XXX: quick and dirty.
      if(abs(score) < T(4))
        cstats(i, j) = score;
      else
        cstats(i, j) = T(8);
    }
  }
  
  cerr << "OK, sorting phrases." << flush;
  vector<int>           phrases;
  vector<pair<T, int> > work;
  vector<vector<int> >  idxs;
  auto                  cstatsw(cstats);
  for(int ii = 0; ii < cstat.size() && work.size() < cstat.size() / depth; ii ++) {
    vector<pair<T, int> > cidxs;
    for(int i = 0; i < cstatsw.rows(); i ++)
      if(!binary_search(phrases.begin(), phrases.end(), i)) {
        T lscore(0);
        for(int j = 0; j < cstatsw.rows(); j ++)
          lscore += cstatsw(j, i);
        cidxs.push_back(make_pair(lscore, i));
      }
    sort(cidxs.begin(), cidxs.end());
    const int& i(cidxs[0].second);
    phrases.push_back(i);
    sort(phrases.begin(), phrases.end());
    
    idxs.push_back(vector<int>());
    vector<pair<T, int> > scores;
    for(int j = 0; j < cstat.size(); j ++)
      if(!binary_search(phrases.begin(), phrases.end(), j))
        scores.push_back(make_pair(cstatsw(i, j), j));
    sort(scores.begin(), scores.end());
    T score(0);
    for(int j = 0; j < min(depth, int(scores.size())); j ++) {
      idxs[ii].push_back(scores[j].second);
      phrases.push_back(scores[j].second);
      score += scores[j].first;
      for(int k = 0; k < cstatsw.rows(); k ++)
        cstatsw(k, scores[j].second) = cstatsw(scores[j].second, k) = T(0);
    }
    sort(phrases.begin(), phrases.end());
    work.push_back(make_pair(score, i));
    for(int k = 0; k < cstatsw.rows(); k ++)
      cstatsw(k, i) = cstatsw(i, k) = T(0);
  }
  sort(work.begin(), work.end());
  
  cerr << "making outputs" << flush;
  U result;
  for(int jj = 0; jj < work.size(); jj ++) {
    const int&         j(work[jj].second);
    const vector<int>& idt(idxs[jj]);
    for(int k = 0; k < idt.size() / Mgather + 1; k ++) {
      corpushl<T, U> cs(cstat[j]);
      for(int l = k * Mgather; l < min((k + 1) * Mgather, int(idt.size())); l ++)
        cs += cstat[idt[l]];
      result += U("<div href=\"#\"><span class=\"small\">");
      result += to_string(work[jj].first) + U("<br/>");
      result += cs.serialize();
      result += U("</span><br/><span class=\"small\">");
      result += U("base : ") + to_string(j) + U(" - ");
      result += cs.reverseLink(cstat0[j]);
      result += U("<br/>");
      for(int l = k * Mgather; l < min((k + 1) * Mgather, int(idt.size())); l ++) {
        result += to_string(idt[l]) + U(" : ");
        result += to_string(cstats(j, idt[l])) + U(" - ");
        result += cs.reverseLink(cstat0[idt[l]]);
        result += U("<br/>");
      }
      result += U("</span></div><br/>");
    }
  }
  return result;
}

template <typename T, typename U> U diff(const U& input, const vector<U>& words, const vector<U>& detail0, const vector<U>& detailtitle0, const vector<U>& detail1, const vector<U>& detailtitle1, const vector<U>& delimiter, const int& szwindow, const T& thresh = T(0), const T& redig = T(1)) {
  cerr << "Diff: preparing inputs..." << endl;
  vector<corpus<T, U> >   cstat0;
  vector<corpus<T, U> >   dstat0;
  vector<corpushl<T, U> > cstat;
  vector<corpushl<T, U> > dstat;
  getDetailed<T, U>(cstat0, cstat, input, words, detailtitle0, detail0, delimiter, szwindow);
  getDetailed<T, U>(dstat0, dstat, input, words, detailtitle1, detail1, delimiter, szwindow);
  
  cerr << " making diffs" << endl;
  vector<corpushl<T, U> > diffs;
  for(int i = 0; i < cstat.size(); i ++) {
    cstat[i].reDig(redig);
    dstat[i].reDig(redig);
    auto diff(cstat[i] - dstat[i]);
    diff.reDig(redig);
    diffs.push_back(diff);
  }
  getAbbreved(diffs, words, detailtitle0, detail0, delimiter, szwindow);
  
  cerr << " making outputs." << flush;
  U result;
  for(int i = 0; i < cstat.size(); i ++)
    if(diffs[i].cdot(diffs[i]) != T(0)) {
      auto& diff(diffs[i]);
      diff    = diff.simpleThresh(thresh);
      result += diff.serialize() + U("<br/>\n");
      result += cstat[i].reverseLink(cstat0[i]) + U("<br/>\n");
      result += diff.reverseLink(cstat0[i]) + U("<br/>\n");
      result += diff.reverseLink(dstat0[i]) + U("<br/><br/><br/>\n");
    }
  result += U("Done");
  return result;
}

#define _CORPUS_
#endif


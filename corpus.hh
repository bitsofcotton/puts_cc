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

#include "sparse.hh"

using std::cerr;
using std::endl;
using std::flush;
using std::string;
using std::vector;
using std::sort;
using std::distance;
using std::upper_bound;
using std::binary_search;
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
using std::move;
using std::replace;

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

template <typename T> bool lessFirst(const T& a, const T& b) {
  return a.first < b.first;
}

template <typename T, typename U> class corpus {
public:
  typedef SimpleSparseVector<T> Vec;
  typedef SimpleSparseMatrix<T> Mat;
  typedef SimpleSparseTensor<T> Tensor;
  
  corpus();
  ~corpus();
  
  corpus<T, U>&    operator = (const corpus<T, U>& other);
  corpus<T, U>&    compute(const U& input, const vector<U>& delimiter, const vector<U>& words);
        U          getAttributed(const vector<U>& highlight) const;
  const U&         getOrig() const;
  const vector<U>& getWords() const;
  const Tensor&    getCorpus() const;
private:
  U                    orig;
  vector<U>            words;
  vector<vector<int> > ptrs;
  vector<int>          uptrs;
  vector<int>          pdelim;
  Tensor               corpust;
  int                  Midx;
   
  void getWordPtrs(const U& input, const vector<U>& delimiter);
  void corpusEach();
};

template <typename T, typename U> corpus<T,U>::corpus() {
  ;
}

template <typename T, typename U> corpus<T,U>::~corpus() {
  // auto called destructors for string.
  ;
}

template <typename T, typename U> corpus<T, U>& corpus<T, U>::operator = (const corpus<T, U>& other) {
  orig    = other.orig;
  words   = other.words;
  ptrs    = other.ptrs;
  uptrs   = other.uptrs;
  pdelim  = other.pdelim;
  corpust = other.corpust;
  Midx    = other.Midx;
  return *this;
}

template <typename T, typename U> corpus<T, U>& corpus<T,U>::compute(const U& input, const vector<U>& delimiter, const vector<U>& words) {
  this->words = words;
  cerr << "c" << flush;
  getWordPtrs(input, delimiter);
  corpusEach();
  return *this;
}

template <typename T, typename U> U corpus<T,U>::getAttributed(const vector<U>& highlight) const {
  U   result;
  int i;
  for(i = 0; i < orig.size(); ) {
    const auto lb(upper_bound(highlight.begin(), highlight.end(), U(&(orig.c_str()[i])), lessEqualStrClip<U>));
    if(highlight.begin() <= lb && lb < highlight.end() && equalStrClip<U>(*lb, U(&(orig.c_str()[i])))) {
      result += U("<font class=\"match\">");
      result += *lb;
      result += U("</font>");
      i      += lb->size();
    } else
      result += orig[i ++];
  }
  return result;
}

template <typename T, typename U> const U& corpus<T,U>::getOrig() const {
  return orig;
}

template <typename T, typename U> const vector<U>& corpus<T,U>::getWords() const {
  return words;
}

template <typename T, typename U> const SimpleSparseTensor<T>& corpus<T,U>::getCorpus() const {
  return corpust;
}

template <typename T, typename U> void corpus<T,U>::getWordPtrs(const U& input, const vector<U>& delimiter) {
  ptrs   = vector<vector<int> >();
  ptrs.resize(words.size(), vector<int>());
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
  orig = U(input);
  int i(0), i0(0);
  for( ; i < input.size(); i ++) {
    work += input[i];
    for(int ii = 0; ii < workd.size(); ii ++) {
      workd[ii]  = workd[ii].substr(1, workd[ii].size() - 1);
      workd[ii] += input[i];
      for(int j = 0; j < delimiter.size(); j ++)
        if(workd[ii] == delimiter[j] && pdelim[pdelim.size() - 1] < i)
          pdelim.push_back(i);
    }
    auto lo(upper_bound(words.begin(), words.end(), work, lessEqualStrClip<U>));
    auto up(upper_bound(words.begin(), words.end(), work, lessNotEqualStrClip<U>));
    bool match(false);
    for(auto itr(lo); itr < up; ++ itr)
      if(equalStrClip<U>(work, *itr)) {
        if(work.size() == itr->size()) {
          matchwidx.push_back(distance(words.begin(), itr));
          matchidxs.push_back(i0);
        } else if(work.size() < itr->size())
          match = true;
      }
    if(match && i < input.size() - 1)
      continue;
    if(matchwidx.size() > 0) {
      const int j(matchwidx.size() - 1);
      ptrs[matchwidx[j]].push_back(matchidxs[j]);
      uptrs.emplace_back(matchwidx[j]);
      Midx = matchidxs[j];
      matchwidx = vector<int>();
      matchidxs = vector<int>();
    }
    if(i == input.size() - 1)
      break;
    i   -= work.size() - 1;
    i0   = i;
    work = U();
  }
  pdelim.push_back(Midx + 2);
  const auto headidx(distance(words.begin(), lower_bound(words.begin(), words.end(), U("^"))));
  const auto tailidx(distance(words.begin(), lower_bound(words.begin(), words.end(), U("$"))));
  assert(0 <= headidx && headidx < words.size());
  assert(0 <= tailidx && tailidx < words.size());
  ptrs[headidx] = pdelim;
  ptrs[tailidx] = pdelim;
  uptrs.emplace_back(headidx);
  uptrs.emplace_back(tailidx);
  sort(uptrs.begin(), uptrs.end());
  uptrs.erase(unique(uptrs.begin(), uptrs.end()), uptrs.end());
  return;
}

template <typename T, typename U> void corpus<T,U>::corpusEach() {
  corpust = Tensor();
  for(auto itr0(uptrs.begin()); itr0 != uptrs.end(); ++ itr0) {
    const int i(*itr0);
    if(!ptrs[i].size()) continue;
    for(auto itr1(itr0); itr1 != uptrs.end(); ++ itr1) {
      const int j(*itr1);
      if(!ptrs[j].size()) continue;
      for(auto itr2(itr1); itr2 != uptrs.end(); ++ itr2) {
        const int k(*itr2);
        if(!ptrs[k].size()) continue;
        // XXX checkme:
        if(words[i] == U("$") || words[j] == U("$") || words[k] == U("$") ||
           words[i] == U("^") || words[j] == U("^") || words[k] == U("^"))
          continue;
        int ctru = 0;
        int ctrv = 0;
        int kk   = 0;
        for(auto itr = ptrs[k].begin(); itr != ptrs[k].end(); ++ itr) {
          while(ctru < ptrs[i].size() && ptrs[i][ctru] < *itr) ctru ++;
          ctru --;
          if(ctru < 0) ctru = 0;
          assert(0 <= ctru && ctru < ptrs[i].size());
          if(*itr <= ptrs[i][ctru])
            continue;
          while(ctrv < ptrs[j].size() && ptrs[j][ctrv] < *itr) ctrv ++;
          if(ptrs[j].size() <= ctrv || ptrs[j][ctrv] <= *itr)
            break;
          assert(0 <= ctrv && ctrv < ptrs[j].size());
          for( ; kk < pdelim.size() - 1; kk ++)
            if(pdelim[kk] <= *itr && *itr < pdelim[kk + 1])
              break;
          assert(0 <= kk && kk < pdelim.size());
          if(ptrs[i][ctru] < pdelim[kk] ||
             ptrs[j][ctrv] < pdelim[kk] ||
               pdelim[min(kk + 1, int(pdelim.size() - 1))] <= ptrs[i][ctru] ||
               pdelim[min(kk + 1, int(pdelim.size() - 1))] <= ptrs[j][ctrv])
            continue;
          // XXX configure me:
          const T buf0(log(T(abs(*itr + .5 - ptrs[i][ctru])) * T(2) * exp(T(1))));
          const T buf1(log(T(abs(*itr + .5 - ptrs[j][ctrv])) * T(2) * exp(T(1))));
          // const T buf0(abs(*itr + .5 - ptrs[i][ctru]));
          // const T buf1(abs(*itr + .5 - ptrs[j][ctrv]));
          const T work(T(1) / (buf0 * buf0 + buf1 * buf1));
          if(isfinite(work))
            corpust[i][j][k] += sqrt(work) / Midx;
        }
      }
    }
  }
  return;
}


template <typename T, typename U> class corpushl {
public:
  typedef SimpleSparseVector<T> Vec;
  typedef SimpleSparseMatrix<T> Mat;
  typedef SimpleSparseTensor<T> Tensor;
  
  corpushl();
  corpushl(const corpus<T, U>&   obj);
  corpushl(corpus<T, U>&&   obj);
  corpushl(const corpushl<T, U>& obj);
  corpushl(corpushl<T, U>&& obj);
  ~corpushl();
  
        corpushl<T, U>& operator += (const corpushl<T, U>& other);
        corpushl<T, U>& operator -= (const corpushl<T, U>& other);
        corpushl<T, U>& operator *= (const T& t);
        corpushl<T, U>& operator /= (const T& t);
        corpushl<T, U>  operator +  (const corpushl<T, U>& other) const;
        corpushl<T, U>  operator -  () const;
        corpushl<T, U>  operator -  (const corpushl<T, U>& other) const;
        corpushl<T, U>  operator *  (const T& t)                  const;
        corpushl<T, U>  operator /  (const T& t)                  const;
        corpushl<T, U>& operator =  (const corpushl<T, U>& other);
        corpushl<T, U>& operator =  (corpushl<T,U>&& other);
        bool            operator == (const corpushl<T, U>& other) const;
        bool            operator != (const corpushl<T, U>& other) const;
        corpushl<T, U>  withDetail(const U& word, const corpushl<T, U>& other, const T& thresh = T(0)) const;
        T               cdot(const corpushl<T, U>& other) const;
        T               absmax() const;
  const T               prej(const corpushl<T, U>& prejs) const;
  const T               prej2(const vector<corpushl<T, U> >& prej0, const vector<corpushl<T, U> >& prej1, const T& thresh) const;
        corpushl<T, U>& invertInsist();
  const corpushl<T, U>  conflictPart() const;
  const vector<U>&      getWords()  const;
  const Tensor&         getCorpus() const;
        U               serialize() const;
        corpushl<T, U>  abbrev(const U& word, const corpushl<T, U>& work, const T& thresh = T(0)) const;
        vector<U>       reverseLink() const;
        U               reverseLink(const corpus<T, U>& orig) const;
        pair<T, T>      compareStructure(const corpushl<T, U>& src, const T& thresh = T(1e-4), const T& thresh2 = T(.125)) const;
        corpushl<T, U>& reDig(const T& ratio);
        corpushl<T, U>  simpleThresh(const T& ratio) const;

private:
  U            serializeSub(const vector<int>& idxs) const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> singularValues() const;
  vector<int>  countIdx(const T& thresh = T(0)) const;
  void         merge5(Tensor& d, const int& i, const int& ki, const int& kk, const int& kj, const int& j, const T& intensity) const;
  
  vector<U>    words;
  Tensor       corpust;
};

template <typename T, typename U> corpushl<T,U>::corpushl() {
  // already initialized by compiler, do not need to initialize same.
/*
  words   = vector<U>();
  corpust = Tensor();
*/
  ;
}

template <typename T, typename U> corpushl<T,U>::~corpushl() {
  // auto called destructors for string.
  ;
}

template <typename T, typename U> corpushl<T,U>::corpushl(const corpus<T, U>& obj) {
  words   = vector<U>(obj.getWords());
  corpust = Tensor(obj.getCorpus());
  *this  /= cdot(*this);
}

template <typename T, typename U> corpushl<T,U>::corpushl(corpus<T, U>&& obj) {
  words   = move(obj.getWords());
  corpust = move(obj.getCorpus());
  *this  /= cdot(*this);
}

template <typename T, typename U> corpushl<T,U>::corpushl(const corpushl<T, U>& obj) {
  *this = obj;
}

template <typename T, typename U> corpushl<T,U>::corpushl(corpushl<T, U>&& obj) {
  *this = obj;
}

template <typename T, typename U> corpushl<T, U>& corpushl<T, U>::operator = (const corpushl<T, U>& other) {
  words   = vector<U>(other.words);
  corpust = Tensor(other.corpust);
  return *this;
}

template <typename T, typename U> corpushl<T, U>& corpushl<T, U>::operator = (corpushl<T, U>&& other) {
  words   = move(other.words);
  corpust = move(other.corpust);
  return *this;
}

template <typename T, typename U> bool corpushl<T, U>::operator == (const corpushl<T, U>& other) const {
  return ! (*this != other);
}

template <typename T, typename U> bool corpushl<T, U>::operator != (const corpushl<T, U>& other) const {
#if defined(_STRICT_WORD_ASSERT_)
  assert(words == other.words);
#else
  assert(words.size() == other.words.size());
#endif
  return corpust != other.corpust;
}

template <typename T, typename U> corpushl<T, U>& corpushl<T, U>::operator += (const corpushl<T, U>& other) {
#if defined(_STRICT_WORD_ASSERT_)
  assert(words == other.words);
#else
  assert(words.size() == other.words.size());
#endif
  corpust += other.corpust;
  return *this;
}

template <typename T, typename U> corpushl<T, U>& corpushl<T, U>::operator -= (const corpushl<T, U>& other) {
#if defined(_STRICT_WORD_ASSERT_)
  assert(words == other.words);
#else
  assert(words.size() == other.words.size());
#endif
  corpust -= other.corpust;
  return *this;
}

template <typename T, typename U> corpushl<T, U>& corpushl<T, U>::operator *= (const T& t) {
  corpust *= t;
  return *this;
}

template <typename T, typename U> corpushl<T, U>& corpushl<T, U>::operator /= (const T& t) {
  corpust /= t;
  return *this;
}

template <typename T, typename U> corpushl<T, U> corpushl<T, U>::operator + (const corpushl<T, U>& other) const {
  corpushl<T, U> result(*this);
  return result += other;
}

template <typename T, typename U> corpushl<T, U> corpushl<T, U>::operator - () const {
  corpushl<T, U> result(*this);
  result.corpust = - result.corpust;
  return result;
}

template <typename T, typename U> corpushl<T, U> corpushl<T, U>::operator - (const corpushl<T, U>& other) const {
  corpushl<T, U> result(*this);
  return result -= other;
}

template <typename T, typename U> corpushl<T, U> corpushl<T, U>::operator * (const T& t) const {
  corpushl<T, U> work(*this);
  return work *= t;
}

template <typename T, typename U> corpushl<T, U> corpushl<T, U>::operator / (const T& t) const {
  corpushl<T, U> work(*this);
  return work /= t;
}

template <typename T, typename U> corpushl<T, U> corpushl<T, U>::withDetail(const U& word, const corpushl<T, U>& other, const T& thresh) const {
#if defined(_STRICT_WORD_ASSERT_)
  assert(words == other.words);
#else
  assert(words.size() == other.words.size());
#endif
  const auto itr(lower_bound(words.begin(), words.end(), word));
  const int  eeidx(distance(words.begin(), itr));
  assert(0 <= eeidx && eeidx < words.size() && *itr == word);
  const auto idxs(countIdx(T(0)));
  if(!binary_search(idxs.begin(), idxs.end(), eeidx))
    return *this;
  cerr << "withDetail : " << word << endl;
  corpushl<T, U> result(*this + other);
  const T x0(corpust[eeidx][eeidx][eeidx]);
  const auto& ci0(other.corpust.iter());
  for(auto itr0(ci0.begin()); itr0 != ci0.end(); ++ itr0) {
    const auto& ci1(itr0->second.iter());
    const auto& ii(itr0->first);
    for(auto itr1(ci1.begin()); itr1 != ci1.end(); ++ itr1) {
      const auto& ci2(itr1->second.iter());
      const auto& jj(itr1->first);
      for(auto itr2(ci2.begin()); itr2 != ci2.end(); ++ itr2) {
        // Sum-up detailed word into result pool without definition row, col.
        const int& kk(itr2->first);
        if(itr2->second == T(0) || !(ii == eeidx || jj == eeidx || kk == eeidx))
          continue;
        // add crossing points
        const auto& ti0(corpust.iter());
        for(auto titr0(ti0.begin()); titr0 != ti0.end(); ++ titr0) {
          const auto& ti1(titr0->second.iter());
          const int& tii(titr0->first);
          if(tii == eeidx) continue;
          // add line points.
          for(auto titr1(ti1.begin()); titr1 != ti1.end(); ++ titr1) {
            const int& tjj(titr1->first);
            if(tjj == eeidx) continue;
            merge5(result.corpust, tii, ii, kk, jj, tjj, titr1->second[eeidx] * itr2->second * x0);
          }
        }
        for(auto titr0(ti0.begin()); titr0 != ti0.end(); ++ titr0) {
          const auto& ti2(titr0->second[eeidx].iter());
          const int& tii(titr0->first);
          if(tii == eeidx) continue;
          for(auto titr2(ti2.begin()); titr2 != ti2.end(); ++ titr2) {
            const int& tkk(titr2->first);
            if(tkk == eeidx) continue;
            merge5(result.corpust, tii, tkk, ii, kk, jj, titr2->second * itr2->second * x0);
          }
        }
        const auto& ti1(corpust[eeidx].iter());
        for(auto titr1(ti1.begin()); titr1 != ti1.end(); ++ titr1) {
          const auto& ti2(titr1->second.iter());
          const int& tjj(titr1->first);
          if(tjj == eeidx) continue;
          for(auto titr2(ti2.begin()); titr2 != ti2.end(); ++ titr2) {
            const int& tkk(titr2->first);
            if(tkk == eeidx) continue;
            merge5(result.corpust, ii, kk, jj, tjj, tkk, titr2->second * itr2->second * x0);
          }
        }
      }
    }
  }
  return result;
}

template <typename T, typename U> T corpushl<T, U>::cdot(const corpushl<T, U>& other) const {
#if defined(_STRICT_WORD_ASSERT_)
  assert(words == other.words);
#else
  assert(words.size() == other.words.size());
#endif
  T res(0);
  const auto& oi0(other.corpust.iter());
  for(auto itr0(oi0.begin()); itr0 != oi0.end(); ++ itr0)
    if(const_cast<const Tensor&>(corpust)[itr0->first].iter().size()) {
      const auto& oi1(itr0->second.iter());
      for(auto itr1(oi1.begin()); itr1 != oi1.end(); ++ itr1)
        if(const_cast<const Tensor&>(corpust)[itr0->first][itr1->first].iter().size()) {
          const auto& oi2(itr1->second.iter());
          for(auto itr2(oi2.begin()); itr2 != oi2.end(); ++ itr2)
            res += itr2->second * (const_cast<const Tensor&>(corpust))[itr0->first][itr1->first][itr2->first];
      }
  }
  return res;
}

template <typename T, typename U> T corpushl<T, U>::absmax() const {
  T res(0);
  const auto& ci0(corpust.iter());
  for(auto itr0(ci0.begin()); itr0 != ci0.end(); ++ itr0) {
    const auto& ci1(itr0->second.iter());
    for(auto itr1(ci1.begin()); itr1 != ci1.end(); ++ itr1) {
      const auto& ci2(itr1->second.iter());
      for(auto itr2(ci2.begin()); itr2 != ci2.end(); ++ itr2)
        res = max(res, abs(itr2->second));
    }
  }
  return res;
}

template <typename T, typename U> const T corpushl<T, U>::prej(const corpushl<T, U>& prejs) const {
#if defined(_STRICT_WORD_ASSERT_)
  assert(words == prejs.words);
#else
  assert(words.size() == prejs.words.size());
#endif
  static bool shown(false);
  if(!shown) {
    cerr << "XXX : confirm me corpushl::prej function." << endl;
    shown = true;
  }
  // XXX confirm me: need some other counting methods?
  const auto n2this(cdot(*this));
  if(n2this == T(0))
    return T(0);
  const auto n2p(prejs.cdot(prejs));
  if(n2p == T(0))
    return T(0);
  // XXX checkme : * 2 ratio.
  return cdot(prejs) / sqrt(n2this * n2p) * T(2);
}

template <typename T, typename U> const T corpushl<T, U>::prej2(const vector<corpushl<T, U> >& prej0, const vector<corpushl<T, U> >& prej1, const T& thresh) const {
  assert(words == prej0.words && words == prej1.words());
  static bool shown(false);
  if(!shown) {
    cerr << "XXX confirm me: corpushl::prej2" << endl;
    shown = true;
  }
  // XXX confirm me: is this correct counting method?
  corpushl<T, U> p0(*this), p1(*this);
  for(int i = 0; i < prej0.size(); i ++)
    p0 = p0.abbrev(string("P") + to_string(i), prej0[i]);
  for(int i = 0; i < prej1.size(); i ++)
    p1 = p1.abbrev(string("Q") + to_string(i), prej1[i]);
  p0 = p0.simpleThresh(thresh);
  p1 = p1.simpleThresh(thresh);
  return T(p0.words.size() - prej0.size()) / T(p1.words.size() - prej1.size());
}

template <typename T, typename U> corpushl<T, U>& corpushl<T, U>::invertInsist() {
  assert(0 && "confirm me: corpushl::invertInsist do not implemented NOT word table.");
  // XXX confirm me: this method cannot calculate in logically correct
  //                 because of it's method.
  return *this;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::conflictPart() const {
  assert(0 && "confirm me: corpushl::conflictPart do not implemented NOT word table.");
  // search conflict parts.
  // dictionary base of the word 'NOT' is needed.
  corpushl<T, U> result;
  return result;
}

template <typename T, typename U> U corpushl<T, U>::serialize() const {
  cerr << "s" << flush;
  auto plus(*this), minus(*this);
  const auto& pi0(plus.corpust.iter());
  for(auto itr0(pi0.begin()); itr0 != pi0.end(); ++ itr0) {
    const auto& pi1(itr0->second.iter());
    for(auto itr1(pi1.begin()); itr1 != pi1.end(); ++ itr1) {
      const auto& pi2(itr1->second.iter());
      for(auto itr2(pi2.begin()); itr2 != pi2.end(); ++ itr2)
        if(itr2->second < T(0))
          plus.corpust[itr0->first][itr1->first][itr2->first] = T(0);
    }
  }
  const auto& mi0(minus.corpust.iter());
  for(auto itr0(mi0.begin()); itr0 != mi0.end(); ++ itr0) {
    const auto& mi1(itr0->second.iter());
    for(auto itr1(mi1.begin()); itr1 != mi1.end(); ++ itr1) {
      const auto& mi2(itr1->second.iter());
      for(auto itr2(mi2.begin()); itr2 != mi2.end(); ++ itr2)
        if(itr2->second > T(0))
          minus.corpust[itr0->first][itr1->first][itr2->first] = T(0);
        else
          minus.corpust[itr0->first][itr1->first][itr2->first] = - itr2->second;
    }
  }
  const auto entire(countIdx(T(0)));
  return plus.serializeSub(entire)  + U(".&nbsp;&nbsp;&nbsp;-&nbsp;&nbsp;&nbsp;") +
         minus.serializeSub(entire) + U(".");
}

template <typename T, typename U> U corpushl<T, U>::serializeSub(const vector<int>& idxs) const {
  cerr << "." << flush;
  if(idxs.size() <= 1) {
    if(idxs.size())
      return words[idxs[0]];
    return U();
  }
  vector<pair<int, int> > cscore;
  cscore.reserve(idxs.size());
  // N.B. i0 - i1 - i2 is stored in corpust[i0][i2][i1].
  for(int i = 0; i < idxs.size(); i ++) {
    int lscore(0);
    for(int j = 0; j < idxs.size(); j ++) if((const_cast<const Tensor&>(corpust))[idxs[j]].iter().size())
      for(int k = 0; k < idxs.size(); k ++)
        if((const_cast<const Tensor&>(corpust))[idxs[j]][idxs[k]][idxs[i]] != T(0))
          lscore --;
    cscore.push_back(make_pair(lscore, idxs[i]));
  }
  sort(cscore.begin(), cscore.end());
  for(int si = 0; si < cscore.size(); si ++) {
    vector<int> middle, left, right;
    left.reserve(idxs.size());
    middle.reserve(idxs.size());
    right.reserve(idxs.size());
    middle.push_back(cscore[si].second);
    if(!cscore[si].first)
      goto symmetric;
    vector<pair<int, int> > score;
    for(int i = 0; i < idxs.size(); i ++) {
      if(idxs[i] == middle[0]) continue;
      int lscore(0);
      for(int j = 0; j < idxs.size(); j ++)
        if(idxs[j] != middle[0]) {
          if((const_cast<const Tensor&>(corpust))[idxs[i]][idxs[j]][middle[0]] != T(0))
            lscore --;
          if((const_cast<const Tensor&>(corpust))[idxs[j]][idxs[i]][middle[0]] != T(0))
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
    if((middle.size() && (left.size() || right.size())) || (left.size() && right.size()))
      return serializeSub(left) + serializeSub(middle) + serializeSub(right);
    // XXX checkme with some speed matter.
    break;
  }
 symmetric:
  U result;
  for(int i = 0; i < cscore.size(); i ++)
    if(cscore[i].first != T(0))
      result += words[cscore[i].second];
  return result;
}

template <typename T, typename U> corpushl<T, U> corpushl<T, U>::abbrev(const U& word, const corpushl<T, U>& work, const T& thresh) const {
#if defined(_STRICT_WORD_ASSERT_)
  assert(words == work.words);
#else
  assert(words.size() == work.words.size());
#endif
  const T tn(     cdot(work));
  const T td(work.cdot(work));
  if(td <= T(0))
    return *this;
  cerr << "abbrev: " << word << " : fixme ratio." << endl;
  auto result((*this * td - work * tn) / td);
  const int widx(distance(result.words.begin(), lower_bound(result.words.begin(), result.words.end(), word)));
  assert(0 <= widx && widx < result.words.size() && result.words[widx] == word);
  result.corpust[widx][widx][widx] += (tn < T(0) ? - T(1) : T(1)) * sqrt(abs(tn));
  const auto& ci0(result.corpust.iter());
  Mat c_ij, c_jk, c_ik;
  for(auto itr0(ci0.begin()); itr0 != ci0.end(); ++ itr0) {
    const auto& ci1(itr0->second.iter());
    for(auto itr1(ci1.begin()); itr1 != ci1.end(); ++ itr1) {
      const auto& ci2(itr1->second.iter());
      for(auto itr2(ci2.begin()); itr2 != ci2.end(); ++ itr2) {
        const int& i1(itr0->first);
        const int& j1(itr1->first);
        const int& k1(itr2->first);
        if(0 <= i1 && 0 <= j1)
          c_ij[i1][j1] += itr2->second;
        if(0 <= j1 && 0 <= k1)
          c_jk[j1][k1] += itr2->second;
        if(0 <= k1 && 0 <= i1)
          c_ik[i1][k1] += itr2->second;
      }
    }
  }
  const auto okidx(result.countIdx(T(0)));
  for(int i = 0; i < okidx.size(); i ++) {
    const auto ii(okidx[i]);
    if(ii == widx) continue;
    for(int j = 0; j < okidx.size(); j ++) {
      const auto jj(okidx[j]);
      if(jj == widx) continue;
      for(int k = 0; k < okidx.size(); k ++) {
        const auto kk(okidx[k]);
        if(kk == widx) continue;
        // XXX fixme ratio.
        const T score((const_cast<const Tensor&>(corpust))[ii][jj][kk] * (c_ij[ii][jj] + c_jk[jj][kk] + c_ik[ii][kk]) / result.words.size());
        result.corpust[widx][jj][kk] += score / T(3);
        result.corpust[ii][widx][kk] += score / T(3);
        result.corpust[ii][jj][widx] += score / T(3);
        result.corpust[ii][jj][kk]   -= score;
      }
    }
  }
  return result.simpleThresh(thresh);
}

template <typename T, typename U> vector<U> corpushl<T, U>::reverseLink() const {
  vector<U> res;
  const auto idx(countIdx(T(0)));
  res.reserve(idx.size());
  for(int i = 0; i < idx.size(); i ++)
    res.push_back(words[idx[i]]);
  return res;
}

template <typename T, typename U> U corpushl<T, U>::reverseLink(const corpus<T, U>& orig) const {
  return orig.getAttributed(reverseLink());
}

template <typename T, typename U> pair<T, T> corpushl<T, U>::compareStructure(const corpushl<T, U>& src, const T& thresh, const T& thresh2) const {
  // get H-SVD singular values for each of them and sort:
  const auto s0(singularValues()), s1(src.singularValues());
  
  // get compared.
  pair<T, T> result;
  result.first = result.second = T(0);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> S0(s0.size(), s0.size()), S1(s1.size(), s1.size());
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
  Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd0(S0, 0);
  Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd1(S1, 0);
  const auto ss0(svd0.singularValues());
  const auto ss1(svd1.singularValues());
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

template <typename T, typename U> corpushl<T, U>& corpushl<T, U>::reDig(const T& ratio) {
  auto& ci0(corpust.iter());
  for(auto itr0(ci0.begin()); itr0 != ci0.end(); ++ itr0) {
    auto& ci1(itr0->second.iter());
    for(auto itr1(ci1.begin()); itr1 != ci1.end(); ++ itr1) {
      auto& ci2(itr1->second.iter());
      for(auto itr2(ci2.begin()); itr2 != ci2.end(); ++ itr2)
        itr2->second = (itr2->second < T(0) ? - T(1) : T(1)) * exp(log(abs(itr2->second)) * ratio);
    }
  }
  return *this;
}

template <typename T, typename U> corpushl<T, U> corpushl<T, U>::simpleThresh(const T& ratio) const {
  assert(0 <= ratio);
  const auto thisabsmax(absmax());
  const auto okidx(countIdx(ratio * thisabsmax));
  corpushl<T, U> result;
  result.words   = words;
  result.corpust = Tensor();
  for(int i = 0; i < okidx.size(); i ++) {
    const auto ii(okidx[i]);
    if((const_cast<const Tensor&>(corpust))[okidx[i]].iter().size())
      for(int j = 0; j < okidx.size(); j ++) {
        const auto jj(okidx[j]);
        if((const_cast<const Tensor&>(corpust))[okidx[i]][okidx[j]].iter().size())
          for(int k = 0; k < okidx.size(); k ++) {
            const auto kk(okidx[k]);
            if(ratio * thisabsmax < abs((const_cast<const Tensor&>(corpust))[ii][jj][kk]))
              result.corpust[ii][jj][kk] = (const_cast<const Tensor&>(corpust))[ii][jj][kk];
          }
      }
  }
  return result;
}

template <typename T, typename U> vector<int> corpushl<T,U>::countIdx(const T& thresh) const {
  vector<int> okidx;
  const auto& ci0(corpust.iter());
  for(auto itr0(ci0.begin()); itr0 != ci0.end(); ++ itr0) {
    const auto& ci1(itr0->second.iter());
    for(auto itr1(ci1.begin()); itr1 != ci1.end(); ++ itr1) {
      const auto& ci2(itr1->second.iter());
      for(auto itr2(ci2.begin()); itr2 != ci2.end(); ++ itr2) {
        if(thresh < abs(itr2->second)) {
          okidx.push_back(itr0->first);
          okidx.push_back(itr1->first);
          okidx.push_back(itr2->first);
        }
      }
    }
  }
  sort(okidx.begin(), okidx.end());
  okidx.erase(unique(okidx.begin(), okidx.end()), okidx.end());
  return okidx;
}

template <typename T, typename U> void corpushl<T,U>::merge5(Tensor& d, const int& i, const int& ki, const int& kk, const int& kj, const int& j, const T& intensity) const {
  if(intensity == T(0)) return;
  d[ i][kk][ki] += intensity;
  d[ i][kj][ki] += intensity;
  d[ i][ j][ki] += intensity;
  d[ i][kj][kk] += intensity;
  d[ i][ j][kk] += intensity;
  d[ i][ j][kj] += intensity;
  d[ki][kj][kk] += intensity;
  d[ki][ j][kk] += intensity;
  d[ki][ j][kj] += intensity;
  d[kk][ j][kj] += intensity;
  return;
}

template <typename T, typename U> Eigen::Matrix<T, Eigen::Dynamic, 1> corpushl<T, U>::singularValues() const {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> planes(words.size(), words.size());
  for(int i = 0; i < words.size(); i ++) {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> buf(words.size(), words.size());
    for(int j = 0; j < words.size(); j ++) {
      for(int k = 0; k < words.size(); k ++)
        if(isfinite((const_cast<const Tensor&>(corpust))[i][j][k]))
          buf(k, j) = (const_cast<const Tensor&>(corpust))[i][j][k];
        else {
          cerr << "nan" << flush;
          buf(k, j) = T(0);
        }
      for(int k = words.size(); k < buf.rows(); k ++)
        buf(k, j) = T(0);
    }
    Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(buf, 0);
    planes.col(i) = svd.singularValues();
  }
  Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(planes, 0);
  return svd.singularValues();
}

template <typename T, typename U> const vector<U>& corpushl<T, U>::getWords() const {
  return words;
}

template <typename T, typename U> const SimpleSparseTensor<T>& corpushl<T, U>::getCorpus() const {
  return corpust;
}



template <typename T, typename U> void getAbbreved(corpushl<T, U>& cstat, const vector<U>& words, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& delimiter, const int& szwindow) {
  assert(detailtitle.size() == detail.size());
  for(int i = 0; i < detail.size(); i ++)
    cstat = cstat.abbrev(detailtitle[i], corpushl<T, U>(corpus<T, U>().compute(detail[i], delimiter, words)));
  return;
}

template <typename T, typename U> U getCut(const U& input, const int& idx, const int& szwindow) {
  if(input.size() < idx * szwindow / 2)
    return "XXX: getCut";
  return input.substr(idx * szwindow / 2, szwindow);
}

template <typename T, typename U> U getTagged(const U& name, const corpus<T, U>& cstat0, const corpushl<T, U>& cstat, const int& idx, const T& score, const U& input, const int& szwindow) {
  U tagged(U("<span id=\"") + name + to_string(idx) + U("\">"));
  tagged += U("score: ") + to_string(score) + U(" : ");
  tagged += getCut<T, U>(input, idx, szwindow);
  tagged += U("<input class=\"gatherdetail\" type=\"checkbox\"><div class=\"gatherdetail\">");
  tagged += cstat.reverseLink(cstat0) + U("<br/>");
  tagged += cstat0.getAttributed(cstat0.getWords());
  tagged += U("</div></span>");
  return tagged;
}

template <typename T, typename U> bool getDetailed(corpus<T, U>& cstat0, corpushl<T, U>& cstat, const U& input, const int& idx, const vector<U>& words, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& delimiter, const int& szwindow, const T& thresh) {
  assert(detailtitle.size() == detail.size());
  cstat0 = corpus<T, U>();
  cstat  = corpushl<T, U>();
  if(! (0 <= idx && idx < input.size() / szwindow * 2 + 1))
    return false;
  cstat0.compute(input.substr(idx * szwindow / 2, szwindow), delimiter, words);
  cstat  = corpushl<T, U>(cstat0).simpleThresh(thresh);
  for(int i = 0; i < detail.size(); i ++)
    cstat = cstat.withDetail(detailtitle[i], corpushl<T, U>(corpus<T, U>().compute(detail[i], delimiter, words)).simpleThresh(thresh), thresh).simpleThresh(thresh);
  return true;
}

template <typename T, typename U> U preparedTOC(const U& input, const vector<U>& words, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& topictitle, const vector<U>& topics, const vector<U>& delimiter, const int& szwindow, const int& depth, const T& threshin, const T& redig = T(1), const bool& reverse = false) {
  assert(detailtitle.size() == detail.size());
  assert(topictitle.size()  == topics.size());
  cerr << "preparedTOC..." << flush;
  U result;
  if(!topics.size())
    result += U("zero input.<br/>");
  vector<corpus<T, U> >   cistats;
  vector<corpushl<T, U> > istats;
  cistats.push_back(corpus<T, U>());
  istats.push_back(corpushl<T, U>());
  for(int j = 0; getDetailed<T, U>(cistats[cistats.size() - 1], istats[istats.size() - 1], input, j, words, detailtitle, detail, delimiter, szwindow, threshin); j ++) {
    istats[j].reDig(redig);
    cistats.push_back(corpus<T, U>());
    istats.push_back(corpushl<T, U>());
  }
  for(int i = 0; i < topics.size(); i ++) {
    vector<pair<T, pair<int, int> > > topicidx;
    vector<corpus<T, U> >   cstats;
    vector<corpushl<T, U> > stats;
    cstats.push_back(corpus<T, U>());
    stats.push_back(corpushl<T, U>());
    for(int j = 0; getDetailed<T, U>(cstats[cstats.size() - 1], stats[stats.size() - 1], topics[i], j, words, detailtitle, detail, delimiter, szwindow, threshin); j ++) {
      stats[j].reDig(redig);
      cstats.push_back(corpus<T, U>());
      stats.push_back(corpushl<T, U>());
    }
    for(int j = 0; j < istats.size() - 1; j ++) {
      const auto& cstat0(cistats[j]);
      const auto& stat0(istats[j]);
      int idx(0);
      T   score(0);
      for(int k = 0; k < stats.size() - 1; k ++) {
        const auto& cstat1(cstats[k]);
        const auto& stat1(stats[k]);
        const T lscore(reverse ? T(1) / abs(stat0.prej(stat1))
                               :            stat0.prej(stat1) );
        if(isfinite(lscore) && score <= lscore) {
          idx   = k;
          score = lscore;
        }
      }
      topicidx.push_back(make_pair(- score, make_pair(j, idx)));
    }
    sort(topicidx.begin(), topicidx.end());
    for(int j = 0; j < min(int(topicidx.size()), depth); j ++) {
      const auto& cstat0(cistats[topicidx[j].second.first]);
      const auto& cstat1(cstats[topicidx[j].second.second]);
      const auto& stat0(istats[topicidx[j].second.first]);
      const auto& stat1(stats[topicidx[j].second.second]);
      result += topictitle[i] + U(" ");
      result += getTagged<T,U>(U("prepTOC_") + to_string(i) + U("-") + to_string(j) + U("-1"), cstat0, stat1, topicidx[j].second.first, topicidx[j].first, input, szwindow);
      result += U("<br/>");
      result += getTagged<T,U>(U("prepTOC_") + to_string(i) + U("-") + to_string(j) + U("-2"), cstat1, stat0, topicidx[j].second.second, topicidx[j].first, topics[i], szwindow);
      result += U("<br/><br/>");
    }
    result += U("<br/><br/>");
  }
  return result;
}

template <typename T, typename U> U optimizeTOC(const U& input, const vector<U>& words, const vector<U>& detail, const vector<U>& detailtitle, const vector<U>& delimiter, const int& szwindow, const int& depth, const T& threshin, const T& redig = T(1), const bool& countnum = false, const U& notcheck = U("")) {
  assert(notcheck == U(""));
  cerr << "optimizeTOC..." << endl;
  SimpleSparseMatrix<T> scores;
  int Midx(0);
  vector<corpus<T, U> >   cstats;
  vector<corpushl<T, U> > stats;
  cstats.push_back(corpus<T, U>());
  stats.push_back(corpushl<T, U>());
  for(int i = 0; getDetailed<T, U>(cstats[cstats.size() - 1], stats[stats.size() - 1], input, i, words, detailtitle, detail, delimiter, szwindow, threshin); i ++) {
    stats[i].reDig(redig);
    cstats.push_back(corpus<T, U>());
    stats.push_back(corpushl<T, U>());
  }
  for(int i = 0; i < cstats.size() - 1; i ++) {
    const auto& stat0(stats[i]);
    for(int j = i + 1; j < cstats.size() - 1; j ++) {
      const auto& stat1(stats[j]);
      scores[i][j] = - stat0.prej(stat1);
      Midx = max(Midx, j);
    }
  }
  cerr << "OK" << flush;
  U   result;
  int idx(0);
  vector<int> phrases;
  for(int ii = 0; phrases.size() < Midx; ii ++) {
    vector<vector<pair<T, pair<int, int> > > > lscore;
    int  lidx(0);
    T    Mscore(0);
    bool fixed(false);
    for(int i = 0; i < Midx; i ++)
      if(!binary_search(phrases.begin(), phrases.end(), i)) {
        vector<pair<T, pair<int, int> > > llscore;
        for(int j = 0; j < Midx; j ++)
          if(!binary_search(phrases.begin(), phrases.end(), j))
            llscore.push_back(make_pair(scores[min(i, j)][max(i, j)], make_pair(i, j)));
        sort(llscore.begin(), llscore.end());
        lscore.push_back(llscore);
        bool ok(false);
        if(countnum) {
          ok = Mscore < llscore.size();
          Mscore = max(Mscore, T(llscore.size()));
        } else {
          T lllscore(0);
          for(int j = 0; j < min(depth, int(llscore.size())); j ++)
            lllscore -= llscore[j].first;
          ok = Mscore < lllscore;
          if(fixed)
            Mscore = max(Mscore, lllscore);
          else {
            Mscore = lllscore;
            fixed  = true;
          }
        }
        if(ok)
          lidx  = i;
      } else
        lscore.push_back(vector<pair<T, pair<int, int> > >());
    if(lscore[lidx].size() <= 0)
      break;
    phrases.push_back(lscore[lidx][0].second.first);
    sort(phrases.begin(), phrases.end());
    const auto& cstat0(cstats[lscore[lidx][0].second.first]);
    const auto& stat0(stats[lscore[lidx][0].second.first]);
    auto cs(stat0);
    for(int i = 0; i < min(depth, int(lscore[lidx].size())); i ++) {
      cs += stats[lscore[lidx][i].second.second];
      phrases.push_back(lscore[lidx][i].second.second);
    }
    sort(phrases.begin(), phrases.end());
    result += U("<form action=\"../../../../puts.php\" method=\"POST\"><div>");
    result += U("base : ");
    result += getTagged<T,U>(U("optTOC_") + to_string(lscore[lidx][0].second.first), cstat0, cs, lscore[lidx][0].second.first, lscore[lidx][0].first, input, szwindow) + U("<br/>");
    result += U("<br/>Show/Hide : <input class=\"gather\" type=\"checkbox\"><div class=\"gather\">");
    for(int i = 0; i < min(depth, int(lscore[lidx].size())); i ++) {
      const auto& cstat0(cstats[lscore[lidx][i].second.second]);
      const auto& stat0(stats[lscore[lidx][i].second.second]);
      result += getTagged<T,U>(U("optTOC_") + to_string(lscore[lidx][i].second.second), cstat0, cs, lscore[lidx][i].second.second, lscore[lidx][i].first, input, szwindow) + U("<br/>");
    }
    result += U("</div></div><textarea name=\"entry\">");
    result += getCut<T,U>(input, lscore[lidx][0].second.first, szwindow) + U("\n");
    for(int i = 0; i < min(depth, int(lscore[lidx].size())); i ++)
      result += getCut<T,U>(input, lscore[lidx][i].second.second, szwindow) + U("\n");
    result += U("</textarea>");
    result += U("<input type=\"hidden\" name=\"name\" value=\"append\" />");
    result += U("<input type=\"hidden\" name=\"adddict\" value=\"\" />");
    result += U("<input type=\"submit\" value=\"Append\" />");
    result += U("</form><br/>");
  }
  return result;
}

template <typename T, typename U> U diff(const U& input, const vector<U>& words, const vector<U>& detail0, const vector<U>& detailtitle0, const vector<U>& detail1, const vector<U>& detailtitle1, const vector<U>& delimiter, const int& szwindow, const T& threshin, const int& depth = T(20), const T& redig = T(1)) {
  cerr << "diff..." << flush;
  corpus<T, U>   cstat0, dstat0;
  corpushl<T, U> cstat,  dstat;
  U result;
  vector<pair<T, int> > scores;
  for(int i = 0; ; i ++) {
    if(!getDetailed<T, U>(cstat0, cstat, input, i, words, detailtitle0, detail0, delimiter, szwindow, threshin) ||
       !getDetailed<T, U>(dstat0, dstat, input, i, words, detailtitle1, detail1, delimiter, szwindow, threshin))
      break;
    getAbbreved<T, U>(cstat, words, detailtitle1, detail1, delimiter, szwindow);
    getAbbreved<T, U>(dstat, words, detailtitle0, detail0, delimiter, szwindow);
    cstat.reDig(redig);
    dstat.reDig(redig);
    const auto score(abs(cstat.cdot(dstat)) / sqrt(cstat.cdot(cstat) * dstat.cdot(dstat)) - T(1));
    if(isfinite(score))
      scores.push_back(make_pair(score, i));
  }
  sort(scores.begin(), scores.end());
  for(int ii = 0; ii < min(depth, int(scores.size())); ii ++) {
    const T&   score(scores[ii].first);
    const int& i(scores[ii].second);
    if(!getDetailed<T, U>(cstat0, cstat, input, i, words, detailtitle0, detail0, delimiter, szwindow, threshin) ||
       !getDetailed<T, U>(dstat0, dstat, input, i, words, detailtitle1, detail1, delimiter, szwindow, threshin))
      break;
    getAbbreved<T, U>(cstat, words, detailtitle1, detail1, delimiter, szwindow);
    getAbbreved<T, U>(dstat, words, detailtitle0, detail0, delimiter, szwindow);
    cstat.reDig(redig);
    dstat.reDig(redig);
    auto diff(cstat - dstat);
    diff.reDig(redig);
    result += U("(") + to_string(score) + U(") : ");
    result += diff.serialize() + U("<br/>\n");
    result += getTagged<T,U>(U("optTOC_") + to_string(i), cstat0, cstat, i, i, input, szwindow) + U("<br/>");
    result += getTagged<T,U>(U("optTOC_") + to_string(i), dstat0, dstat, i, i, input, szwindow) + U("<br/>");
    result += getTagged<T,U>(U("optTOC_") + to_string(i), cstat0, diff, i, i, input, szwindow) + U("<br/>");
    result += getTagged<T,U>(U("optTOC_") + to_string(i), dstat0, diff, i, i, input, szwindow) + U("<br/>");
  }
  return result;
}

#define _CORPUS_
#endif


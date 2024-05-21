/* BSD 3-Clause License:
 * Copyright (c) 2018-2021, bitsofcotton.
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

extern std::vector<std::string> words;

using std::vector;
using std::map;
using std::lower_bound;
using std::sort;
using std::unique;
using std::cerr;
using std::endl;
using std::flush;
using std::max;
using std::min;
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

template <typename T> static inline bool equalStrClip(const T& a, const T& b) {
  int cmp(0), jidx(0);
  for( ; !cmp && jidx < min(a.size(), b.size()); jidx ++)
    cmp = a[jidx] ^ b[jidx];
  return !cmp && min(a.size(), b.size()) <= jidx;
}

template <typename T> static inline bool lessEqualStrClip(const T& a, const T& b) {
  return a < b ||  equalStrClip<T>(a, b);
}

template <typename T> static inline bool lessNotEqualStrClip(const T& a, const T& b) {
  return a < b && !equalStrClip<T>(a, b);
}

template <typename T> class gram_t;
template <typename T> static inline bool lessCount(const gram_t<T>& dst, const gram_t<T>& src) {
  return (dst.rptr.size() < src.rptr.size()) ||
    (dst.rptr.size() == src.rptr.size() && dst.str < src.str);
}

template <typename T> class gram_t {
public:
  T           str;
  vector<int> rptr;
  inline gram_t() {
    this->str  = T();
    this->rptr = vector<int>();
  }
  inline ~gram_t() { ; }
  inline gram_t(const gram_t<T>& x) { *this = x; }
  inline gram_t& operator = (const gram_t<T>& x) {
    str  = x.str;
    rptr = x.rptr;
    return *this;
  }
  inline bool operator == (const gram_t<T>& x) const {
    return ! (*this != x);
  }
  inline bool operator != (const gram_t<T>& x) const {
    return str != x.str || rptr != x.rptr;
  }
  inline bool operator < (const gram_t<T>& x1) const {
    return str < x1.str;
  }
};

template <typename T, typename U> class lword {
public:
  inline lword(const int& loop) {
    this->dicts.resize(loop, vector<gram_t<U> >());
  }
  inline ~lword() { ; }
  
  vector<gram_t<U> > compute(const U& input);

private:
  inline bool       isin(const U& key) {
    assert(key.size() < dicts.size());
    const vector<gram_t<U> >& dict(dicts[key.size()]);
    gram_t<U> key0;
    key0.str = key;
    auto p(lower_bound(dict.begin(), dict.end(), key0));
    return dict.begin() <= p && p < dict.end() && p->str == key;
  }

  inline gram_t<U>& find(const U& key) {
    static gram_t<U> dummy;
    assert(key.size() < dicts.size());
    vector<gram_t<U> >& dict(dicts[key.size()]);
    gram_t<U> key0;
    key0.str = key;
    auto p(lower_bound(dict.begin(), dict.end(), key0));
    if(p < dict.begin() || dict.end() <= p || p->str != key) {
      assert(0 && "slipping find.");
      return dummy;
    }
    return *p;
  }

  inline void       assign(const gram_t<U>& val) {
    assert(val.str.size() < dicts.size());
    vector<gram_t<U> >& dict(dicts[val.str.size()]);
    auto p(lower_bound(dict.begin(), dict.end(), val));
    if(val.rptr.size()) {
      // delete duplicates:
      gram_t<U> work;
      work.str = val.str;
      auto& vptr(work.rptr = val.rptr);
      std::sort(vptr.begin(), vptr.end());
      vptr.erase(std::unique(vptr.begin(), vptr.end()), vptr.end());
      if(p < dict.begin() || dict.end() <= p || p->str != work.str) {
        dict.emplace_back(work);
        sort(dict.begin(), dict.end());
      } else
        *p = work;
    } else if(dict.begin() <= p && p < dict.end() && p->str == val.str)
      dict.erase(p);
    return;
  }

  vector<T>                   dict0;
  vector<vector<gram_t<U> > > dicts;
};

template <typename T, typename U> vector<gram_t<U> > lword<T, U>::compute(const U& input) {
  cerr << "lword(" << dicts.size() << ", " << input.size() << ")" << endl;
  // bigram
  map<U, vector<int> > mapw;
  for(int i = 1; i < input.size(); i ++) {
    U work;
    work += T(input[i - 1]);
    work += T(input[i]);
    mapw[work].emplace_back(i);
    if(!binary_search(dict0.begin(), dict0.end(), input[i])) {
      dict0.emplace_back(input[i]);
      sort(dict0.begin(), dict0.end());
    }
  }
  for(auto itr = mapw.begin(); itr != mapw.end(); ++ itr) {
    gram_t<U> work;
    work.str = itr->first;
    work.rptr.insert(work.rptr.end(), itr->second.begin(), itr->second.end());
    assign(work);
  }
  // construct
  int i;
  for(i = 2; i < dicts.size() - 1; i ++) {
    map<U, vector<int> > amap;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(auto itr = dicts[i].begin(); itr < dicts[i].end(); ++ itr) {
      const gram_t<U>& idxkey(*itr);
      U key2;
      for(int j = 1; j < idxkey.str.size(); j ++)
        key2 += idxkey.str[j];
      key2 += T(' ');
      for(auto itr2 = dict0.begin(); itr2 != dict0.end(); ++ itr2) {
        key2[key2.size() - 1] = T(*itr2);
        if(!isin(key2))
          continue;
        const gram_t<U>& idxkey2(find(key2));
        U workkey(idxkey.str);
        workkey += *itr2;
        vector<int> idxwork;
        int tt = 0;
        int ss = 0;
        while(tt < idxkey2.rptr.size() && ss < idxkey.rptr.size()) {
          if(idxkey2.rptr[tt] == idxkey.rptr[ss] + 1) {
            idxwork.emplace_back(idxkey.rptr[ss]);
            ss ++;
            tt ++;
          } else if(idxkey2.rptr[tt] >= idxkey.rptr[ss])
            ss ++;
          else
            tt ++;
        }
        assert(idxwork.size() <= min(idxkey.rptr.size(), idxkey2.rptr.size()));
        if(idxwork.size() < 2) continue;
#if defined(_OPENMP)
#pragma omp critical
#endif
        {
          amap[workkey].insert(amap[workkey].end(), idxwork.begin(), idxwork.end());
        }
      }
    }
    // construct next stage.
    for(auto itr = amap.begin(); itr != amap.end(); ++ itr) {
      gram_t<U> work;
      work.str  = itr->first;
      work.rptr.insert(work.rptr.end(), itr->second.begin(), itr->second.end());
      assign(work);
    }
  }
  
  // count up.
  vector<gram_t<U> > words;
  for(int i = 0; i < dicts.size(); i ++)
    words.insert(words.end(), dicts[i].begin(), dicts[i].end());
  dict0 = vector<T>();
  dicts = vector<vector<gram_t<U> > >();
  sort(words.begin(), words.end());
  vector<gram_t<U> > result;
  result.reserve(words.size());
  for(int i = 0; i < words.size() - 1; i ++) {
    const auto mw(min(words[i].str.size(), words[i + 1].str.size()));
    if(! (words[i].str.substr(0, mw) == words[i + 1].str.substr(0, mw) &&
          words[i].rptr == words[i + 1].rptr) ) {
      if(words[i].str.size() < words[i + 1].str.size() &&
         words[i].str.substr(1, mw) == words[i + 1].str.substr(0, mw) &&
         words[i].rptr.size() == words[i + 1].rptr.size()) {
        int j;
        for(j = 0; j < words[i].rptr.size(); j ++)
          if(words[i].rptr[j] != words[i + 1].rptr[j] + 1)
            break;
        if(j != words[i].rptr.size())
          result.emplace_back(move(words[i]));
      } else if(words[i].str.size() > words[i + 1].str.size() &&
         words[i].str.substr(0, mw) == words[i + 1].str.substr(1, mw) &&
         words[i].rptr.size() == words[i + 1].rptr.size()) {
        int j;
        for(j = 0; j < words[i].rptr.size(); j ++)
          if(words[i].rptr[j] + 1 != words[i + 1].rptr[j])
            break;
        if(j != words[i].rptr.size())
          result.emplace_back(move(words[i]));
      } else
        result.emplace_back(move(words[i]));
   }
  }
  result.emplace_back(move(words[words.size() - 1]));
  return result;
}


template <typename T, typename U> class corpus {
public:
  typedef SimpleSparseVector<T> Vec;
  typedef SimpleSparseMatrix<T> Mat;
  typedef SimpleSparseTensor<T> Tensor;
  
  corpus(const U& input, const vector<U>& delimiter);
  
  inline corpus() { ; }
  inline corpus(const corpus<T, U>& other) { *this = other; }
  inline corpus(corpus<T, U>&& other) { *this = other; }
  inline ~corpus() { ; }
  
  inline U getAttributed(const vector<U>& highlight) const {
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
  inline pair<vector<U>, U> reverseLink() const {
    pair<vector<U>, U> res;
    const auto idx(countIdx(T(0)));
    res.first.reserve(idx.size());
    for(int i = 0; i < idx.size(); i ++)
      res.first.emplace_back(words[idx[i]]);
    res.second = getAttributed(res.first);
    return res;
  }
  inline corpus<T, U>& operator += (const corpus<T, U>& other) {
    orig    += U("+") + other.orig;
    corpust += other.corpust;
    return *this;
  }
  inline corpus<T, U>& operator -= (const corpus<T, U>& other) {
    orig    += U("-") + other.orig;
    corpust -= other.corpust;
    return *this;
  }
  inline corpus<T, U>& operator *= (const T& t) {
    orig    += U("*") + U(to_string(t));
    corpust *= t;
    return *this;
  }
  inline corpus<T, U>& operator /= (const T& t) {
    orig    += U("/") + U(to_string(t));
    corpust /= t;
    return *this;
  }
  inline corpus<T, U>  operator +  (const corpus<T, U>& other) const {
    auto result(*this);
    return result += other;
  }
  inline corpus<T, U>  operator -  () const {
    auto result(*this);
    result.orig    = U("-") + result.orig;
    result.corpust = - result.corpust;
    return result;
  }
  inline corpus<T, U>  operator -  (const corpus<T, U>& other) const {
    auto result(*this);
    return result -= other;
  }
  inline corpus<T, U>  operator *  (const T& t)                  const {
    auto result(*this);
    return result *= t;
  }
  inline corpus<T, U>  operator /  (const T& t)                  const {
    auto result(*this);
    return result /= t;
  }
  inline corpus<T, U>& operator =  (const corpus<T, U>& other) {
    corpust = other.corpust;
    orig    = other.orig;
    return *this;
  }
  inline corpus<T, U>& operator =  (corpus<T, U>&& other) {
    corpust = move(other.corpust);
    orig    = move(other.orig);
    return *this;
  }
  inline bool          operator == (const corpus<T, U>& other) const {
    return ! (*this != other);
  }
  inline bool          operator != (const corpus<T, U>& other) const {
    return corpust != other.corpust;
  }
  T             cdot(const corpus<T, U>& other) const {
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
  T             absmax() const {
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
  corpus<T, U>& reDig(const T& ratio) {
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
  corpus<T, U> simpleThresh(const T& ratio) const {
    assert(0 <= ratio);
    const auto thisabsmax(absmax());
    const auto okidx(countIdx(ratio * thisabsmax));
    corpus<T, U> result;
    result.orig = orig;
    for(int i = 0; i < okidx.size(); i ++) {
      const auto& ii(okidx[i]);
      if((const_cast<const Tensor&>(corpust))[okidx[i]].iter().size())
        for(int j = 0; j < okidx.size(); j ++) {
          const auto& jj(okidx[j]);
          if((const_cast<const Tensor&>(corpust))[okidx[i]][okidx[j]].iter().size())
            for(int k = 0; k < okidx.size(); k ++) {
              const auto& kk(okidx[k]);
              if(ratio * thisabsmax < abs((const_cast<const Tensor&>(corpust))[ii][jj][kk]))
                result.corpust[ii][jj][kk] = (const_cast<const Tensor&>(corpust))[ii][jj][kk];
            }
        }
    }
    return result;
  }
  inline SimpleVector<T> singularValues(const SimpleMatrix<T>& m) const {
    const auto SV(m.SVD() * m);
    SimpleVector<T> w(SV.rows());
    for(int i = 0; i < w.size(); i ++)
      w[i] = sqrt(SV.row(i).dot(SV.row(i)));
    return w;
  }
  vector<int>  countIdx(const T& thresh = T(0)) const {
    vector<int> okidx;
    const auto& ci0(corpust.iter());
    for(auto itr0(ci0.begin()); itr0 != ci0.end(); ++ itr0) {
      const auto& ci1(itr0->second.iter());
      for(auto itr1(ci1.begin()); itr1 != ci1.end(); ++ itr1) {
        const auto& ci2(itr1->second.iter());
        for(auto itr2(ci2.begin()); itr2 != ci2.end(); ++ itr2) {
          if(thresh < abs(itr2->second)) {
            okidx.emplace_back(itr0->first);
            okidx.emplace_back(itr1->first);
            okidx.emplace_back(itr2->first);
          }
        }
      }
    }
    sort(okidx.begin(), okidx.end());
    okidx.erase(unique(okidx.begin(), okidx.end()), okidx.end());
    return okidx;
  }
  const T       prej(const corpus<T, U>& prejs) const;
  const T       prej2(const vector<corpus<T, U> >& prej0, const vector<corpus<T, U> >& prej1, const T& thresh) const;
  corpus<T, U>& invertInsist();
  corpus<T, U>  conflictPart() const;
  U             serialize() const;
  corpus<T, U>  withDetail(const U& word, const corpus<T, U>& other, const T& thresh = T(0)) const;
  corpus<T, U>  abbrev(const U& word, const corpus<T, U>& work, const T& thresh = T(0)) const;
  pair<T, T>    compareStructure(const corpus<T, U>& src, const T& thresh = T(1e-4), const T& thresh2 = T(.125)) const;
  corpus<T, U>& absfy();

  Tensor corpust;
private:
  SimpleVector<T> singularValues() const;
  U     serializeSub(const vector<int>& idxs) const;
  void  merge5(Tensor& d, const int& i, const int& ki, const int& kk, const int& kj, const int& j, const T& intensity) const;
  
  U     orig;
};

template <typename T, typename U> corpus<T,U>::corpus(const U& input, const vector<U>& delimiter) {
  // get word ptrs.
  vector<vector<int> > ptrs;
  vector<int>          uptrs;
  vector<int>          pdelim;
  ptrs.resize(words.size(), vector<int>());
  pdelim.emplace_back(0);
  U work;
  vector<int> matchwidx;
  vector<int> matchidxs;
  int dM(0);
  for(int i = 0; i < delimiter.size(); i ++)
    dM = max(dM, int(delimiter[i].size()));
  vector<U> workd;
  for(int i = 0; i < dM; i ++) {
    workd.emplace_back(U(""));
    for(int j = i; j < dM; j ++)
      workd[i] += U(" ");
  }
  orig = U(input);
  int i(0), i0(0), Midx(0);
  for( ; i < orig.size(); i ++) {
    work += orig[i];
    for(int ii = 0; ii < workd.size(); ii ++) {
      workd[ii]  = workd[ii].substr(1, workd[ii].size() - 1);
      workd[ii] += orig[i];
      for(int j = 0; j < delimiter.size(); j ++)
        if(workd[ii] == delimiter[j] && pdelim[pdelim.size() - 1] < i)
          pdelim.emplace_back(i);
    }
    auto lo(upper_bound(words.begin(), words.end(), work, lessEqualStrClip<U>));
    auto up(upper_bound(words.begin(), words.end(), work, lessNotEqualStrClip<U>));
    bool match(false);
    for(auto itr(lo); itr < up; ++ itr)
      if(equalStrClip<U>(work, *itr)) {
        if(work.size() == itr->size()) {
          matchwidx.emplace_back(distance(words.begin(), itr));
          matchidxs.emplace_back(i0);
        } else if(work.size() < itr->size())
          match = true;
      }
    if(match && i < orig.size() - 1)
      continue;
    if(matchwidx.size() > 0) {
      const int j(matchwidx.size() - 1);
      ptrs[matchwidx[j]].emplace_back(matchidxs[j]);
      uptrs.emplace_back(matchwidx[j]);
      Midx = matchidxs[j];
      matchwidx.resize(0);
      matchidxs.resize(0);
    }
    i0 = (i -= work.size() - 1) + 1;
    if(i == orig.size() - 1)
      break;
    work = U();
  }
  pdelim.emplace_back(Midx + 2);
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
  
  // corpus each
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
          if(ptrs[i][ctru] < pdelim[kk] && pdelim[kk] <= ptrs[j][ctrv])
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
}

template <typename T, typename U> const T corpus<T, U>::prej(const corpus<T, U>& prejs) const {
  static bool shown(false);
  if(!shown) {
    cerr << "XXX : confirm me corpus::prej function." << endl;
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

template <typename T, typename U> const T corpus<T, U>::prej2(const vector<corpus<T, U> >& prej0, const vector<corpus<T, U> >& prej1, const T& thresh) const {
  static bool shown(false);
  if(!shown) {
    cerr << "XXX confirm me: corpus::prej2" << endl;
    shown = true;
  }
  // XXX confirm me: is this correct counting method?
  corpus<T, U> p0(*this), p1(*this);
  for(int i = 0; i < prej0.size(); i ++)
    p0 = p0.abbrev(string("P") + to_string(i), prej0[i]);
  for(int i = 0; i < prej1.size(); i ++)
    p1 = p1.abbrev(string("Q") + to_string(i), prej1[i]);
  p0 = p0.simpleThresh(thresh);
  p1 = p1.simpleThresh(thresh);
  return T(words.size() - prej0.size()) / T(words.size() - prej1.size());
}

template <typename T, typename U> corpus<T, U>& corpus<T, U>::invertInsist() {
  assert(0 && "confirm me: corpus::invertInsist do not implemented NOT word table.");
  // XXX confirm me: this method cannot calculate in logically correct
  //                 because of it's method.
  return *this;
}

template <typename T, typename U> corpus<T, U> corpus<T, U>::conflictPart() const {
  assert(0 && "confirm me: corpus::conflictPart do not implemented NOT word table.");
  // search conflict parts.
  // dictionary base of the word 'NOT' is needed.
  corpus<T, U> result;
  return result;
}

template <typename T, typename U> corpus<T, U>& corpus<T, U>::absfy() {
  const auto& pi0(corpust.iter());
  for(auto itr0(pi0.begin()); itr0 != pi0.end(); ++ itr0) {
    const auto& pi1(itr0->second.iter());
    for(auto itr1(pi1.begin()); itr1 != pi1.end(); ++ itr1) {
      const auto& pi2(itr1->second.iter());
      for(auto itr2(pi2.begin()); itr2 != pi2.end(); ++ itr2)
        if(itr2->second < T(0)) {
          // N.B. we estimate minus sign on the tensor as to be reverse order.
          corpust[itr1->first][itr0->first][itr2->first] -= itr2->second;
          corpust[itr0->first][itr1->first][itr2->first]  = T(int(0));
        }
    }
  }
  return *this = simpleThresh(T(int(0)));
}

template <typename T, typename U> U corpus<T, U>::serialize() const {
  cerr << "s" << flush;
  auto plus(*this);
  return plus.absfy().serializeSub(plus.countIdx(T(int(0))));
}

template <typename T, typename U> U corpus<T, U>::serializeSub(const vector<int>& idxs) const {
  cerr << "." << flush;
  if(idxs.size() <= 1) {
    if(idxs.size())
      return words[idxs[0]];
    return U();
  }
  vector<pair<int, int> > score;
  score.reserve(idxs.size());
  // N.B. i0 - i1 - i2 is stored in corpust[i0][i2][i1].
  for(int i = 0; i < idxs.size(); i ++) {
    int lscore(0);
    for(int j = 0; j < idxs.size(); j ++)
      if(const_cast<const Tensor&>(corpust)[idxs[j]].iter().size()) {
        for(int k = 0; k < idxs.size(); k ++)
          if(const_cast<const Tensor&>(corpust)[idxs[j]][idxs[i]][idxs[k]] != T(0))
            lscore --;
      }
    const auto& ii(const_cast<const Tensor&>(corpust)[idxs[i]].iter());
    if(ii.size()) for(int j = 0; j < idxs.size(); j ++)
      if(const_cast<const Tensor&>(corpust)[idxs[i]][idxs[j]].iter().size()) {
        for(int k = 0; k < idxs.size(); k ++)
          if(const_cast<const Tensor&>(corpust)[idxs[i]][idxs[j]][idxs[k]] != T(0))
            lscore ++;
      }
    // XXX: middle data ignored.
    score.emplace_back(make_pair(lscore, idxs[i]));
  }
  sort(score.begin(), score.end());
  vector<int> left, right;
  left.reserve(idxs.size());
  right.reserve(idxs.size());
  int i(0);
  // N.B.: we only focus orders itself.
  for( ; i < score.size() / 2; i ++)
    left.emplace_back(score[i].second);
  for( ; i < score.size(); i ++)
    right.emplace_back(score[i].second);
  return serializeSub(left) + serializeSub(right);
}

template <typename T, typename U> pair<T, T> corpus<T, U>::compareStructure(const corpus<T, U>& src, const T& thresh, const T& thresh2) const {
  // get H-SVD singular values for each of them and sort:
  const auto s0(singularValues()), s1(src.singularValues());
  
  // get compared.
  pair<T, T> result;
  result.first = result.second = T(0);
  SimpleMatrix<T> S0(s0.size(), s0.size());
  SimpleMatrix<T> S1(s1.size(), s1.size());
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
  const auto ss0(singularValues(S0));
  const auto ss1(singularValues(S1));
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

template <typename T, typename U> SimpleVector<T> corpus<T, U>::singularValues() const {
  SimpleMatrix<T> planes(words.size(), words.size());
  for(int i = 0; i < words.size(); i ++) {
    SimpleMatrix<T> buf(words.size(), words.size());
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
    planes.col(i) = singularValues(buf);
  }
  return singularValues(planes);
}

template <typename T, typename U> corpus<T, U> corpus<T, U>::withDetail(const U& word, const corpus<T, U>& other, const T& thresh) const {
  const auto itr(lower_bound(words.begin(), words.end(), word));
  const int  eeidx(distance(words.begin(), itr));
  assert(0 <= eeidx && eeidx < words.size() && *itr == word);
  const auto idxs(countIdx(T(0)));
  if(!binary_search(idxs.begin(), idxs.end(), eeidx))
    return *this;
  cerr << "withDetail : " << word << endl;
  corpus<T, U> result(*this);
  // XXX:
  // corpus<T, U> result(*this + other);
  const T x0(const_cast<const Tensor&>(corpust)[eeidx][eeidx][eeidx]);
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
            merge5(result.corpust, tii, ii, kk, jj, tjj, titr1->second[eeidx] * itr2->second * (x0 + T(1)));
          }
        }
        for(auto titr0(ti0.begin()); titr0 != ti0.end(); ++ titr0) {
          const auto& ti2(const_cast<const Mat&>(titr0->second)[eeidx].iter());
          const int& tii(titr0->first);
          if(tii == eeidx) continue;
          for(auto titr2(ti2.begin()); titr2 != ti2.end(); ++ titr2) {
            const int& tkk(titr2->first);
            if(tkk == eeidx) continue;
            merge5(result.corpust, tii, tkk, ii, kk, jj, titr2->second * itr2->second * (x0 + T(1)));
          }
        }
        const auto& ti1(const_cast<const Tensor&>(corpust)[eeidx].iter());
        for(auto titr1(ti1.begin()); titr1 != ti1.end(); ++ titr1) {
          const auto& ti2(titr1->second.iter());
          const int& tjj(titr1->first);
          if(tjj == eeidx) continue;
          for(auto titr2(ti2.begin()); titr2 != ti2.end(); ++ titr2) {
            const int& tkk(titr2->first);
            if(tkk == eeidx) continue;
            merge5(result.corpust, ii, kk, jj, tkk, tjj, titr2->second * itr2->second * (x0 + T(1)));
          }
        }
      }
    }
  }
  return result;
}

template <typename T, typename U> corpus<T, U> corpus<T, U>::abbrev(const U& word, const corpus<T, U>& work, const T& thresh) const {
  const T tn(     cdot(work));
  const T td(work.cdot(work));
  if(td <= T(0))
    return *this;
  cerr << "abbrev: " << word << " : fixme ratio." << endl;
  auto result(*this);
  // XXX:
  //auto result((*this * td - work * tn) / td);
  const int widx(distance(words.begin(), lower_bound(words.begin(), words.end(), word)));
  assert(0 <= widx && widx < words.size() && words[widx] == word);
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
    const auto& ii(okidx[i]);
    if(ii == widx) continue;
    for(int j = 0; j < okidx.size(); j ++) {
      const auto& jj(okidx[j]);
      if(jj == widx) continue;
      for(int k = 0; k < okidx.size(); k ++) {
        const auto& kk(okidx[k]);
        if(kk == widx) continue;
        const auto  denom(c_ij[ii][jj] + c_jk[jj][kk] + c_ik[ii][kk]);
        // XXX:
        if(denom == T(0)) continue;
        const auto& score((const_cast<const Tensor&>(corpust))[ii][jj][kk]);
        result.corpust[widx][jj][kk] += score * c_jk[jj][kk] / denom;
        result.corpust[ii][widx][kk] += score * c_ik[ii][kk] / denom;
        result.corpust[ii][jj][widx] += score * c_ij[ii][jj] / denom;
        result.corpust[ii][jj][kk]   -= score;
      }
    }
  }
  return result.simpleThresh(thresh);
}

template <typename T, typename U> void corpus<T,U>::merge5(Tensor& d, const int& i, const int& ki, const int& kk, const int& kj, const int& j, const T& intensity) const {
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


template <typename T, typename U> static inline U getCut(const U& input, const int& idx, const int& szwindow) {
  if(idx < 0 || input.size() <= idx * szwindow / 2)
    return U("XXX: getCut");
  return input.substr(idx * szwindow / 2, szwindow);
}

template <typename T, typename U> static inline std::ostream& outTagged(std::ostream& os, const U& name, const corpus<T, U>& cstat, const int& idx, const T& score, const U& input, const int& szwindow) {
  os << "<span id=\"" << name << idx << "\">";
  os << "score: " << score << " : ";
  os << getCut<T, U>(input, idx, szwindow);
  os << "<input class=\"gatherdetail\" type=\"checkbox\"><div class=\"gatherdetail\">";
  os << cstat.reverseLink().second << "<br/>";
  os << cstat.getAttributed(words);
  os << "</div></span>";
  return os;
}

template <typename T, typename U> static inline void getAbbreved(corpus<T, U>& cstat, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& delimiter) {
  assert(detailtitle.size() == detail.size());
  for(int i = 0; i < detail.size(); i ++)
    cstat = cstat.abbrev(detailtitle[i], corpus<T, U>(detail[i], delimiter));
  return;
}

template <typename T, typename U> static inline bool getDetailed(corpus<T, U>& cstat, const U& input, const int& idx, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& delimiter, const int& szwindow, const T& thresh) {
  assert(detailtitle.size() == detail.size());
  if(idx < 0 || input.size() <= idx * szwindow / 2)
    return false;
  assert(0 <= idx && idx * szwindow / 2 < input.size());
  cstat = corpus<T, U>(getCut<T>(input, idx, szwindow), delimiter).simpleThresh(thresh);
  for(int i = 0; i < detail.size(); i ++)
    cstat = cstat.withDetail(detailtitle[i], corpus<T, U>(detail[i], delimiter), thresh).simpleThresh(thresh);
  return true;
}

template <typename T, typename U> std::ostream& preparedTOC(std::ostream& os, const U& input, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& topictitle, const vector<U>& topics, const vector<U>& delimiter, const int& szwindow, const int& outblock, const int& nrwords, const T& redig = T(1), const bool& reverse = false) {
  assert(detailtitle.size() == detail.size());
  assert(topictitle.size()  == topics.size());
  os << "prepTOC: " << flush;
  if(!topics.size())
    return os << "zero input. <br />";
  vector<corpus<T, U> > istats;
  T threshin(int(0));
  for(int i = - int(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) );
          i <= 0; i ++) {
    threshin = T(int(1)) - pow(T(int(2)), - T(abs(i)));
    vector<int> idx;
    istats.resize(0);
    istats.resize(input.size() / (szwindow / 2));
    for(int j = 0; j < istats.size(); j ++) {
      getDetailed<T, U>(istats[j], input, j, detailtitle, detail, delimiter, szwindow, threshin);
      istats[j].reDig(redig);
      istats[j].absfy();
      auto lidx(istats[j].countIdx());
      idx.insert(idx.end(), lidx.begin(), lidx.end());
    }
    sort(idx.begin(), idx.end());
    idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
    cerr << threshin << " : " << idx.size() << endl;
    if(nrwords <= idx.size()) break;
  }
  for(int i = 0; i < topics.size(); i ++) {
    vector<pair<T, pair<int, int> > > topicidx;
    vector<corpus<T, U> > stats;
    stats.resize(input.size() / (szwindow / 2));
    for(int j = 0; j < stats.size(); j ++) {
      getDetailed<T, U>(stats[j], topics[i], j, detailtitle, detail, delimiter, szwindow, threshin);
      stats[j].reDig(redig);
      stats[j].absfy();
    }
    for(int j = 0; j < istats.size(); j ++) {
      int idx(0);
      T   score(0);
      for(int k = 0; k < stats.size(); k ++) {
        const auto  lscore(reverse ? T(1) / abs(istats[j].prej(stats[k]))
                                   :            istats[j].prej(stats[k]) );
        if(isfinite(lscore) && score <= lscore) {
          idx   = k;
          score = lscore;
        }
      }
      topicidx.emplace_back(make_pair(- score, make_pair(j, idx)));
    }
    sort(topicidx.begin(), topicidx.end());
    for(int j = 0; j < topicidx.size(); j ++) {
      const auto& stat0(istats[topicidx[j].second.first]);
      const auto& stat1(stats[topicidx[j].second.second]);
      if(outblock < j)
        break;
      os << topictitle[i] << " ";
      outTagged<T,U>(os, U("prepTOC0_") + to_string(i) + U("_"), stat0, topicidx[j].second.first, topicidx[j].first, input, szwindow);
      os << "<br/>";
      outTagged<T,U>(os, U("prepTOC_")  + to_string(i) + U("_"), stat1, topicidx[j].second.second, topicidx[j].first, topics[i], szwindow);
      os << "<br/><br/>" << endl;
    }
    os << "<br/><br/>" << endl;
  }
  return os;
}

template <typename T, typename U> std::ostream& optimizeTOC(std::ostream& os, const U& input, const vector<U>& detail, const vector<U>& detailtitle, const vector<U>& delimiter, const int& szwindow, const int& outblock, const int& nrwords, const T& redig = T(1), const bool& countnum = false, const U& notcheck = U("")) {
  assert(notcheck == U(""));
  os << "optTOC: " << flush;
  SimpleSparseMatrix<T> scores;
  int Midx(0);
  vector<corpus<T, U> > stats;
  T   threshin(int(0));
  for(int i = - int(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) );
          i <= 0; i ++) {
    threshin = T(int(1)) - pow(T(int(2)), - T(abs(i)));
    vector<int> idx;
    stats.resize(0);
    stats.resize(input.size() / (szwindow / 2));
    for(int j = 0; j < stats.size(); j ++) {
      getDetailed<T, U>(stats[j], input, j, detailtitle, detail, delimiter, szwindow, threshin);
      stats[j].reDig(redig);
      stats[j].absfy();
      auto lidx(stats[j].countIdx());
      idx.insert(idx.end(), lidx.begin(), lidx.end());
      // XXX: ordinary C compiler freezes here with exhaust of memory, then,
      //      abort trap exit.
    }
    sort(idx.begin(), idx.end());
    idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
    cerr << threshin << " : " << idx.size() << endl;
    if(nrwords <= idx.size()) break;
  }
  for(int i = 0; i < stats.size(); i ++)
    for(int j = i + 1; j < stats.size(); j ++) {
      scores[i][j] = - stats[i].prej(stats[j]);
      Midx = max(Midx, j);
    }
  os << "." << flush;
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
          if(!binary_search(phrases.begin(), phrases.end(), j) &&
             llscore.size() <= outblock)
            llscore.emplace_back(make_pair(scores[min(i, j)][max(i, j)], make_pair(i, j)));
        sort(llscore.begin(), llscore.end());
        lscore.emplace_back(llscore);
        bool ok(false);
        if(countnum) {
          T lllscore(0);
          for(int j = 0; j < llscore.size(); j ++)
            lllscore += T(1);
          ok = Mscore < lllscore;
          Mscore = max(Mscore, lllscore);
        } else {
          T lllscore(0);
          for(int j = 0; j < llscore.size(); j ++)
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
        lscore.emplace_back(vector<pair<T, pair<int, int> > >());
    if(lscore[lidx].size() <= 0) {
      for(int i = 0; i < lscore.size(); i ++)
        if(lscore[i].size()) {
          lidx = i;
          break;
        }
      if(! lscore[lidx].size())
        break;
    }
    phrases.emplace_back(lscore[lidx][0].second.first);
    sort(phrases.begin(), phrases.end());
    const auto& cs(stats[lscore[lidx][0].second.first]);
    for(int i = 0; i < lscore[lidx].size(); i ++)
      phrases.emplace_back(lscore[lidx][i].second.second);
    sort(phrases.begin(), phrases.end());
    os << "<form action=\"../../../../puts.php\" method=\"POST\"><div>";
    os << "base : ";
    outTagged<T,U>(os, U("optTOC0_"), cs, lscore[lidx][0].second.first, lscore[lidx][0].first, input, szwindow) << "<br/>";
    os << "<br/>Show/Hide : <input class=\"gather\" type=\"checkbox\"><div class=\"gather\">";
    for(int i = 0; i < lscore[lidx].size(); i ++)
      outTagged<T,U>(os, U("optTOC_"), stats[lscore[lidx][i].second.second], lscore[lidx][i].second.second, lscore[lidx][i].first, input, szwindow) << "<br/>";
    os << "</div></div><textarea name=\"entry\">";
    os << getCut<T,U>(input, lscore[lidx][0].second.first, szwindow);
    os << "\n";
    for(int i = 0; i < lscore[lidx].size(); i ++)
      os << getCut<T,U>(input, lscore[lidx][i].second.second, szwindow) << "\n";
    os << "</textarea>";
    os << "<input type=\"hidden\" name=\"name\" value=\"append\" />";
    os << "<input type=\"hidden\" name=\"adddict\" value=\"\" />";
    os << "<input type=\"submit\" value=\"Append\" />";
    os << "</form><br/>" << endl;
  }
  return os << endl;
}

template <typename T, typename U> std::ostream& diff(std::ostream& os, const U& input, const vector<U>& detail0, const vector<U>& detailtitle0, const vector<U>& detail1, const vector<U>& detailtitle1, const vector<U>& delimiter, const int& szwindow, const int& outblock, const int& nrwords, const T& redig = T(1), const bool& same = false) {
  assert(detail0.size() == detailtitle0.size() &&
         detail1.size() == detailtitle1.size());
  os << "diff:" << flush;
  T threshin(int(0));
  for(int i = - int(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) );
          i <= 0; i ++) {
    threshin = T(int(1)) - pow(T(int(2)), - T(abs(i)));
    corpus<T, U> stat;
    vector<int> idx;
    for(int j = 0; j < input.size() / (szwindow / 2); j ++) {
      getDetailed<T, U>(stat, input, j, detailtitle0, detail0, delimiter, szwindow, threshin);
      stat.reDig(redig);
      stat.absfy();
      auto lidx(stat.countIdx());
      idx.insert(idx.end(), lidx.begin(), lidx.end());
    }
    for(int j = 0; j < input.size() / (szwindow / 2); j ++) {
      getDetailed<T, U>(stat, input, j, detailtitle1, detail1, delimiter, szwindow, threshin);
      stat.reDig(redig);
      stat.absfy();
      auto lidx(stat.countIdx());
      idx.insert(idx.end(), lidx.begin(), lidx.end());
    }
    sort(idx.begin(), idx.end());
    idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
    cerr << threshin << " : " << idx.size() << endl;
    if(nrwords <= idx.size()) break;
  }
  corpus<T, U> cstat, dstat;
  vector<pair<T, int> > scores;
  for(int i = 0; ; i ++) {
    if(!getDetailed<T, U>(cstat, input, i, detailtitle0, detail0, delimiter, szwindow, threshin) ||
       !getDetailed<T, U>(dstat, input, i, detailtitle1, detail1, delimiter, szwindow, threshin))
      break;
    getAbbreved<T, U>(cstat, detailtitle1, detail1, delimiter);
    getAbbreved<T, U>(dstat, detailtitle0, detail0, delimiter);
    cstat.reDig(redig);
    dstat.reDig(redig);
    const auto score(cstat.prej(dstat));
    os << score << ":" << flush;
    if(isfinite(score))
      scores.emplace_back(make_pair(same ? - score : score, i));
  }
  sort(scores.begin(), scores.end());
  for(int ii = 0; ii < scores.size() && ii <= outblock; ii ++) {
    const T&   score(scores[ii].first);
    const int& i(scores[ii].second);
    if(!getDetailed<T, U>(cstat, input, i, detailtitle0, detail0, delimiter, szwindow, threshin) ||
       !getDetailed<T, U>(dstat, input, i, detailtitle1, detail1, delimiter, szwindow, threshin))
      break;
    getAbbreved<T, U>(cstat, detailtitle1, detail1, delimiter);
    getAbbreved<T, U>(dstat, detailtitle0, detail0, delimiter);
    cstat.reDig(redig);
    dstat.reDig(redig);
    auto diff(cstat - dstat);
    diff.reDig(redig);
    os << "score: " << score << " : " << diff.serialize() << "<br/>" << endl;
    outTagged<T,U>(os, U("optTOC_src"), cstat, i, score, input, szwindow) << "<br/>" << endl;
    outTagged<T,U>(os, U("optTOC_dst"), dstat, i, score, input, szwindow) << "<br/>" << endl;
    outTagged<T,U>(os, U("optTOC_diff"), diff, i, score, input, szwindow) << "<br/><br/>" << endl;
  }
  return os << endl;
}

template <typename T> static inline vector<T> cutText(const T& input, const vector<T>& eliminate, const vector<T>& delimiter, const bool& f_sort = false) {
  vector<T> result;
  T         workbuf;
  for(int i = 0; i < input.size(); i ++) {
    workbuf += input[i];
    for(int j = 0; j < delimiter.size(); j ++)
      if(workbuf.size() >= delimiter[j].size() &&
         workbuf.substr(workbuf.size() - delimiter[j].size(), delimiter[j].size()) == delimiter[j]) {
        if(workbuf.size() - delimiter[j].size())
          result.emplace_back(workbuf.substr(0, workbuf.size() - delimiter[j].size()));
        workbuf = T();
        goto next;
      }
    for(int j = 0; j < eliminate.size(); j ++)
      if(workbuf.size() >= eliminate[j].size() &&
        workbuf.substr(workbuf.size() - eliminate[j].size(), eliminate[j].size()) == eliminate[j]) {
        workbuf = workbuf.substr(0, workbuf.size() - eliminate[j].size());
        break;
      }
   next:
    ;
  }
  if(workbuf.size())
    result.emplace_back(workbuf);
  if(f_sort)
    sort(result.begin(), result.end());
  return result;
}

template <typename T, typename U> static inline SimpleVector<T> countWords(const U& orig, const vector<U>& words) {
  SimpleVector<T> result(words.size());
  for(int i = 0; i < result.size(); i ++)
    result[i] = T(0);
  for(int i = 0; i < orig.size(); i ++) {
    const U work(orig.substr(i, orig.size() - i - 1));
    for(int j = 0; j < words.size(); j ++)
      if(equalStrClip<U>(work, words[j]) && work.size() >= words[j].size())
        result[j] = T(1);
  }
  return result;
}

template <typename T, typename U> vector<int> pseudoWordsBalance(const vector<U>& orig, const vector<U>& words, int nloop = - 1) {
  cerr << "pseudoWordsBalance : " << orig.size() << ", " << words.size() << flush;
  SimpleMatrix<T> result(words.size(), orig.size());
  for(int i = 0; i < orig.size(); i ++)
    result.setCol(i, countWords<T, U>(orig[i], words));
  
  vector<int> vres;
  if(nloop <= 0)
    nloop = result.cols();
  
  for(int i = 0; i < min(int(result.cols()), nloop); i ++) {
    vector<pair<T, int> > scores;
    for(int j = 0; j < result.cols(); j ++)
      scores.emplace_back(make_pair(- result.col(j).dot(result.col(j)), j));
    sort(scores.begin(), scores.end());
    if(scores[0].first == T(0))
      break;
    const int& idx(scores[0].second);
    vres.emplace_back(idx);
    for(int k = 0; k < result.rows(); k ++)
      if(result(k, idx) != T(0))
        result.row(k) *= T(0);
    T sum(0);
    for(int j = 0; j < result.rows(); j ++)
      sum += result.row(j).dot(result.row(j));
    cerr << ":" << sum << flush;
  }
  cerr << endl;
  return vres;
}

template <typename T, typename U> std::ostream& predTOC(std::ostream& os, const U& input, const vector<U>& detailtitle, const vector<U>& detail, const vector<U>& delimiter, const int& szwindow, const int& nrwords0, const T& redig = T(1) ) {
  assert(detailtitle.size() == detail.size());
  const int nrwords(sqrt(T(nrwords0)));
  os << "predTOC: " << flush;
  vector<SimpleSparseTensor<T> > in;
  T threshin(int(0));
  vector<int> idx;
  for(int i = - int(- log(SimpleMatrix<T>().epsilon()) / log(T(int(2))) ) * 2;
          i <= 0; i ++) {
    vector<corpus<T, U> > istats;
    threshin = T(int(1)) - pow(T(int(2)), - T(abs(i)));
    idx.resize(0);
    istats.resize(input.size() / (szwindow / 2));
    for(int j = 0; j < istats.size(); j ++) {
      getDetailed<T, U>(istats[j], input, j, detailtitle, detail, delimiter, szwindow, threshin);
      istats[j].reDig(redig);
      istats[j].absfy();
      auto lidx(istats[j].countIdx());
      idx.insert(idx.end(), lidx.begin(), lidx.end());
    }
    sort(idx.begin(), idx.end());
    idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
    cerr << threshin << " : " << idx.size() << endl;
    if(nrwords <= idx.size()) {
      in.resize(istats.size());
      for(int j = 0; j < istats.size(); j ++)
        in[j] = move(istats[j].corpust);
      break;
    }
  }
  auto p(predSTen<T>(in, idx));
  p.first.insert(p.first.end(), p.second.begin(), p.second.end());
  vector<string> hist;
  hist.reserve(p.first.size());
  for(int i = 0; i < p.first.size(); i ++) {
    corpus<T, U> pstats;
    pstats.corpust = p.first[i];
    getAbbreved<T>(pstats, detailtitle, detail, delimiter);
    auto serial(pstats.simpleThresh(threshin).serialize());
    if(binary_search(hist.begin(), hist.end(), serial)) continue;
    os << serial << "<br /><br />" << endl;
    hist.emplace_back(move(serial));
    sort(hist.begin(), hist.end());
    if(i == p.first.size() / 2 - 1) os << "<hr />" << endl;
  }
  return os;
}

#define _CORPUS_
#endif


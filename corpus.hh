/* BSD 3-Clause License:
 * Copyright (c) 2018-2024, bitsofcotton.
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

// N.B. will frequently update when NOT word table is in adjust.

using std::move;
using std::cerr;
using std::endl;
using std::flush;
using std::string;
using std::to_string;
using std::lower_bound;
using std::upper_bound;
using std::binary_search;
using std::pair;
using std::make_pair;
using std::vector;
using std::unique;
using std::distance;
using std::sort;
using std::isfinite;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::max;
using std::min;

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
    assert(T(int(0)) <= ratio);
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
  int i(0), i0(0), Midx(0), lastlen(0);
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
          lastlen = max(lastlen, int(work.size()));
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
      i0 = (i -= work.size() - lastlen - 1) + 1;
      lastlen = 0;
    } else
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
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
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
          if(isfinite(work)) {
#if defined(_OPENMP)
#pragma omp critical
#endif
            {
              corpust[i][j][k] += sqrt(work) / T(int(Midx));
            }
          }
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

#define _CORPUS_
#endif


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
#if !defined(_LWORD_)

#include "corpus.hh"

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

template <typename T> class word_t {
public:
  T   str;
  int count;
  word_t() {
    str   = T();
    count = 0;
  }
  ~word_t() {
  }
  word_t(const word_t<T>& other) {
    *this = other;
  }
  word_t& operator = (const word_t<T>& other) {
    str    = other.str;
    count  = other.count;
    return *this;
  }
  bool operator < (const word_t<T>& x1) const {
    return (count < x1.count) || (count == x1.count && str < x1.str);
  }
};

template <typename T> class gram_t {
public:
  T           str;
  vector<int> ptr0;
  vector<int> ptr1;
  gram_t() {
    this->str  = T();
    this->ptr0 = vector<int>();
    this->ptr1 = vector<int>();
  }
  ~gram_t() {
    ;
  }
  gram_t(const gram_t<T>& x) {
    *this = x;
  }
  gram_t& operator = (const gram_t<T>& x) {
    str  = x.str;
    ptr0 = x.ptr0;
    ptr1 = x.ptr1;
    return *this;
  }
  bool operator < (const gram_t<T>& x1) const {
    return str < x1.str;
  }
};

template <typename T, typename U> class lword {
public:
  lword();
  ~lword();
  void init(const int& loop, const int& mthresh, const int& Mthresh);
  
  const vector<word_t<U> >& compute(const U& input);
private:
  vector<T>                   dict0;
  vector<vector<gram_t<U> > > dicts;
  vector<word_t<U> >          words;
  
  int mthresh;
  int Mthresh;
  
  void makeBigram(const U& input);
  void constructNwords();
  
  bool       isin(const U& key);
  gram_t<U>& find(const U& key);
  void       assign(const gram_t<U>& val);
};

template <typename T, typename U> lword<T, U>::lword() {
  init(60, 2, 2);
}

template <typename T, typename U> lword<T, U>::~lword() {
  // already freed in another places.
  ;
}

template <typename T, typename U> void lword<T, U>::init(const int& loop, const int& mthresh, const int& Mthresh) {
  this->mthresh = mthresh;
  this->Mthresh = Mthresh;
/*
  this->dict0   = vector<T>();
  this->dicts   = vector<vector<gram_t<U> > >();
*/
  this->dicts.resize(loop, vector<gram_t<U> >());
//  this->words   = vector<word_t<U> >();
  return;
}

template <typename T, typename U> bool lword<T, U>::isin(const U& key) {
  assert(key.size() < dicts.size());
  const vector<gram_t<U> >& dict(dicts[key.size()]);
  gram_t<U> key0;
  key0.str = key;
  auto p(lower_bound(dict.begin(), dict.end(), key0));
  return dict.begin() <= p && p < dict.end() && p->str == key;
}

template <typename T, typename U> gram_t<U>& lword<T, U>::find(const U& key) {
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

template <typename T, typename U> void lword<T, U>::assign(const gram_t<U>& val) {
  assert(val.str.size() < dicts.size());
  vector<gram_t<U> >& dict(dicts[val.str.size()]);
  auto p(lower_bound(dict.begin(), dict.end(), val));
  if(val.ptr0.size()) {
    // delete duplicates:
    gram_t<U> work;
    work.str = val.str;
    // XXX can spped up with sorting.
    for(int i = 0; i < val.ptr0.size(); i ++) {
      bool flag(false);
      for(int j = 0; j < work.ptr0.size(); j ++)
        if(work.ptr0[j] == val.ptr0[i] &&
           work.ptr1[j] == val.ptr1[i]) {
          flag = true;
          break;
        }
      if(flag)
        continue;
      work.ptr0.push_back(val.ptr0[i]);
      work.ptr1.push_back(val.ptr1[i]);
    }
    sort(work.ptr0.begin(), work.ptr0.end());
    sort(work.ptr1.begin(), work.ptr1.end());
    if(p < dict.begin() || dict.end() <= p || p->str != work.str) {
      dict.push_back(work);
      sort(dict.begin(), dict.end());
    } else
      *p = work;
  } else if(dict.begin() <= p && p < dict.end() && p->str == val.str)
    dict.erase(p);
  return;
}

template <typename T, typename U> const vector<word_t<U> >& lword<T, U>::compute(const U& input) {
  makeBigram(input);
  constructNwords();
  for(auto itr = dicts.begin(); itr != dicts.end(); ++ itr)
    for(auto itr2 = itr->begin(); itr2 != itr->end(); ++ itr2) {
      word_t<U> work;
      work.count = itr2->ptr0.size();
      if(!work.count)
        continue;
      work.str   = itr2->str;
      bool flag(false);
      for(int i = 0; i < words.size(); i ++)
        if(work.str == words[i].str) {
          flag = true;
          break;
        }
      if(flag) continue;
      words.push_back(work);
    }
  sort(words.begin(), words.end());
  return words;
}

template <typename T, typename U> void lword<T, U>::makeBigram(const U& input) {
  map<U, vector<int> > mapw;
  for(int i = 1; i < input.size(); i ++) {
    U work;
    work += T(input[i - 1]);
    work += T(input[i]);
    mapw[work].push_back(i);
    if(!binary_search(dict0.begin(), dict0.end(), input[i])) {
      dict0.push_back(input[i]);
      sort(dict0.begin(), dict0.end());
    }
  }
  for(auto itr = mapw.begin(); itr != mapw.end(); ++ itr) {
    gram_t<U> work;
    work.str = itr->first;
    for(auto itr2 = itr->second.begin(); itr2 != itr->second.end(); ++ itr2) {
      work.ptr0.push_back(*itr2 - 1);
      work.ptr1.push_back(*itr2);
    }
    assign(work);
  }
  return;
}

template <typename T, typename U> void lword<T, U>::constructNwords() {
  int i;
  for(i = 2; i < dicts.size(); i ++) {
    cerr << i << flush;
    map<U, vector<int> > map0;
    map<U, vector<int> > map1;
    map<U, vector<int> > dmap;
    for(auto itr = dicts[i].begin(); itr != dicts[i].end(); ++ itr) {
      cerr << "." << flush;
      const gram_t<U>& idxkey(*itr);
      for(auto itr2 = dict0.begin(); itr2 != dict0.end(); ++ itr2) {
        U key2;
        for(int j = 1; j < idxkey.str.size(); j ++)
          key2 += idxkey.str[j];
        key2 += T(*itr2);
        if(!isin(key2))
          continue;
        const gram_t<U>& idxkey2(find(key2));
        U workkey(idxkey.str);
        workkey += *itr2;
        if(isin(workkey))
          continue;
        vector<int> idxwork[4];
        for(int j = 0; j < 4; j ++)
          idxwork[j] = vector<int>();
        int tt = 0;
        int ss = 0;
        while(tt < idxkey2.ptr0.size() && ss < idxkey.ptr0.size()) {
          if(idxkey2.ptr1[tt] == idxkey.ptr1[ss] + 1) {
            idxwork[0].push_back(idxkey.ptr0[ss]);
            idxwork[1].push_back(idxkey2.ptr1[tt]);
            idxwork[2].push_back(ss);
            idxwork[3].push_back(tt);
            ss ++;
            tt ++;
          } else if(idxkey2.ptr1[tt] > idxkey.ptr1[ss])
            ss ++;
          else
            tt ++;
        }
        if(idxwork[0].size() < Mthresh)
          continue;
        const int diff(min( idxkey.ptr0.size(), idxkey2.ptr0.size()) - idxwork[0].size());
        const int diffM(max(idxkey.ptr0.size(), idxkey2.ptr0.size()) - idxwork[0].size());
        if(diff < 0 || mthresh < diffM)
          continue;
        map0[workkey].insert(map0[workkey].end(), idxwork[0].begin(), idxwork[0].end());
        map1[workkey].insert(map1[workkey].end(), idxwork[1].begin(), idxwork[1].end());
        dmap[idxkey.str].insert(dmap[idxkey.str].end(), idxwork[2].begin(), idxwork[2].end());
        dmap[key2].insert(dmap[key2].end(), idxwork[3].begin(), idxwork[3].end());
      }
    }
    // delete this stage garbage.
    for(auto itr = dmap.begin(); itr != dmap.end(); ++ itr) {
      sort(itr->second.begin(), itr->second.end());
      const gram_t<U>& before(find(itr->first));
            gram_t<U>  after(before);
      after.ptr0 = vector<int>();
      after.ptr1 = vector<int>();
      for(int j = 0; j < before.ptr1.size(); j ++)
        if(!binary_search(itr->second.begin(), itr->second.end(), j)) {
          after.ptr0.push_back(before.ptr0[j]);
          after.ptr1.push_back(before.ptr1[j]);
        }
      assign(after);
    }
    // construct next stage.
    for(auto itr = map0.begin(); itr != map0.end(); ++ itr) {
      gram_t<U> work;
      work.str  = itr->first;
      work.ptr0 = itr->second;
      work.ptr1 = map1[itr->first];
      assign(work);
    }
    cerr << endl;
  }
  return;
}


// some functions.
template <typename T> vector<T> cutText(const T& input, const vector<T>& eliminate0, const vector<T>& delimiter0, const bool& f_sort = false) {
  vector<T> delimiter(delimiter0);
  vector<T> eliminate(eliminate0);
  sort(delimiter.begin(), delimiter.end());
  sort(eliminate.begin(), eliminate.end());
  delimiter.erase(unique(delimiter.begin(), delimiter.end()), delimiter.end());
  eliminate.erase(unique(eliminate.begin(), eliminate.end()), eliminate.end());
  vector<T> result;
  T         workbuf;
  for(int i = 0; i < input.size(); i ++) {
    workbuf += input[i];
    for(int j = 0; j < delimiter.size(); j ++)
      if(workbuf.size() >= delimiter[j].size() &&
         workbuf.substr(workbuf.size() - delimiter[j].size(), delimiter[j].size()) == delimiter[j]) {
        if(workbuf.size() - delimiter[j].size())
          result.push_back(workbuf.substr(0, workbuf.size() - delimiter[j].size()));
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
    result.push_back(workbuf);
  if(f_sort)
    sort(result.begin(), result.end());
  return result;
}

template <typename T, typename U> Eigen::Matrix<T, Eigen::Dynamic, 1> countWords(const U& orig, const vector<U>& words) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> result(words.size());
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

template <typename T, typename U> vector<int> pseudoWordsBalance(const vector<U>& orig, const vector<U>& words0, int nloop = - 1) {
  cerr << "pseudoWordsBalance : " << orig.size() << ", " << words0.size() << flush;
  vector<U> words(words0);
  sort(words.begin(), words.end());
  words.erase(unique(words.begin(), words.end()), words.end());
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(words.size(), orig.size());
  for(int i = 0; i < orig.size(); i ++)
    result.col(i) = countWords<T, U>(orig[i], words);
  
  vector<int> vres;
  if(nloop <= 0)
    nloop = result.cols();
  
  for(int i = 0; i < min(int(result.cols()), nloop); i ++) {
    vector<pair<T, int> > scores;
    for(int j = 0; j < result.cols(); j ++)
      scores.push_back(make_pair(- result.col(j).dot(result.col(j)), j));
    sort(scores.begin(), scores.end());
    if(scores[0].first == T(0))
      break;
    const int& idx(scores[0].second);
    vres.push_back(idx);
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

#define _LWORD_
#endif


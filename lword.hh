#if !defined(_LWORD_)

#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <cstring>

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

template <typename T> vector<T> cutText(const T& input, const vector<T>& eliminate0, const vector<T>& delimiter0);

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
  
  const vector<word_t<U> >& compute(const T* input);
private:
  vector<T>                   dict0;
  vector<vector<gram_t<U> > > dicts;
  vector<word_t<U> >          words;
  
  int mthresh;
  int Mthresh;
  
  void makeBigram(const T* input);
  void constructNwords();
  void makeWords();
  
  bool       isin(const U key);
  gram_t<U>& find(const U key);
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
  this->dict0   = vector<T>();
  this->dicts   = vector<vector<gram_t<U> > >();
  for(int i = 0; i < loop; i ++)
    this->dicts.push_back(vector<gram_t<U> >());
  this->words   = vector<word_t<U> >();
  return;
}

template <typename T, typename U> bool lword<T, U>::isin(const U key) {
  if(dicts.size() < key.size()) {
    cerr << "dictionary small" << endl;
    return false;
  }
  const vector<gram_t<U> >& dict(dicts[key.size()]);
  gram_t<U> key0;
  key0.str = key;
  auto p(lower_bound(dict.begin(), dict.end(), key0));
  return !(p < dict.begin() || dict.end() <= p) && p->str == key;
}

template <typename T, typename U> gram_t<U>& lword<T, U>::find(const U key) {
  if(dicts.size() < key.size()) {
    cerr << "XXX: find slips." << endl;
    static gram_t<U> dummy;
    return dummy;
  }
  vector<gram_t<U> >& dict(dicts[key.size()]);
  gram_t<U> key0;
  key0.str = key;
  auto p(lower_bound(dict.begin(), dict.end(), key0));
  if(p < dict.begin() || dict.end() <= p || p->str != key) {
    cerr << "XXX: slipping find." << endl;
    static gram_t<U> dummy;
    return dummy;
  }
  return *p;
}

template <typename T, typename U> void lword<T, U>::assign(const gram_t<U>& val) {
  if(dicts.size() < val.str.size()) {
    cerr << "assign : dictionary small." << endl;
    return;
  }
  vector<gram_t<U> >& dict(dicts[val.str.size()]);
  auto p(lower_bound(dict.begin(), dict.end(), val));
  if(val.ptr0.size()) {
    // delete duplicates:
    gram_t<U> work;
    work.str = val.str;
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

template <typename T, typename U> const vector<word_t<U> >& lword<T, U>::compute(const T* input) {
  cerr << "bi-gramming..." << endl;
  makeBigram(input);
  cerr << "constructing longest words..." << endl;
  constructNwords();
  cerr << "making word tables..." << endl;
  makeWords();
  return words;
}

template <typename T, typename U> void lword<T, U>::makeBigram(const T* input) {
  map<T, T> map0;
  map<U, vector<int> > map;
  for(int i = 1; input[i]; i ++) {
    U work;
    work += T(input[i - 1]);
    work += T(input[i]);
    map[work].push_back(i);
    map0[input[i]] = input[i];
  }
  for(auto itr = map0.begin(); itr != map0.end(); ++ itr)
    dict0.push_back(itr->first);
  sort(dict0.begin(), dict0.end());
  for(auto itr = map.begin(); itr != map.end(); ++ itr) {
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
  for(i = 0; i < dicts.size(); i ++) {
    cerr << "constructing " << i + 2 << " table";
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
        workkey += T(*itr2);
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
        if(0 < diff && mthresh < diffM)
          continue;
        for(int i = 0; i < idxwork[0].size(); i ++) {
          map0[U(workkey)].push_back(idxwork[0][i]);
          map1[U(workkey)].push_back(idxwork[1][i]);
          dmap[U(idxkey.str)].push_back(idxwork[2][i]);
          dmap[U(key2)].push_back(idxwork[3][i]);
        }
      }
    }
    // delete this stage garbage.
    for(auto itr = dmap.begin(); itr != dmap.end(); ++ itr) {
      sort(itr->second.begin(), itr->second.end());
      const gram_t<U>& before(find(itr->first));
            gram_t<U>  after(before);
      after.ptr0 = vector<int>();
      after.ptr1 = vector<int>();
      for(int j = 0; j < before.ptr1.size(); j ++) {
        bool flag = false;
        for(auto itr2 = itr->second.begin(); itr2 != itr->second.end(); ++ itr2)
          if(*itr2 == j) {
            flag = true;
            break;
          }
        if(flag)
          continue;
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

template <typename T, typename U> void lword<T, U>::makeWords() {
  for(auto itr = dicts.begin(); itr != dicts.end(); ++ itr) {
    for(auto itr2 = itr->begin(); itr2 != itr->end(); ++ itr2) {
      word_t<U> work;
      work.count = itr2->ptr0.size();
      if(!work.count)
        continue;
      work.str   = itr2->str;
      // XXX fixme not here.
      bool flag(false);
      for(int i = 0; i < words.size(); i ++)
        if(work.str == words[i].str) {
          flag = true;
          break;
        }
      if(flag) continue;
      words.push_back(work);
    }
  }
  sort(words.begin(), words.end());
  return;
}


// some functions.
template <typename T> vector<T> cutText(const T& input, const vector<T>& eliminate0, const vector<T>& delimiter0) {
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
    bool flag(false);
    for(int j = 0; j < delimiter.size(); j ++)
      if(workbuf.size() >= delimiter[j].size() &&
         workbuf.substr(workbuf.size() - delimiter[j].size(), delimiter[j].size()) == delimiter[j]) {
        if(workbuf.size() - delimiter[j].size()) {
          result.push_back(workbuf.substr(0, workbuf.size() - delimiter[j].size()));
          sort(result.begin(), result.end());
          result.erase(unique(result.begin(), result.end()), result.end());
        }
        workbuf = T();
        flag = true;
        break;
      }
    if(flag) continue;
    for(int j = 0; j < eliminate.size(); j ++)
      if(workbuf.size() >= eliminate[j].size() &&
        workbuf.substr(workbuf.size() - eliminate[j].size(), eliminate[j].size()) == eliminate[j]) {
        workbuf = workbuf.substr(0, workbuf.size() - eliminate[j].size());
        break;
      }
  }
  if(workbuf.size())
    result.push_back(workbuf);
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

// XXX very pseudo, this cannot be applied with simple DP.
template <typename T, typename U> vector<int> pseudoWordsBalance(const vector<U>& orig, const vector<U>& words0, int nloop = - 1) {
  cerr << "pseudoWordsBalance initializing matrix with initial words " << orig.size() << ", " << words0.size() << flush;
  vector<U> words(words0);
  sort(words.begin(), words.end());
  words.erase(unique(words.begin(), words.end()), words.end());
  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(words.size(), orig.size());
  for(int i = 0; i < orig.size(); i ++)
    result.col(i) = countWords<T, U>(orig[i], words);
  
  vector<int> vres;
  if(nloop <= 0)
    nloop = result.cols();
  
  cerr << "pseudo-balancing..." << flush;
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
    cerr << sum << ":" << flush;
  }
  return vres;
}

#define _LWORD_
#endif


#if !defined(_CORPUS_)

#include <Eigen/Core>
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
using std::pair;
using std::make_pair;
using std::to_string;
using std::isfinite;
using std::sqrt;
using std::min;
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
  
  void init(const U* words, const int& nthresh, const int& Nthresh, const int& Mwords = 150);
  corpus<T, U>&         operator = (const corpus<T, U>& other);
  const void            compute(const U* input, const vector<string>& delimiter = vector<string>());
  const vector<string>& getWords() const;
  const Tensor&         getCorpus() const;
private:
  vector<string>       words0;
  vector<string>       words;
  vector<vector<int> > ptrs0;
  vector<vector<int> > ptrs;
  vector<int>          pdelim;
  Tensor               corpust;
  int                  nthresh;
  int                  Nthresh;
  int                  Mwords;
  T                    Midx;
   
  void getWordPtrs(const U* input, const vector<string>& delimiter);
  void corpusEach();
};

template <typename T, typename U> corpus<T,U>::corpus() {
  ;
}

template <typename T, typename U> corpus<T,U>::~corpus() {
  // auto called destructors for string.
  ;
}

template <typename T, typename U> void corpus<T,U>::init(const U* words, const int& nthresh, const int& Nthresh, const int& Mwords) {
  string buf;
  for(int i = 0; words[i]; i ++) {
    if(words[i] == ',' || words[i] == '\n' || words[i] == '\r') {
      if(buf.size()) {
        switch(buf[buf.size() - 1]) {
        case '\n': case '\r': case ' ': case '\t': case ',':
          this->words0.push_back(buf.substr(0, buf.size() - 1));
          break;
        default:
          this->words0.push_back(buf);
        }
        this->ptrs0.push_back(vector<int>());
      }
      buf = string();
      continue;
    } else if(words[i] == ' ' || words[i] == '\t')
      continue;
    buf += words[i];
  }
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

template <typename T, typename U> const void corpus<T,U>::compute(const U* input, const vector<string>& delimiter) {
  cerr << "Getting word pointers." << endl;
  getWordPtrs(input, delimiter);
  if(Mwords < words.size()) {
    cerr << "exceeds Mwords." << endl;
    return;
  }
  cerr << "Corpus..." << endl;
  corpusEach();
  return;
}

template <typename T, typename U> const vector<string>& corpus<T,U>::getWords() const {
  return words;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic>& corpus<T,U>::getCorpus() const {
  return corpust;
}

template <typename T, typename U> void corpus<T,U>::getWordPtrs(const U* input, const vector<string>& delimiter) {
  sort(words0.begin(), words0.end());
  pdelim = vector<int>();
  pdelim.push_back(0);
  string work;
  vector<int> matchwidx;
  vector<int> matchidxs;
  int dM(0);
  for(int i = 0; i < delimiter.size(); i ++)
    dM = max(dM, int(delimiter[i].size()));
  vector<string> workd;
  for(int i = 0; i < dM; i ++) {
    workd.push_back(string(""));
    for(int j = i; j < dM; j ++)
      workd[i] += string(" ");
  }
  int i(0);
  for( ; input[i]; i ++) {
    work += input[i];
    for(int ii = 0; ii < workd.size(); ii ++) {
      workd[ii]  = workd[ii].substr(1, workd[ii].size() - 1);
      workd[ii] += input[i];
      for(int j = 0; j < delimiter.size(); j ++)
        if(workd[ii] == delimiter[j] && i != pdelim[pdelim.size() - 1])
          pdelim.push_back(i);
    }
    auto lo(words0.begin() + distance(words0.begin(), upper_bound(words0.begin(), words0.end(), work, lessEqualStrClip<string>)));
    auto up(words0.begin() + distance(words0.begin(), upper_bound(words0.begin(), words0.end(), work, lessNotEqualStrClip<string>)));
    bool match(false);
    for(auto itr(lo); itr < up; ++ itr) {
      if(equalStrClip<string>(work, *itr)) {
        if(work.size() >= itr->size()) {
          matchwidx.push_back(distance(words0.begin(), itr));
          matchidxs.push_back(i);
          match = true;
        } else if(work.size() < itr->size())
          match = true;
      }
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
    i   -= work.size() - 1;
    work = string();
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
  words = vector<string>();
  words.push_back(string("^"));
  vector<int> head, tail;
  head.push_back(0);
  ptrs.push_back(head);
  for(auto itr = words0.begin(); itr != words0.end(); ++ itr) {
    const int idx = distance(words0.begin(), itr);
    if(ptrs0[idx].size()) {
      words.push_back(*itr);
      ptrs.push_back(ptrs0[idx]);
    }
  }
  words.push_back(string("$"));
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
    if(!ptrs[i].size())
      continue;
    for(int j = 0; j < words.size(); j ++) {
      int kk(0);
      if(!ptrs[j].size())
        continue;
      corpust(i, j) = Vec(words.size());
      for(int k = 0; k < corpust(i, j).size(); k ++)
        corpust(i, j)[k] = 0;
      for(int k = 0; k < words.size(); k ++) {
        if(!ptrs[k].size())
          continue;
        int ctru = 0;
        int ctrv = 0;
        for(auto itr = ptrs[k].begin(); itr != ptrs[k].end(); ++ itr) {
          while(ctru < ptrs[i].size() && ptrs[i][ctru] < *itr) ctru ++;
          if(ctru >= ptrs[i].size())
            continue;
          while(ctrv < ptrs[j].size() && ptrs[j][ctrv] < *itr) ctrv ++;
          if(ptrs[j][ctrv] < *itr)
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
          const T work(buf0 * buf0 + buf1 * buf1);
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
  
  const corpushl<T, U>  cast(const vector<string>& words) const;
  const corpushl<T, U>& operator += (const corpushl<T, U>& other);
  const corpushl<T, U>& operator -= (const corpushl<T, U>& other);
  const corpushl<T, U>& operator *= (const T& t);
  const bool            operator <  (const corpushl<T, U>& other) const;
  const corpushl<T, U>  operator +  (const corpushl<T, U>& other) const;
  const corpushl<T, U>  operator -  () const;
  const corpushl<T, U>  operator -  (const corpushl<T, U>& other) const;
  const corpushl<T, U>  operator *  (const T& t)                  const;
  const corpushl<T, U>& operator =  (const corpushl<T, U>& other);
  const corpushl<T, U>  withDetail(const string& word, const corpushl<T, U>& other);
  const T               distanceInUnion(const corpushl<T, U>& other) const;
  const corpushl<T, U>  match2relPseudo(const corpushl<T, U>& other) const;
  const string          toc(const vector<corpushl<T, U> >& base, const vector<string>& words, const vector<corpushl<T, U> >& meanings, const string& mwords, const int& nloop, const T& thresh) const;
  const string          toc(const vector<vector<corpushl<T, U> > >& base, const vector<string>& words, const vector<corpushl<T, U> >& meanings, const string& mwords, const int& nloop, const T& thresh) const;
  const vector<T>       prej(const vector<corpushl<T, U> >& prejs) const;
  const T               culturalConflicts(const corpushl<T, U>& base) const;
  const corpushl<T, U>  conflictPart();
  const string          optToc();
  const vector<string>& getWords() const;
  const Tensor&         getCorpus() const;
  const string          serialize(const T& thresh, const T& threshall, const T& threshm = - T(1)) const;
  const string          summary(const vector<string>& words, const vector<corpushl<T, U> >& meanings, const T& thresh) const;
  const corpushl<T, U>  abbrev(const string word, const corpushl<T, U>& mean);
  const vector<string>  reverseLink(const corpushl<T, U>& orig) const;
  const corpushl<T, U>  reDig(const T& ratio);
  const corpushl<T, U>  simpleThresh(const T& ratio);
private:
  const string          serializeSub(const vector<pair<T, pair<pair<int, int>, int> > >& idxs, const bool& sign) const;
  vector<string>        words;
  Tensor                corpust;
   
  vector<string> gatherWords(const vector<string>& in0, const vector<string>& in1, vector<int>& ridx0, vector<int>& ridx1) const;
  const Tensor   prepareDetail(const corpushl<T, U>& other, const vector<string>& workwords, const int& eidx, const vector<int>& ridx0, const vector<int>& ridx1, const vector<int>& ridx2);
};

template <typename T, typename U> corpushl<T,U>::corpushl() {
  ;
}

template <typename T, typename U> corpushl<T,U>::~corpushl() {
  // auto called destructors for string.
  ;
}

template <typename T, typename U> corpushl<T,U>::corpushl(const corpus<T, U>& obj) {
  words   = vector<string>(obj.getWords());
  corpust = Tensor(obj.getCorpus());
}

template <typename T, typename U> corpushl<T,U>::corpushl(const corpushl<T, U>& obj) {
  *this = obj;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::cast(const vector<string>& words) const {
  cerr << "cast" << endl;
  corpushl<T, U> result;
  vector<string> sword(words);
  vector<int>    idxs;
  sort(sword.begin(), sword.end());
  for(int i = 0; i < sword.size(); i ++) {
    auto p(upper_bound(words.begin(), words.end(), sword[i]));
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
  words   = vector<string>(other.words);
  corpust = Tensor(other.corpust);
  return *this;
}

template <typename T, typename U> const corpushl<T, U>& corpushl<T, U>::operator += (const corpushl<T, U>& other) {
  return *this = *this + other;
}

template <typename T, typename U> const corpushl<T, U>& corpushl<T, U>::operator -= (const corpushl<T, U>& other) {
  return *this = *this - other;
}

template <typename T, typename U> const corpushl<T, U>& corpushl<T, U>::operator *= (const T& t) {
  cerr << "operator *=" << endl;
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
  cerr << "operator +" << endl;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.corpust.rows(); i ++)
    for(int j = 0; j < result.corpust.cols(); j ++) {
      result.corpust(i, j) = Vec(result.words.size());
      assert((0 <= ridx0[i] && 0 <= ridx0[j]) ||
             (0 <= ridx1[i] && 0 <= ridx1[j]));
      for(int k = 0; k < result.corpust(i, j).size(); k ++) {
        result.corpust(i, j)[k] = 0.;
        if(0 <= ridx0[i] && 0 <= ridx0[j] && 0 <= ridx0[k])
          result.corpust(i, j)[k] += corpust(ridx0[i], ridx0[j])[ridx0[k]];
        if(0 <= ridx1[i] && 0 <= ridx1[j] && 0 <= ridx1[k])
          result.corpust(i, j)[k] += other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
      }
    }
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator - () const {
  cerr << "operator -" << endl;
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

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::withDetail(const string& word, const corpushl<T, U>& other) {
  cerr << "withDetail : enter" << endl;
  if(words.size() <= 0 || other.words.size() <= 0)
    return *this;
  auto itr(upper_bound(words.begin(), words.end(), word));
  int fidx(distance(words.begin(), itr));
  if(!(0 <= fidx && fidx < words.size() && *itr == word))
    return *this;
  std::cout << *itr << ", " << word << ": " << itr->size() << " / " << word.size() << endl;
  corpushl<T, U> result;
  vector<int>    ridx0, ridx1, ridx2;
  vector<string> workwords(gatherWords(words, other.words, ridx0, ridx1));
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
  if(!flag1 || eidx < 0 || ridx1[eidx] < 0)
    return result;
  result.corpust = Tensor(workwords.size() - 1, workwords.size() - 1);
  Tensor work(prepareDetail(other, workwords, eidx, ridx0, ridx1, ridx2));
  for(int i = 0; i < workwords.size(); i ++)
    if(ridx0[i] >= 0 && ridx2[i] >= 0)
      result.words.push_back(words[ridx0[i]]);
    else if(ridx1[i] >= 0 && ridx2[i] >= 0)
      result.words.push_back(other.words[ridx1[i]]);
#if defined(_OPENMP)
#pragma omp paralell for schedule(static, 1)
#endif
  for(int i = 0; i < workwords.size(); i ++) if(ridx2[i] >= 0)
    for(int j = 0; j < workwords.size(); j ++) if(ridx2[j] >= 0) {
      result.corpust(ridx2[i], ridx2[j]) = Vec(workwords.size() - 1);
      for(int k = 0; k < result.corpust(ridx2[i], ridx2[j]).size(); k ++) if(ridx2[k] >= 0) {
        if(ridx0[i] >= 0 && ridx0[j] >= 0 && ridx0[k] >= 0)
          result.corpust(ridx2[i], ridx2[j])[ridx2[k]] = corpust(ridx0[i], ridx0[j])[ridx0[k]];
        else
          result.corpust(ridx2[i], ridx2[j])[ridx2[k]] = T(0);
      }
    }
  cerr << "withDetail: adding detail table." << endl;
  result.corpust += work;
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::match2relPseudo(const corpushl<T, U>& other) const {
  cerr << "match2relPseudo : 0/2" << endl;
  corpushl<T, U> result(*this);
  vector<int>    ridx0, ridx1;
  vector<string> words(gatherWords(result.words, other.words, ridx0, ridx1));
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
        mul(ridx0[i], ridx0[j])[ridx0[k]] = 1.;
        mul(ridx0[j], ridx0[k])[ridx0[i]] = 1.;
        mul(ridx0[k], ridx0[i])[ridx0[j]] = 1.;
      }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.corpust.rows(); i ++)
    for(int j = 0; j < result.corpust.cols(); j ++)
      for(int k = 0; k < result.corpust(i, j).size(); k ++)
        result.corpust(i, j)[k] *= mul(i, j)[k];
  return result;
}

template <typename T, typename U> const T corpushl<T, U>::distanceInUnion(const corpushl<T, U>& other) const {
  T res(0);
  vector<int>    ridx0, ridx1;
  vector<string> drop(gatherWords(words, other.words, ridx0, ridx1));
  for(int i = 0; i < drop.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0)
    for(int j = 0; j < drop.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < drop.size(); k ++) if(ridx0[k] >= 0 && ridx1[k] >= 0)
        res += corpust(ridx0[i], ridx0[j])[ridx0[k]] * other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
  return res;
}

template <typename T, typename U> const string corpushl<T, U>::toc(const vector<corpushl<T, U> >& base, const vector<string>& words, const vector<corpushl<T, U> >& meanings, const string& mwords, const int& nloop, const T& thresh) const {
  vector<vector<corpushl<T, U> > > basev;
  for(int i = 0; i < base.size(); i ++) {
    vector<corpushl<T, U> > work;
    work.push_back(base[i]);
    basev.push_back(work);
  }
  return toc(basev, words, meanings, mwords, nloop, thresh);
}

template <typename T, typename U> const string corpushl<T, U>::toc(const vector<vector<corpushl<T, U> > >& base, const vector<string>& words, const vector<corpushl<T, U> >& meanings, const string& mwords, const int& nloop, const T& thresh) const {
  string result;
  corpushl<T, U> work0(*this);
  result += mwords + string(" : ");
  vector<pair<T, string> > mlist;
  for(int j = 0; j < meanings.size(); j ++) {
    corpushl<T, U> work(this->match2relPseudo(meanings[j]));
    work0 -= work;
    for(int k = 0; k < base.size(); k ++)
      for(int l = 0; l < base.size(); l ++)
        work = work.abbrev(words[k], base[k][l]);
    mlist.push_back(make_pair(work.distanceInUnion(work) /
                              distanceInUnion(*this),
                              work.serialize(thresh, pow(thresh, T(8)))));
  }
  sort(mlist.begin(), mlist.end());
  result += to_string(mlist[mlist.size() - 1].first) + string(":") +
            mlist[mlist.size() - 1].second + string("\n");
  return result;
}

template <typename T, typename U> const vector<T> corpushl<T, U>::prej(const vector<corpushl<T, U> >& prejs) const {
  vector<T> result;
  for(int i = 0; i < prejs.size(); i ++)
    // XXX confirm me: need some other counting methods?
    // XXX fixme: amount of the word that is not said in the context is
    //            also important.
    result.push_back(distanceInUnion(match2relPseudo(prejs[i])) / distanceInUnion(*this));
  return result;
}

template <typename T, typename U> const T corpushl<T, U>::culturalConflicts(const corpushl<T, U>& base) const {
  corpushl<T, U> work(this->match2relPseudo(base));
  corpushl<T, U> workc(base.match2relPseudo(*this));
  corpushl<T, U> dwork(*this - work);
  corpushl<T, U> dworkc(base - workc);
  // XXX this is broken method. DO NOT USE THIS ONLY.: fixme...
  return (dwork.distanceInUnion(dwork) - dworkc.distanceInUnion(dworkc) + work.distanceInUnion(work) + workc.distanceInUnion(workc)) / (this->distanceInUnion(*this) + base.distanceInUnion(base));
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::conflictPart() {
  corpushl<T, U> result;
  // search conflict parts.
  // dictionary base of the word 'NOT' is needed.
  return;
}

template <typename T, typename U> const string corpushl<T, U>::optToc() {
  string result;
  // optimize corpust to make minimized toc.
  // stub.
  return result;
}

template <typename T, typename U> const string corpushl<T, U>::serialize(const T& thresh, const T& threshall, const T& threshm) const {
  cerr << "serialize" << endl;
  string result;
  string resultm;
  T th(thresh);
  T MM(0.);
  vector<int> looked;
  if(corpust.rows() < 1 || corpust.cols() < 1) {
    result = string("*NULL*");
    return result;
  }
  for( ; threshall < th; th *= thresh) {
    cerr << "." << flush;
    vector<pair<T, int> > MMM;
    for(int kk = 0; kk < corpust(0, 0).size(); kk ++) {
      MMM.push_back(make_pair(0, kk));
      for(int i = 0; i < corpust.rows(); i ++)
        for(int j = 0; j < corpust.cols(); j ++)
          MMM[kk].first = max(MMM[kk].first, abs(corpust(i, j)[kk]));
    }
    sort(MMM.begin(), MMM.end());
    MM = max(MM, MMM[MMM.size() - 1].first);
    for(int kk = 0; kk < MMM.size(); kk ++) {
      vector<pair<T, pair<pair<int, int>, int> > > idxs;
      const int k(MMM[MMM.size() - 1 - kk].second);
      if(binary_search(looked.begin(), looked.end(), k))
        continue;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int i = 0; i < corpust.rows(); i ++) if(!binary_search(looked.begin(), looked.end(), i))
        for(int j = 0; j < corpust.cols(); j ++) if(!binary_search(looked.begin(), looked.end(), j))
          if(threshm <= abs(corpust(i, j)[k]) &&
             MM * th <= abs(corpust(i, j)[k]) &&
             abs(corpust(i, j)[k]) <= MM * th / thresh)
            idxs.push_back(make_pair(abs(corpust(i, j)[k]),
                             make_pair(make_pair(i, j), k)));
      sort(idxs.begin(), idxs.end());
      if(!idxs.size()) continue;
      string buf;
      buf = serializeSub(idxs, true);
      if(buf.size())
        result  += buf;
      buf = serializeSub(idxs, false);
      if(buf.size())
        resultm += buf;
      for(int i = 0; i < idxs.size(); i ++) {
        looked.push_back(idxs[i].second.first.first);
        looked.push_back(idxs[i].second.first.second);
        looked.push_back(idxs[i].second.second);
      }
      sort(looked.begin(), looked.end());
      looked.erase(unique(looked.begin(), looked.end()), looked.end());
    }
  }
  return result + "\n-" + resultm;
}

template <typename T, typename U> const string corpushl<T, U>::serializeSub(const vector<pair<T, pair<pair<int, int>, int> > >& idxs, const bool& sign) const {
  string result;
  if(idxs.size() <= 0)
    return result;
  vector<pair<T, int> > iwords;
  vector<pair<T, int> > lwords;
  for(int i = 0; i < idxs.size(); i ++) {
    const int& ii(idxs[idxs.size() - 1 - i].second.first.first);
    // N.B. added order.
    const int& kk(idxs[idxs.size() - 1 - i].second.first.second);
    const int& jj(idxs[idxs.size() - 1 - i].second.second);
    if(sign ? corpust(ii, kk)[jj] <= T(0) :
              corpust(ii, kk)[jj] >= T(0))
      continue;
    int idx(- 1);
    for(int j = 0; j < iwords.size(); j ++)
      if(iwords[j].second == ii) {
        idx = j;
        break;
      }
    T val(0);
    if(jj < corpust.cols() && ii < corpust(0, jj).size())
      for(int k = 0; k < corpust.rows(); k ++) {
        if(corpust(k, jj).size())
          val += corpust(k, jj)[ii];
      }
    else
      val = T(1);
    if(idx < 0)
      iwords.push_back(make_pair(corpust(ii, kk)[jj] * val, ii));
    else
      iwords[idx].first += corpust(ii, kk)[jj] * val;
    idx = - 1;
    for(int j = 0; j < lwords.size(); j ++)
      if(lwords[j].second == kk) {
        idx = j;
        break;
      }
    val = T(0);
    if(jj < corpust.rows() && kk < corpust(jj, 0).size())
      for(int k = 0; k < corpust.cols(); k ++) {
        if(corpust(jj, k).size())
          val += corpust(jj, k)[kk];
      }
    else
      val = T(1);
    if(idx < 0)
      lwords.push_back(make_pair(corpust(ii, kk)[jj] * val, kk));
    else
      lwords[idx].first += corpust(ii, kk)[jj] * val;
  }
  sort(iwords.begin(), iwords.end());
  sort(lwords.begin(), lwords.end());
  vector<int> id, ld;
  for(int k = 0; k < iwords.size(); k ++)
    for(int kk = 0; kk < lwords.size(); kk ++)
      if(iwords[k].second == lwords[kk].second) {
        if(iwords[k].first > lwords[kk].first)
          ld.push_back(lwords.size() - 1 - kk);
        else
          id.push_back(k);
      }
  // XXX is this needed?
  for(int k = 0; k < iwords.size(); k ++)
    if(iwords[k].second == idxs[0].second.second)
      id.push_back(k);
  for(int k = 0; k < lwords.size(); k ++)
    if(lwords[k].second == idxs[0].second.second)
      ld.push_back(k);
  sort(id.begin(), id.end());
  sort(ld.begin(), ld.end());
  id.erase(unique(id.begin(), id.end()), id.end());
  ld.erase(unique(ld.begin(), ld.end()), ld.end());
  for(int k = 0; k < ld.size(); k ++)
    lwords.erase(lwords.begin() + ld[k] - k);
  for(int k = 0; k < id.size(); k ++)
    iwords.erase(iwords.begin() + id[k] - k);
  for(int k = 0; k < iwords.size(); k ++)
    result += words[iwords[iwords.size() - 1 - k].second];
  result += words[idxs[0].second.second];
  for(int k = 0; k < lwords.size(); k ++)
    result += words[lwords[k].second];
  result += string(".");
  return result;
}

template <typename T, typename U> const string corpushl<T, U>::summary(const vector<string>& words, const vector<corpushl<T, U> >& meanings, const T& thresh) const {
  string    result;
  vector<T> score;
  for(int i = 0; i < words.size(); i ++) {
    corpushl<T, U> mn(*this);
    mn = mn.abbrev(words[i], meanings[i]) - *this;
    score.push_back(mn.distanceInUnion(mn));
  }
  vector<T> sscore(score);
  sort(sscore.begin(), sscore.end());
  for(int i = 0; i < sscore.size(); i ++)
    for(int j = 0; j < score.size(); j ++)
      if(sscore[sscore.size() - i - 1] == score[j] && sscore[sscore.size() - i - 1] > thresh) {
        result += words[j];
        result += string("(") + to_string(score[j]) + string(")");
      }
  result += string("\n");
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::abbrev(const string word, const corpushl<T, U>& mean) {
  cerr << "abbrev 0/2" << endl;
  corpushl<T, U> work(this->match2relPseudo(mean));
  // XXX fixme: need to add abbreved word for work.
  cerr << "abbrev 1/2" << endl;
  const T t(distanceInUnion(work) / work.distanceInUnion(work));
  cerr << "abbrev 2/2" << endl;
  if(isfinite(t))
    return *this - (work * t);
  // can't abbrev.
  cerr << "abbrev : XXX : nan" << endl;
  return *this;
}

template <typename T, typename U> const vector<string> corpushl<T, U>::reverseLink(const corpushl<T, U>& orig) const {
  vector<string> res;
  vector<int>    ridx0, ridx1;
  corpushl<T, U> work(this->match2relPseudo(orig));
  vector<string> rwords(gatherWords(words, work.words, ridx0, ridx1));
  res.push_back(string(" (") + to_string(work.distanceInUnion(work)) + string(")"));
  for(int i = 0; i < rwords.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0)
    for(int j = 0; j < rwords.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < rwords.size(); k ++)
        if(ridx0[k] >= 0 && ridx1[k] >= 0 &&
           abs(work.corpust(ridx1[i], ridx1[j])[ridx1[k]]) != T(0) )
          res.push_back(work.words[i] + string("-") +
                        work.words[j] + string("-") +
                        work.words[k] + string(" @ ") +
                        to_string(work.corpust(ridx1[i], ridx1[j])[ridx1[k]]));
  return res;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::reDig(const T& ratio) {
#if defined(_OPENMP)
#pragma omp paralell for schedule(static, 1)
#endif
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      for(int k = 0; k < corpust(i, j).size(); k ++)
        corpust(i, j)[k] = (corpust(i, j)[k] < T(0) ? - T(1) : T(1)) * exp(log(abs(corpust(i, j)[k])) * ratio);
  return *this;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::simpleThresh(const T& ratio) {
  T absmax(0);
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      for(int k = 0; k < corpust(i, j).size(); k ++)
        if(absmax < corpust(i, j)[k])
          absmax = abs(corpust(i, j)[k]);
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      for(int k = 0; k < corpust(i, j).size(); k ++)
        if(abs(corpust(i, j)[k]) < absmax * ratio)
          corpust(i, j)[k] = T(0);
  return *this;
}

template <typename T, typename U> vector<string> corpushl<T, U>::gatherWords(const vector<string>& in0, const vector<string>& in1, vector<int>& ridx0, vector<int>& ridx1) const {
  cerr << "gatherWords" << endl;
  vector<string> result;
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
  vector<string> sin0(in0), sin1(in1);
  sort(sin0.begin(), sin0.end());
  sort(sin1.begin(), sin1.end());
  int rbufsize(in0.size() + in1.size() + 1);
  for(int i = 0; i < rbufsize; i ++) {
    ridx0.push_back(- 1);
    ridx1.push_back(- 1);
  }
  int j = 0;
  for(int i = 0; i < sin0.size(); i ++) {
    for(; i < sin0.size() && j < sin1.size(); j ++) {
      if(sin0[i] == sin1[j]) {
        result.push_back(string(sin0[i]));
        const int d0(distance(in0.begin(), equal_range(in0.begin(), in0.end(), sin0[i]).first));
        const int d1(distance(in1.begin(), equal_range(in1.begin(), in1.end(), sin1[j]).first));
        if(0 <= d0 && d0 < in0.size())
          ridx0[result.size() - 1] = d0;
        if(0 <= d1 && d1 < in1.size())
          ridx1[result.size() - 1] = d1;
        i ++;
        continue;
      } else if(sin0[i] > sin1[j]) {
        result.push_back(string(sin1[j]));
        const int d1(distance(in1.begin(), equal_range(in1.begin(), in1.end(), sin1[j]).first));
        if(0 <= d1 && d1 < in1.size())
          ridx1[result.size() - 1] = d1;
        continue;
      } else
        j --;
      j ++;
      break;
    }
    if(sin0.size() <= i) break;
    result.push_back(string(sin0[i]));
    const int d0(distance(in0.begin(), equal_range(in0.begin(), in0.end(), sin0[i]).first));
    if(0 <= d0 && d0 < in0.size())
      ridx0[result.size() - 1] = distance(in0.begin(), equal_range(in0.begin(), in0.end(), sin0[i]).first);
  }
  return result;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic> corpushl<T, U>::prepareDetail(const corpushl<T, U>& other, const vector<string>& workwords, const int& eidx, const vector<int>& ridx0, const vector<int>& ridx1, const vector<int>& ridx2) {
  cerr << "prepareDetail 0/2" << endl;
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
  for(int i = 0; i < workwords.size(); i ++) if(ridx2[i] >= 0)
    for(int j = 0; j < workwords.size(); j ++) if(ridx2[j] >= 0) {
      for(int k = 0; k < workwords.size(); k ++) {
        const T ratio(corpust(ridx0[i], ridx0[j])[ridx0[k]]);
        if(ridx2[k] >= 0 && ridx1[k] >= 0 && ridx1[i] >= 0 && ridx1[j] >= 0)
          res(ridx2[i], ridx2[j])[ridx2[k]] +=
            ratio * other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
      }
    }

  // Sum-up words defined in definition row, col without crossing point.
  {
    const int i(eidx);
    for(int j = 0; j < workwords.size(); j ++) if(ridx2[j] >= 0) {
      cerr << "prepareDetail 1/2 " << i << "/" << workwords.size() << endl;
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
        } else if(ridx1[j] >= 0 && ridx1[k] >= 0) {
          res(ridx2[i], ridx2[j])[ridx2[k]] += other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
          res(ridx2[j], ridx2[i])[ridx2[k]] += other.corpust(ridx1[j], ridx1[i])[ridx1[k]];
        }
      }
    }
  }
  return res;
}

template <typename T, typename U> const vector<string>& corpushl<T, U>::getWords() const {
  return words;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic>& corpushl<T, U>::getCorpus() const {
  return corpust;
}

#define _CORPUS_
#endif


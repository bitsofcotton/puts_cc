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
  vector<int> head, tail;
  words = vector<string>();
  words.push_back(string("^"));
  head.push_back(0);
  ptrs.push_back(head);
  for(auto itr = words0.begin(); itr != words0.end(); ++ itr) {
    const int idx = distance(words0.begin(), itr);
    if(itr->size() && ptrs0[idx].size()) {
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
          // const T buf0(abs(*itr + .5 - ptrs[i][ctru]));
          // const T buf1(abs(*itr + .5 - ptrs[j][ctrv]));
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
  const string          toc(const vector<corpushl<T, U> >& base, const vector<string>& words, const string& mwords) const;
  const string          toc(const vector<vector<corpushl<T, U> > >& base, const vector<string>& words, const string& mwords) const;
  const vector<T>       prej(const vector<corpushl<T, U> >& prejs) const;
  const T               culturalConflicts(const corpushl<T, U>& base) const;
  const corpushl<T, U>  conflictPart();
  const vector<string>& getWords() const;
  const Tensor&         getCorpus() const;
  const string          serialize() const;
  const string          summary(const vector<string>& words, const vector<corpushl<T, U> >& meanings, const T& thresh) const;
  const corpushl<T, U>  abbrev(const string word, const corpushl<T, U>& mean);
  const vector<string>  reverseLink(const corpushl<T, U>& orig) const;
  const pair<T, T>      compareStructure(const corpushl<T, U>& src, const T& thresh = T(1e-4), const T& thresh2 = T(.125)) const;
  const corpushl<T, U>  reDig(const T& ratio);
  const corpushl<T, U>  simpleThresh(const T& ratio);
private:
  const string          serializeSub(const vector<int>& idxs) const;
  const Vec             singularValues() const;
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
      if(0 <= ridx0[i] && 0 <= ridx0[j])
        for(int k = 0; k < result.corpust(i, j).size(); k ++) {
          result.corpust(i, j)[k] = 0.;
          if(0 <= ridx0[k])
            result.corpust(i, j)[k] += corpust(ridx0[i], ridx0[j])[ridx0[k]];
        }
      if(0 <= ridx1[i] && 0 <= ridx1[j])
        for(int k = 0; k < result.corpust(i, j).size(); k ++)
          if(0 <= ridx1[k])
            result.corpust(i, j)[k] += other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
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
  if(!flag1 || eidx < 0 || ridx0[eidx] < 0)
    return result;
  result.corpust = Tensor(workwords.size() - 1, workwords.size() - 1);
  for(int i = 0; i < workwords.size(); i ++)
    if(ridx0[i] >= 0 && ridx2[i] >= 0)
      result.words.push_back(words[ridx0[i]]);
    else if(ridx1[i] >= 0 && ridx2[i] >= 0)
      result.words.push_back(other.words[ridx1[i]]);
#if defined(_OPENMP)
#pragma omp paralell for schedule(static, 1)
#endif
  for(int i = 0; i < workwords.size(); i ++) if(0 <= ridx2[i])
    for(int j = 0; j < workwords.size(); j ++) if(0 <= ridx2[j]) {
      result.corpust(ridx2[i], ridx2[j]) = Vec(workwords.size() - 1);
      for(int k = 0; k < result.corpust(ridx2[i], ridx2[j]).size(); k ++) if(ridx2[k] >= 0) {
        if(ridx0[i] >= 0 && ridx0[j] >= 0 && ridx0[k] >= 0)
          result.corpust(ridx2[i], ridx2[j])[ridx2[k]] = corpust(ridx0[i], ridx0[j])[ridx0[k]];
        else
          result.corpust(ridx2[i], ridx2[j])[ridx2[k]] = T(0);
      }
    }
  cerr << "withDetail: adding detail table." << endl;
  result.corpust += prepareDetail(other, workwords, eidx, ridx0, ridx1, ridx2);
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

template <typename T, typename U> const string corpushl<T, U>::toc(const vector<corpushl<T, U> >& base, const vector<string>& words, const string& mwords) const {
  vector<vector<corpushl<T, U> > > basev;
  for(int i = 0; i < base.size(); i ++) {
    vector<corpushl<T, U> > work;
    work.push_back(base[i]);
    basev.push_back(work);
  }
  return toc(basev, words, mwords);
}

template <typename T, typename U> const string corpushl<T, U>::toc(const vector<vector<corpushl<T, U> > >& base, const vector<string>& words, const string& mwords) const {
  string result;
  result += mwords + string(" : ");
  vector<pair<T, string> > mlist;
  for(int k = 0; k < base.size(); k ++) {
    corpushl<T, U> work(*this);
    for(int l = 0; l < base[k].size(); l ++)
      work = work.abbrev(words[k], base[k][l]);
    pair<T, T> cr(compareStructure(work));
    const T ncr(sqrt(cr.first * cr.first + cr.second * cr.second));
    cr.first  /= ncr;
    cr.second /= ncr;
    mlist.push_back(make_pair(work.distanceInUnion(*this) /
                              sqrt(distanceInUnion(*this)),
                              work.serialize() + string("(sr:)") +
                              to_string(cr.first) + string(", ") +
                              to_string(cr.second) ) );
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
  corpushl<T, U> work(match2relPseudo(base));
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

template <typename T, typename U> const string corpushl<T, U>::serialize() const {
  cerr << "serialize" << endl;
  if(corpust.rows() < 1 || corpust.cols() < 1)
    return string("*NULL*");
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
  return plus.serializeSub(entire)  + string(".\n-") +
         minus.serializeSub(entire) + string(".");
}

template <typename T, typename U> const string corpushl<T, U>::serializeSub(const vector<int>& idxs) const {
  cerr << "." << flush;
  vector<pair<T, int> > cscore;
  for(int i = 0; i < idxs.size(); i ++)
    if(0 <= idxs[i] && idxs[i] < corpust.rows()) {
      T tl(0), tc(0), tr(0);
      for(int j = 0; j < idxs.size(); j ++)
        // corpust(i0, i2)[i1].
        // left-wins: sum corpust(i0, k)[i1] > sum corpust(k, i0)[i1].
        for(int k = 0; k < idxs.size(); k ++) {
          if(0 <= idxs[k] && idxs[k] < corpust.cols() &&
             0 <= idxs[j] && idxs[j] < corpust(idxs[i], idxs[k]).size())
            tl += pow(corpust(idxs[i], idxs[k])[idxs[j]], T(2));
          if(0 <= idxs[j] && idxs[j] < corpust.rows() &&
             0 <= idxs[k] && idxs[k] < corpust.cols() &&
             0 <= idxs[i] && idxs[i] < corpust(idxs[j], idxs[k]).size())
            tc += pow(corpust(idxs[j], idxs[k])[idxs[i]], T(2));
          if(0 <= idxs[k] && idxs[k] < corpust.rows() &&
             0 <= idxs[j] && idxs[j] < corpust(idxs[k], idxs[i]).size())
            tr += pow(corpust(idxs[k], idxs[i])[idxs[j]], T(2));
        }
      const T lscore(max(tr, tl) / (tc + 1));
      if(isfinite(lscore))
        cscore.push_back(make_pair(lscore, idxs[i]));
    }
  if(cscore.size() <= 0)
    return string();
  if(cscore.size() <= 1)
    return words[cscore[0].second];
  if(cscore.size() <= 2)
    return words[cscore[0].second] + words[cscore[1].second];
  sort(cscore.begin(), cscore.end());
  vector<int> left, right;
  const int& middle(cscore[0].second);
  vector<pair<T, int> > score;
  for(int i = 0; i < idxs.size(); i ++)
    if(idxs[i] != middle && 0 <= idxs[i] && idxs[i] < corpust.rows()) {
      T lscore(0);
      for(int j = 0; j < idxs.size(); j ++)
        if(idxs[j] != middle && 0 <= idxs[j] && idxs[j] < corpust.cols()) {
          lscore += corpust(idxs[i], idxs[j])[middle];
          lscore -= corpust(idxs[j], idxs[i])[middle];
        }
      if(isfinite(lscore))
        score.push_back(make_pair(lscore, idxs[i]));
    }
  sort(score.begin(), score.end());
  for(int i = 0; i < score.size() / 2; i ++)
    left.push_back(score[i].second);
  for(int i = score.size() / 2; i < score.size(); i ++)
    right.push_back(score[i].second);
  return serializeSub(left) + words[middle] + serializeSub(right);
}

template <typename T, typename U> const string corpushl<T, U>::summary(const vector<string>& words, const vector<corpushl<T, U> >& meanings, const T& thresh) const {
  string                   result;
  vector<pair<T, string> > score;
  for(int i = 0; i < words.size(); i ++) {
    corpushl<T, U> mn(*this);
    mn = mn.abbrev(words[i], meanings[i]) - *this;
    const T lscore(mn.distanceInUnion(mn));
    if(isfinite(lscore))
      score.push_back(make_pair(lscore, words[i]));
  }
  sort(score.begin(), score.end());
  for(int i = 0; i < score.size(); i ++) {
    const int ii(score.size() - i - 1);
    if(score[ii].first > thresh) {
      result += score[ii].second;
      result += string("(") + to_string(score[ii].first) + string(")");
    }
  }
  result += string("\n");
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::abbrev(const string word, const corpushl<T, U>& mean) {
  cerr << "abbrev 0/2" << endl;
  corpushl<T, U> work(match2relPseudo(mean));
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
  corpushl<T, U> work(match2relPseudo(orig));
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

template <typename T, typename U> const Eigen::Matrix<T, Eigen::Dynamic, 1> corpushl<T, U>::singularValues() const {
  Mat planes(corpust.rows(), corpust.cols());
  for(int i = 0; i < corpust.rows(); i ++) {
    Mat buf(corpust.rows(), corpust.cols());
    for(int j = 0; j < corpust.cols(); j ++) {
      for(int k = 0; k < corpust(i, j).size(); k ++) {
        if(!isfinite(corpust(i, j)[k])) {
          cerr << "singularValues() : nan" << endl;
          for(int i = 0; i < planes.rows(); i ++)
            for(int j = 0; j < planes.cols(); j ++)
              planes(i, j) = T(0);
          Eigen::JacobiSVD<Mat> svd(planes, 0);
          return svd.singularValues();
        }
        buf(k, j) = corpust(i, j)[k];
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

template <typename T, typename U> const vector<string>& corpushl<T, U>::getWords() const {
  return words;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic>& corpushl<T, U>::getCorpus() const {
  return corpust;
}




template <typename T, typename U> const string optimizeTOC(const string& input, const U* words, const vector<string>& dict, vector<string>& detailwords, const vector<string>& delimiter, const int& szwindow, const int& depth, const T& redig = T(1)) {
  // prepare dictionaries:
  cerr << "optimizeToc: parsing input" << flush;
  vector<vector<corpushl<T, U> > > details;
  assert(dict.size() == detailwords.size());
  for(int i = 0; i < dict.size(); i ++) {
    cerr << "." << flush;
    details.push_back(vector<corpushl<double, char> >());
    for(int j = 0; j < dict[i].size() / szwindow + 1; j ++) {
      corpus<double, char> cstat;
      cstat.init(words, 0, 120);
      cstat.compute(dict[i].substr(j * szwindow, szwindow).c_str(), delimiter);
      details[details.size() - 1].push_back(corpushl<double, char>(cstat));
    }
  }
  vector<corpushl<double, char> > cstat0;
  for(int i = 0; i < input.size() / szwindow + 1; i ++) {
    corpus<double, char> cstat;
    cstat.init(words, 0, 120);
    cstat.compute(input.substr(i * szwindow, szwindow).c_str(), delimiter);
    cstat0.push_back(corpushl<double, char>(cstat));
  }
  
  cerr << "analysing input text." << flush;
  for(int i = 0; i < details.size(); i ++) {
    cerr << "." << flush;
    for(int j = 0; j < cstat0.size(); j ++)
      for(int k = 0; k < details[i].size(); k ++)
        cstat0[j] = cstat0[j].withDetail(detailwords[i], details[i][k]);
  }
  vector<vector<pair<double, int> > > cstats;
  for(int i = 0; i < cstat0.size(); i ++) {
    cerr << "." << flush;
    cstats.push_back(vector<pair<double, int> >());
    for(int j = 0; j < i; j ++)
      cstats[i].push_back(cstats[j][i]);
    for(int j = i; j < cstat0.size(); j ++) {
      cerr << "." << flush;
      cstats[i].push_back(make_pair(cstat0[i].distanceInUnion(cstat0[j]), j));
    }
  }
  for(int i = 0; i < cstats.size(); i ++)
    sort(cstats[i].begin(), cstats[i].end());
  
  cerr << "OK, making outputs." << flush;
  vector<int> phrases;
  cerr << "." << flush;
  vector<pair<double, int> > work;
  vector<vector<int> >       idxs;
  for(int j = 0; j < cstats.size(); j ++) {
    double score(0);
    idxs.push_back(vector<int>());
    for(int k = 0, kk = 0; k < cstats.size() && kk < cstats.size() / depth; k ++) {
      const int k0(cstats.size() - k - 1);
      if(!binary_search(phrases.begin(), phrases.end(), cstats[j][k0].second)) {
        score += cstats[j][k0].first;
        idxs[j].push_back(cstats[j][k0].second);
        phrases.push_back(cstats[j][k0].second);
        sort(phrases.begin(), phrases.end());
        kk ++;
      }
    }
    work.push_back(make_pair(score, j));
  }
  sort(work.begin(), work.end());
  string result;
  for(int j = 0; j < work.size(); j ++) {
    const int jj(work[work.size() - j - 1].second);
    if(idxs[jj].size() <= 0)
      break;
    corpushl<double, char> cs(cstat0[jj]);
    for(int j = 0; j < idxs[jj].size(); j ++)
      cs += cstat0[idxs[jj][idxs[jj].size() - j - 1]];
    cs.reDig(redig);
    result += cs.serialize() + string("\n");
  }
  return result;
}

#define _CORPUS_
#endif


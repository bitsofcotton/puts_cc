#if !defined(_CORPUS_HL_)

#include <Eigen/Core>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iterator>
#include <iostream>
#include "corpus.hh"

using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::sort;
using std::find;
using std::distance;
using std::pair;
using std::to_string;
using std::isfinite;
using std::sqrt;

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
  const corpushl<T, U>  operator +  (const corpushl<T, U>& other) const;
  const corpushl<T, U>  operator -  () const;
  const corpushl<T, U>  operator -  (const corpushl<T, U>& other) const;
  const corpushl<T, U>  operator *  (const T& t)                  const;
  const corpushl<T, U>& operator =  (const corpushl<T, U>& other);
  const corpushl<T, U>  withDetail(const string& word, const corpushl<T, U>& other);
  const T               distanceInUnion(const corpushl<T, U>& other) const;
  const corpushl<T, U>  match2relPseudo(const corpushl<T, U>& other) const;
  const string          toc(const vector<corpushl<T, U> >& base, const vector<string>& words, const vector<corpushl<T, U> >& meanings, const vector<string>& mwords, const int& nloop, const T& thresh) const;
  const vector<T>       prej(const vector<corpushl<T, U> >& prejs) const;
  const T               culturalConflicts(const corpushl<T, U>& base) const;
  const corpushl<T, U>  conflictPart();
  const string          optToc();
  const vector<string>& getWords() const;
  const Tensor&         getCorpus() const;
  const string          summary(const vector<string>& words, const vector<corpushl<T, U> >& meanings, const T& thresh) const;
  const string          summary0(const T& thresh) const;
  const corpushl<T, U>  abbrev(const string word, const corpushl<T, U>& mean);
  const vector<string>  reverseLink(const corpushl<T, U>& orig) const;
  const corpushl<T, U>  reDig(const corpushl<T, U>& src, const T& ratio);
  const corpushl<T, U>  simpleThresh(const T& ratio);
private:
  vector<string>                   words;
  Tensor                           corpust;
  typedef vector<string>::iterator vsitr;
   
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

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator + (const corpushl<T, U>& other) const {
  corpushl<T, U> result;
  vector<int>    ridx0, ridx1;
  result.words   = gatherWords(words, other.words, ridx0, ridx1);
  result.corpust = Tensor(result.words.size(), result.words.size());
  cerr << "operator +" << endl;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < result.corpust.rows(); i ++) {
    for(int j = 0; j < result.corpust.cols(); j ++) {
      result.corpust(i, j) = Vec(result.words.size());
      if(!((0 <= ridx0[i] && 0 <= ridx0[j]) || (0 <= ridx1[i] && 0 <= ridx1[j])))
        for(int k = 0; k < result.corpust(i, j).size(); k ++)
          result.corpust(i, j)[k] = 0;
      else
        for(int k = 0; k < result.corpust(i, j).size(); k ++) {
          result.corpust(i, j)[k] = 0.;
          if(0 <= ridx0[i] && 0 <= ridx0[j] && 0 <= ridx0[k])
            result.corpust(i, j)[k] += corpust(ridx0[i], ridx0[j])[ridx0[k]];
          if(0 <= ridx1[i] && 0 <= ridx1[j] && 0 <= ridx1[k])
            result.corpust(i, j)[k] += other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
        }
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
  corpushl<T, U> result;
  cerr << "withDetail : enter" << endl;
  if(words.size() <= 0 || other.words.size() <= 0)
    return result;
  auto itr(find(words.begin(), words.end(), word));
  int fidx(distance(words.begin(), itr));
  if(!(0 <= fidx && fidx < words.size()))
    return result;
  std::cout << *itr << ", " << word << ": " << itr->size() << " / " << word.size() << endl;
  if(*itr != word)
    return result;
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
  // XXX : only counts certain tables.
  for(int i = 0; i < drop.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0) {
    for(int j = 0; j < drop.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < drop.size(); k ++) if(ridx0[i] >= 0 && ridx1[k] >= 0)
        res += corpust(ridx0[i], ridx0[j])[ridx0[k]] * other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
  }
  return res;
}

template <typename T, typename U> const string corpushl<T, U>::toc(const vector<corpushl<T, U> >& base, const vector<string>& words, const vector<corpushl<T, U> >& meanings, const vector<string>& mwords, const int& nloop, const T& thresh) const {
  string result;
  corpushl<T, U> work0(*this);
  for(int i = 0; i < meanings.size(); i ++) {
    corpushl<T, U> work(this->match2relPseudo(meanings[i]));
    work0 -= work;
    for(int j = 0; j < base.size(); j ++)
      work = work.abbrev(words[j], base[j]);
    result += mwords[i] + string("(") + to_string(work.distanceInUnion(work) / distanceInUnion(*this))+ string(")") + string(" : ") + work.summary0(thresh) + string("\n");
  }
  for(int j = 0; j < base.size(); j ++)
    work0 = work0.abbrev(words[j], base[j]);
  result += work0.summary0(thresh);
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

template <typename T, typename U> const string corpushl<T, U>::summary0(const T& thresh) const {
  if(corpust.rows() < 1 || corpust.cols() < 1)
    return string("null");
  string result;
  vector<pair<T, int> > scores;
  for(int i = 0; i < corpust(0, 0).size(); i ++) {
    pair<T, int> work;
    work.first  = 0;
    work.second = i;
    scores.push_back(work);
  }
  // XXX fixme:
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      for(int k = 0; k < corpust(i, j).size(); k ++)
        scores[k].first += corpust(i, j)[k];
  sort(scores.begin(), scores.end());
  for(int i = 0; i < scores.size(); i ++) {
    const int idx = scores.size() - i - 1;
    result += words[scores[idx].second];
    result += string("(") + to_string(scores[idx].first) + string(")");
  }
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::abbrev(const string word, const corpushl<T, U>& mean) {
  cerr << "abbrev 0/2" << endl;
  corpushl<T, U> work(this->match2relPseudo(mean));
  // XXX fixme: need to add abbreved word for work.
  vector<int>    ridx0, ridx1;
  vector<string> words(gatherWords(this->words, work.words, ridx0, ridx1));
  Vec orig(words.size() * words.size() * words.size());
  Vec delta(orig.size());
  cerr << "abbrev 1/2" << endl;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < words.size(); i ++) if(ridx0[i] >= 0 || ridx1[i] >= 0)
    for(int j = 0; j < words.size(); j ++) if(ridx0[j] >= 0 || ridx1[j] >= 0)
      for(int k = 0; k < words.size(); k ++) {
        const int idx(i * words.size() * words.size() + j * words.size() + k);
        if(ridx0[k] >= 0 || ridx1[k] >= 0) {
          if(ridx0[i] >= 0 && ridx0[j] >= 0 && ridx0[k] >= 0)
            orig[idx]  = corpust(ridx0[i], ridx0[j])[ridx0[k]];
          else
            orig[idx]  = 0;
          if(ridx1[i] >= 0 && ridx1[j] >= 0 && ridx1[k] >= 0)
            delta[idx] = work.corpust(ridx1[i], ridx1[j])[ridx1[k]];
          else
            delta[idx] = 0;
        } else
          orig[idx] = delta[idx] = 0;
      }
  const T t(orig.dot(orig) / sqrt(delta.dot(delta)));
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

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::reDig(const corpushl<T, U>& src, const T& ratio) {
  vector<int>    ridx0, ridx1;
  vector<string> words(gatherWords(words, src.words, ridx0, ridx1));
#if defined(_OPENMP)
#pragma omp paralell for schedule(static, 1)
#endif
  for(int i = 0; i < words.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0)
    for(int j = 0; j < words.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < words.size(); k ++) if(ridx0[k] >= 0 && ridx1[k] >= 0)
        corpust(ridx0[i], ridx0[j])[ridx0[k]] = (exp(corpust(ridx0[i], ridx0[j])[ridx0[k]]) - exp(T(1))) / exp(T(1));
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
  for(int i = 0; i < rbufsize; i ++)
    ridx0.push_back(- 1);
  for(int i = 0; i < rbufsize; i ++)
    ridx1.push_back(- 1);
  int j = 0;
  for(int i = 0; i < sin0.size(); i ++) {
    for(; j < sin1.size(); j ++) {
      if(sin0[i] == sin1[j])
        ridx1[result.size()] = distance(in1.begin(), find(in1.begin(), in1.end(), sin1[j]));
      else if(sin0[i] > sin1[j]) {
        result.push_back(string(sin1[j]));
        ridx1[result.size() - 1] = distance(in1.begin(), find(in1.begin(), in1.end(), sin1[j]));
        continue;
      } else
        j --;
      j ++;
      break;
    }
    const int dist(distance(in0.begin(), find(in0.begin(), in0.end(), sin0[i])));
    if(0 <= dist && dist < in0.size()) {
      result.push_back(string(sin0[i]));
      ridx0[result.size() - 1] = dist;
    }
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

#define _CORPUS_HL_
#endif


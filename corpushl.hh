#if !defined(_CORPUS_HL_)

#include <Eigen/Core>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iterator>
#include <iostream>
#include "corpus.hh"

template <typename T, typename U> class corpushl {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>   Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                Vec;
  typedef Eigen::Matrix<Vec, Eigen::Dynamic, Eigen::Dynamic> Tensor;
  
  corpushl();
  corpushl(const corpus<T, U>&   obj);
  corpushl(const corpushl<T, U>& obj);
  ~corpushl();
  
  const corpushl<T, U>            cast(const std::vector<std::string>& words) const;
  const corpushl<T, U>&           operator += (const corpushl<T, U>& other);
  const corpushl<T, U>&           operator -= (const corpushl<T, U>& other);
  const corpushl<T, U>&           operator *= (const T& t);
  const corpushl<T, U>            operator +  (const corpushl<T, U>& other) const;
  const corpushl<T, U>            operator -  () const;
  const corpushl<T, U>            operator -  (const corpushl<T, U>& other) const;
  const corpushl<T, U>            operator *  (const T& t)                  const;
  const corpushl<T, U>&           operator =  (const corpushl<T, U>& other);
  const corpushl<T, U>            withDetail(const std::string& word, const corpushl<T, U>& other);
  const T                         distanceInUnion(const corpushl<T, U>& other) const;
  const corpushl<T, U>            match2relPseudo(const corpushl<T, U>& other) const;
  const std::string               toc(const std::vector<corpushl<T, U> >& base, const std::vector<std::string>& words, const std::vector<corpushl<T, U> >& meanings, const std::vector<std::string>& mwords, const int& nloop, const T& thresh) const;
  const std::vector<T>            prej(const std::vector<corpushl<T, U> >& prejs) const;
  const T                         culturalConflicts(const corpushl<T, U>& base) const;
  const corpushl<T, U>            conflictPart();
  const std::string               optToc();
  const std::vector<std::string>& getWords() const;
  const Tensor&                   getCorpus() const;
  const std::string               summary(const std::vector<std::string>& words, const std::vector<corpushl<T, U> >& meanings, const T& thresh) const;
  const std::string               summary0(const T& thresh) const;
  const corpushl<T, U>            abbrev(const std::string word, const corpushl<T, U>& mean);
private:
  std::vector<std::string>        words;
  Tensor                          corpust;
  typedef std::vector<std::string>::iterator vsitr;
   
  std::vector<std::string>        gatherWords(const std::vector<std::string>& in0, const std::vector<std::string>& in1, std::vector<int>& ridx0, std::vector<int>& ridx1) const;
};

template <typename T, typename U> corpushl<T,U>::corpushl() {
  ;
}

template <typename T, typename U> corpushl<T,U>::~corpushl() {
  // auto called destructors for string.
  ;
}

template <typename T, typename U> corpushl<T,U>::corpushl(const corpus<T, U>& obj) {
  words   = std::vector<std::string>(obj.getWords());
  corpust = Tensor(obj.getCorpus());
}

template <typename T, typename U> corpushl<T,U>::corpushl(const corpushl<T, U>& obj) {
  words   = std::vector<std::string>(obj.words);
  corpust = Tensor(obj.corpust);
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::cast(const std::vector<std::string>& words) const {
  std::cerr << "cast 0/1" << std::endl;
  corpushl<T, U> result;
  std::vector<std::string> sword;
  std::sort(sword.begin(), sword.end());
  std::vector<int>         idxs;
  for(int i = 0; i < sword.size(); i ++) {
    auto p(std::find(words.begin(), words.end(), sword[i]));
    if(*p == sword[i]) {
      idxs.push_back(std::distance(words.begin(), p));
      result.words.push_back(*p);
    }
  }
  result.corpust = Tensor(idxs.size(), idxs.size());
  for(int i = 0; i < idxs.size(); i ++)
    for(int j = 0; j < idxs.size(); j ++)
      result.corpust(i, j) = corpust(idxs[i], idxs[j]);
  return result;
}

template <typename T, typename U> const corpushl<T, U>& corpushl<T, U>::operator = (const corpushl<T, U>& other) {
  words   = std::vector<std::string>(other.words);
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
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      corpust(i, j) *= t;
  return *this;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator + (const corpushl<T, U>& other) const {
  corpushl<T, U> result;
  std::vector<int> ridx0, ridx1;
  result.words   = gatherWords(words, other.words, ridx0, ridx1);
  result.corpust = Tensor(result.words.size(), result.words.size());
  std::cerr << "operator +" << std::endl;
  for(int i = 0; i < result.corpust.rows(); i ++) {
    for(int j = 0; j < result.corpust.cols(); j ++) {
      result.corpust(i, j) = Vec(result.words.size());
      if(!((0 <= ridx0[i] && 0 <= ridx0[j]) || (0 <= ridx1[i] && 0 <= ridx1[j])))
        for(int k = 0; k < result.corpust(i, j).size(); k ++)
          result.corpust(i, j)[k] = 0;
      else
        for(int k = 0; k < result.corpust(i, j).size(); k ++) {
          result.corpust(i, j)[k] = 0.;
          if(!((0 <= ridx0[i] && 0 <= ridx0[j] && 0 <= ridx0[k]) ||
               (0 <= ridx1[i] && 0 <= ridx1[j] && 0 <= ridx1[k])) )
            continue;
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
  std::cerr << "operator -" << std::endl;
  corpushl<T, U> result(*this);
  result.corpust = - result.corpust;
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator - (const corpushl<T, U>& other) const {
  std::cerr << "operator -" << std::endl;
  return (*this) + (- other);
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator * (const T& t) const {
  corpushl<T, U> work(*this);
  return work *= t;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::withDetail(const std::string& word, const corpushl<T, U>& other) {
  std::cerr << "withDetail 1/3" << std::endl;
  if(words.size() <= 0 || other.words.size() <= 0)
    return *this;
  //std::pair<vsitr, vsitr>;
  vsitr itr(std::find(words.begin(), words.end(), word));
  int fidx(std::distance(words.begin(), itr));
  if(!(0 <= fidx && fidx < words.size()))
    return *this;
  std::cout << *itr << ", " << word << ": " << itr->size() << " / " << word.size() << std::endl;
  if(*itr != word)
    return *this;
  corpushl<T, U>   result;
  std::vector<int> ridx0, ridx1;
  std::vector<std::string> workwords(gatherWords(words, other.words, ridx0, ridx1));
  result.corpust = Tensor(workwords.size() - 1, workwords.size() - 1);
  std::vector<int> ridx2;
  int eidx = - 1;
  bool flag1 = false;
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
  if(!flag1 || eidx < 0)
    return *this;
  for(int i = 0; i < workwords.size(); i ++) if(ridx2[i] >= 0)
    for(int j = 0; j < workwords.size(); j ++) if(ridx2[j] >= 0) {
      result.corpust(ridx2[i], ridx2[j]) = Vec(workwords.size() - 1);
      if(ridx0[i] >= 0 && ridx0[j] >= 0) {
        for(int k = 0; k < result.corpust(ridx2[i], ridx2[j]).size(); k ++) if(ridx2[k] >= 0) {
          if(ridx0[k] >= 0)
            result.corpust(ridx2[i], ridx2[j])[ridx2[k]] = corpust(ridx0[i], ridx0[j])[ridx0[k]];
          else
            result.corpust(ridx2[i], ridx2[j])[ridx2[k]] = 0;
        } else
          for(int k = 0; k < result.corpust(ridx2[i], ridx2[j]).size(); k ++)
            result.corpust(ridx2[i], ridx2[j])[k] = 0;
      } else
        for(int k = 0; k < result.corpust(ridx2[i], ridx2[j]).size(); k ++)
          result.corpust(ridx2[i], ridx2[j])[k] = 0;
    }
  for(int i = 0; i < workwords.size(); i ++) if(ridx2[i] >= 0) {
    result.words.push_back(workwords[i]);
    // Sum-up detailed word into result pool without definition itself.
    for(int j = 0; j < workwords.size(); j ++)
      if(ridx2[j] >= 0 && !(eidx == i && eidx == j)) {
        for(int k = 0; k < workwords.size(); k ++)
          if(ridx2[k] >= 0 && ridx1[k] >= 0 && ridx1[i] >= 0 && ridx1[j] >= 0)
            result.corpust(ridx2[i], ridx2[j])[ridx2[k]] += other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
      }
  }
  for(int i = 0; i < workwords.size(); i ++) if(ridx0[i] >= 0 && eidx != i) {
    std::cerr << "withDetail 2/3 : " << i << "/" << workwords.size() << std::endl;
    // Sum-up words defined in both given and detail pool without definition itself.
    for(int j = 0; j < workwords.size(); j ++) if(ridx0[j] >= 0 && eidx != j)
      for(int i1 = 0; i1 < workwords.size(); i1 ++) if(ridx1[i1] >= 0 && ridx2[i1] >= 0)
        for(int j1 = 0; j1 < workwords.size(); j1 ++) if(ridx1[j1] >= 0 && ridx2[j1] >= 0)
          for(int k = 0; k < workwords.size(); k ++) {
            if(ridx1[k] >= 0 && ridx2[k] >= 0)
              // XXX confirm me.
              result.corpust(ridx2[i1], ridx2[j1])[ridx2[k]] += corpust(ridx0[i], ridx0[j])[ridx0[eidx]] * other.corpust(ridx1[i1], ridx1[j1])[ridx1[k]];
          }
  }
  std::cerr << "withDetail 3/3" << std::endl;
  for(int i = 0; i < workwords.size(); i ++) if(ridx1[i] >= 0 && ridx2[i] >= 0) {
    // Sum-up definition itself.
    for(int j = 0; j < workwords.size(); j ++) if(ridx1[j] >= 0 && ridx2[j] >= 0)
      for(int k = 0; k < workwords.size(); k ++) if(ridx1[k] < 0 && ridx2[k] >= 0 && ridx0[k] >= 0)
        result.corpust(ridx2[i], ridx2[j])[ridx2[k]] += corpust(ridx0[eidx], ridx0[eidx])[ridx0[k]];
        // XXX confirm me :
        // * other.corpust(ridx1[i], ridx1[j])[ridx1[k]] * corpust(ridx0[eidx], ridx0[eidx])[ridx0[k]];
  }
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::match2relPseudo(const corpushl<T, U>& other) const {
  std::cerr << "match2relPseudo : 0/2" << std::endl;
  corpushl<T, U> result(*this);
  std::vector<int> ridx0, ridx1;
  std::vector<std::string> words(gatherWords(result.words, other.words, ridx0, ridx1));
  Tensor mul(result.corpust);
  for(int i = 0; i < mul.rows(); i ++)
    for(int j = 0; j < mul.cols(); j ++)
      mul(i, j) *= T(0);
  // XXX fixme:
  for(int i = 0; i < words.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0)
    for(int j = 0; j < words.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < words.size(); k ++) if(ridx0[k] >= 0) {
        mul(ridx0[i], ridx0[j])[ridx0[k]] += 1.;
        mul(ridx0[j], ridx0[k])[ridx0[i]] += 1.;
        mul(ridx0[k], ridx0[i])[ridx0[j]] += 1.;
      }
  for(int i = 0; i < result.corpust.rows(); i ++)
    for(int j = 0; j < result.corpust.cols(); j ++)
      for(int k = 0; k < result.corpust(i, j).size(); k ++)
        result.corpust(i, j)[k] *= mul(i, j)[k];
  return result;
}

template <typename T, typename U> const T corpushl<T, U>::distanceInUnion(const corpushl<T, U>& other) const {
  T res(0);
  std::vector<int> ridx0, ridx1;
  std::vector<std::string> drop(gatherWords(words, other.words, ridx0, ridx1));
  // XXX : only counts certain tables.
  for(int i = 0; i < drop.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0) {
    for(int j = 0; j < drop.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < drop.size(); k ++) if(ridx0[i] >= 0 && ridx1[k] >= 0)
        res += corpust(ridx0[i], ridx0[j])[ridx0[k]] * other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
  }
  return res;
}

template <typename T, typename U> const std::string corpushl<T, U>::toc(const std::vector<corpushl<T, U> >& base, const std::vector<std::string>& words, const std::vector<corpushl<T, U> >& meanings, const std::vector<std::string>& mwords, const int& nloop, const T& thresh) const {
  std::string result;
  corpushl<T, U> work0(*this);
  for(int i = 0; i < meanings.size(); i ++) {
    corpushl<T, U> work(this->match2relPseudo(meanings[i]));
    work0  -= work;
    for(int j = 0; j < base.size(); j ++)
      work = work.abbrev(words[j], base[j]);
    result += mwords[i] + std::string(" : ") + work.summary0(thresh) + std::string("\n");
  }
  for(int j = 0; j < base.size(); j ++)
    work0 = work0.abbrev(words[j], base[j]);
  result += work0.summary0(thresh);
  return result;
}

template <typename T, typename U> const std::vector<T> corpushl<T, U>::prej(const std::vector<corpushl<T, U> >& prejs) const {
  std::vector<T> result;
  for(int i = 0; i < prejs.size(); i ++)
    // XXX confirm me: need some other counting methods?
    // XXX fixme: amount of the word that is not said in the context is
    //            also important.
    result.push_back(distanceInUnion(match2relPseudo(prejs[i])));
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

template <typename T, typename U> const std::string corpushl<T, U>::optToc() {
  std::string result;
  // optimize corpust to make minimized toc.
  // stub.
  return result;
}

template <typename T, typename U> const std::string corpushl<T, U>::summary(const std::vector<std::string>& words, const std::vector<corpushl<T, U> >& meanings, const T& thresh) const {
  std::string result;
  std::vector<T> score;
  for(int i = 0; i < words.size(); i ++) {
    corpushl<T, U> mn(*this);
    mn = mn.abbrev(words[i], meanings[i]) - *this;
    score.push_back(mn.distanceInUnion(mn));
  }
  std::vector<T> sscore(score);
  std::sort(sscore.begin(), sscore.end());
  for(int i = 0; i < sscore.size(); i ++)
    for(int j = 0; j < score.size(); j ++)
      if(sscore[sscore.size() - i - 1] == score[j] && sscore[sscore.size() - i - 1] > thresh) {
        result += words[j];
        result += std::string("(") + std::to_string(score[j]) + std::string(")");
      }
  result += std::string("\n");
  return result;
}

template <typename T> int summarycmp(const std::pair<T, int>& x0, const std::pair<T, int>& x1) {
  return x0.first < x1.first;
}

template <typename T, typename U> const std::string corpushl<T, U>::summary0(const T& thresh) const {
  if(corpust.rows() < 1 || corpust.cols() < 1)
    return std::string("null");
  std::string result;
  std::vector<std::pair<T, int> > scores;
  for(int i = 0; i < corpust(0, 0).size(); i ++) {
    std::pair<T, int> work;
    work.first  = 0;
    work.second = i;
    scores.push_back(work);
  }
  for(int i = 0; i < corpust.rows(); i ++)
    for(int j = 0; j < corpust.cols(); j ++)
      for(int k = 0; k < corpust(i, j).size(); k ++) {
        scores[k].first += corpust(i, j)[k];
        // XXX fixme:
        scores[k].first += corpust(j, k)[i];
        scores[k].first += corpust(k, i)[j];
      }
  std::sort(scores.begin(), scores.end(), summarycmp<T>);
  for(int i = 0; i < scores.size(); i ++) {
    const int idx = scores.size() - i - 1;
    result += words[scores[idx].second];
    result += std::string("(") + std::to_string(scores[idx].first) + std::string(")");
  }
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::abbrev(const std::string word, const corpushl<T, U>& mean) {
  std::cerr << "abbrev 0/2" << std::endl;
  corpushl<T, U>   work(this->match2relPseudo(mean));
  // XXX fixme: need to add abbreved word for work.
  std::vector<int> ridx0, ridx1;
  std::vector<std::string> words = gatherWords(this->words, work.words, ridx0, ridx1);
  Vec orig(words.size() * words.size() * words.size());
  Vec delta(orig.size());
  std::cerr << "abbrev 1/2" << std::endl;
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
  const T t(orig.dot(delta) / std::sqrt(delta.dot(delta)));
  std::cerr << "abbrev 2/2" << std::endl;
  if(std::isfinite(t))
    return *this - (work * t);
  // can't abbrev.
  return *this;
}

template <typename T, typename U> std::vector<std::string> corpushl<T, U>::gatherWords(const std::vector<std::string>& in0, const std::vector<std::string>& in1, std::vector<int>& ridx0, std::vector<int>& ridx1) const {
  std::cerr << "gatherWords" << std::endl;
  std::vector<std::string> result;
  ridx0 = std::vector<int>();
  ridx1 = std::vector<int>();
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
  std::vector<std::string> sin0(in0), sin1(in1);
  std::sort(sin0.begin(), sin0.end());
  std::sort(sin1.begin(), sin1.end());
  int rbufsize(in0.size() + in1.size() + 1);
  for(int i = 0; i < rbufsize; i ++)
    ridx0.push_back(- 1);
  for(int i = 0; i < rbufsize; i ++)
    ridx1.push_back(- 1);
  int j = 0;
  for(int i = 0; i < sin0.size(); i ++) {
    for(; j < sin1.size(); j ++) {
      if(sin0[i] == sin1[j])
        ridx1[result.size()] = std::distance(in1.begin(), std::find(in1.begin(), in1.end(), sin1[j]));
      else if(sin0[i] > sin1[j]) {
        result.push_back(std::string(sin1[j]));
        ridx1[result.size() - 1] = std::distance(in1.begin(), std::find(in1.begin(), in1.end(), sin1[j]));
        continue;
      } else
        j --;
      j ++;
      break;
    }
    const int dist(std::distance(in0.begin(), std::find(in0.begin(), in0.end(), sin0[i])));
    if(0 <= dist && dist < in0.size()) {
      result.push_back(std::string(sin0[i]));
      ridx0[result.size() - 1] = dist;
    }
  }
  return result;
}

template <typename T, typename U> const std::vector<std::string>& corpushl<T, U>::getWords() const {
  return words;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic>& corpushl<T, U>::getCorpus() const {
  return corpust;
}

#define _CORPUS_HL_
#endif


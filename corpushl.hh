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
  const corpushl<T, U>            operator +  (const corpushl<T, U>& other);
  const corpushl<T, U>&           operator =  (const corpushl<T, U>& other);
  const corpushl<T, U>            withDetail(const std::string& word, const corpushl<T, U>& other);
  const T                         distanceInUnion(const corpushl<T, U>& other) const;
  const std::string               toc(const std::vector<corpushl<T, U> >& base);
  const std::string               optToc();
  const std::vector<std::string>& getWords() const;
  const Tensor&                   getCorpus() const;
private:
  std::vector<std::string>        words;
  Tensor                          corpust;
   
  std::vector<std::string> gatherWords(const std::vector<std::string>& in0, const std::vector<std::string>& in1, std::vector<int>& ridx0, std::vector<int>& ridx1);
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
  corpushl<T, U> result;
  std::vector<std::string> sword;
  std::sort(sword.begin(), sword.end());
  std::vector<int>         idxs;
  for(int i = 0; i < sword.size(); i ++) {
    auto p(std::equal_range(words.begin(), words.end(), sword[i]));
    if(*p.first == sword[i]) {
      idxs.push_back(std::distance(words.begin(), p.first));
      result.words.push_back(*p.first);
    }
  }
  result.corpust = Tensor(idxs.size(), idxs.size());
  for(int i = 0; i < idxs.size(); i ++)
    for(int j = 0; j < idxs.size(); j ++)
      result.corpust(i, j) = corpust(i, j);
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

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::operator + (const corpushl<T, U>& other) {
  corpushl<T, U> result;
  std::vector<int> ridx0, ridx1;
  result.words   = gatherWords(words, other.words, ridx0, ridx1);
  result.corpust = Tensor(result.words.size(), result.words.size());
  std::cout << other.words.size();
  for(int i = 0; i < result.corpust.rows(); i ++) {
    std::cerr << "operator + processing row: " << i << " / " << result.corpust.rows() << std::endl;
    for(int j = 0; j < result.corpust.cols(); j ++) {
      result.corpust(i, j) = Vec(result.words.size());
      if(!((0 <= ridx0[i] && 0 <= ridx0[j]) || (0 <= ridx1[i] && 0 <= ridx1[j])))
        for(int k = 0; k < result.corpust(i, j).size(); k ++)
          result.corpust(i, j)[k] = 0;
      else
        for(int k = 0; k < result.corpust(i, j).size(); k ++) {
          if(!((0 <= ridx0[i] && 0 <= ridx0[j] && 0 <= ridx0[k]) ||
               (0 <= ridx1[i] && 0 <= ridx1[j] && 0 <= ridx1[k])) ) {
            result.corpust(i, j)[k] = 0.;
            continue;
          }
          if(0 <= ridx0[i] && 0 <= ridx0[j] && 0 <= ridx0[k])
            result.corpust(i, j)[k] += corpust(ridx0[i], ridx0[j])[ridx0[k]];
          if(0 <= ridx1[i] && 0 <= ridx1[j] && 0 <= ridx1[k])
            result.corpust(i, j)[k] += other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
        }
    }
  }
  return result;
}

template <typename T, typename U> const corpushl<T, U> corpushl<T, U>::withDetail(const std::string& word, const corpushl<T, U>& other) {
  auto itr(std::equal_range(words.begin(), words.end(), word).first);
  if(*itr != word)
    return *this;
  corpushl<T, U>   result;
  std::vector<int> ridx0, ridx1;
  std::vector<std::string> workwords(gatherWords(words, other.words, ridx0, ridx1));
  result.corpust = Tensor(workwords.size() - 1, workwords.size() - 1);
  for(int i = 0; i < result.corpust.rows(); i ++)
    for(int j = 0; j < result.corpust.cols(); j ++)
      result.corpust(i, j) = Vec(workwords.size() - 1);
  std::vector<int> ridx2;
  int eidx = - 1;
  for(int i = 0, ii = 0; i < workwords.size(); i ++) {
    if(workwords[i] == word) {
      ridx2.push_back(- 1);
      eidx = i;
      continue;
    }
    ridx2.push_back(ii ++);
  }
  for(int i = 0; i < workwords.size(); i ++) {
    std::cerr << "withDetail 1/3 : " << i << "/" << workwords.size() << std::endl;
    if(ridx2[i] < 0)
      continue;
    result.words.push_back(workwords[i]);
    for(int j = 0; j < workwords.size(); j ++)
      if(ridx2[j] >= 0) {
        for(int k = 0; k < workwords.size(); k ++)
          if(ridx2[k] >= 0 && ridx0[k] >= 0 && ridx0[i] >= 0 && ridx0[j] >= 0)
            result.corpust(ridx2[i], ridx2[j])[ridx2[k]] += corpust(ridx0[i], ridx0[j])[ridx0[k]];
      }
  }
  for(int i = 0; i < workwords.size(); i ++) if(ridx0[i] >= 0) {
    std::cerr << "withDetail 2/3 : " << i << "/" << workwords.size() << std::endl;
    for(int j = 0; j < workwords.size(); j ++) if(ridx0[j] >= 0)
      for(int i1 = 0; i1 < workwords.size(); i1 ++) if(ridx1[i1] >= 0 && ridx2[i1] >= 0)
        for(int j1 = 0; j1 < workwords.size(); j1 ++) if(ridx1[j1] >= 0 && ridx2[j1] >= 0)
          for(int k = 0, kk = 0; k < workwords.size(); k ++) {
            if(ridx1[k] >= 0 && ridx2[k] >= 0)
              // XXX confirm me.
              result.corpust(ridx2[i1], ridx2[j1])[ridx2[k]] += corpust(ridx0[i], ridx0[j])[ridx0[eidx]] * other.corpust(ridx1[i1], ridx1[j1])[ridx1[k]];
          }
  }
  for(int i = 0; i < workwords.size(); i ++) if(ridx1[i] >= 0 && ridx2[i] >= 0) {
    std::cerr << "withDetail 3/3 : " << i << "/" << workwords.size() << std::endl;
    for(int j = 0; j < workwords.size(); j ++) if(ridx1[j] >= 0 && ridx2[j] >= 0)
      for(int k = 0; k < workwords.size(); k ++) if(ridx1[k] < 0 && ridx2[k] >= 0 && ridx0[k] >= 0)
        result.corpust(ridx2[i], ridx2[j])[ridx2[k]] += corpust(ridx0[eidx], ridx0[eidx])[ridx0[k]];
        // XXX fixme or not.
        // * other.corpust(ridx1[i], ridx1[j])[ridx1[k]] * corpust(ridx0[eidx], ridx0[eidx])[ridx0[k]];
  }
  return result;
}

template <typename T, typename U> const T corpushl<T, U>::distanceInUnion(const corpushl<T, U>& other) const {
  T res(0);
  std::vector<int> ridx0, ridx1;
  std::vector<std::string> drop(gatherWords(words, other.words, ridx0, ridx1));
  for(int i = 0; i < drop.size(); i ++) if(ridx0[i] >= 0 && ridx1[i] >= 0)
    for(int j = 0; j < drop.size(); j ++) if(ridx0[j] >= 0 && ridx1[j] >= 0)
      for(int k = 0; k < drop.size(); k ++) if(ridx0[i] >= 0 && ridx1[k] >= 0)
        res += corpust(ridx0[i], ridx0[j])[ridx0[k]] * other.corpust(ridx1[i], ridx1[j])[ridx1[k]];
  return res;
}

template <typename T, typename U> const std::string corpushl<T, U>::toc(const std::vector<corpushl<T, U> >& base) {
  std::string result;
  for(int i = 0; i < base.size(); i ++) {
    // stub.
    ;
  }
  return result;
}

template <typename T, typename U> const std::string corpushl<T, U>::optToc() {
  std::string result;
  // optimize corpust to make minimized toc.
  // stub.
  return result;
}

template <typename T, typename U> std::vector<std::string> corpushl<T, U>::gatherWords(const std::vector<std::string>& in0, const std::vector<std::string>& in1, std::vector<int>& ridx0, std::vector<int>& ridx1) {
  std::vector<std::string> sin0(in0), sin1(in1), result;
  std::sort(sin0.begin(), sin0.end());
  std::sort(sin1.begin(), sin1.end());
  ridx0 = std::vector<int>();
  ridx1 = std::vector<int>();
  int rbufsize(in0.size() + in1.size() + 1);
  for(int i = 0; i < rbufsize; i ++)
    ridx0.push_back(- 1);
  for(int i = 0; i < rbufsize; i ++)
    ridx1.push_back(- 1);
  int j = 0;
  for(int i = 0; i < sin0.size(); i ++) {
    for(; j < sin1.size(); j ++) {
      if(sin0[i] == sin1[j]) {
        ridx1[result.size() - 1] = std::distance(in1.begin(), std::equal_range(in1.begin(), in1.end(), sin1[j]).first);
        continue;
      } else if(sin0[i] < sin1[j])
        break;
      result.push_back(std::string(sin1[j]));
      ridx1[result.size() - 1] = std::distance(in1.begin(), std::equal_range(in1.begin(), in1.end(), sin1[j]).first);
    }
    result.push_back(std::string(sin0[i]));
    ridx0[result.size() - 1] = std::distance(in0.begin(), std::equal_range(in0.begin(), in0.end(), sin0[i]).first);
  }
  for(; j < sin1.size(); j ++) {
    result.push_back(std::string(sin1[j]));
    ridx1[result.size() - 1] = std::distance(in1.begin(), std::equal_range(in1.begin(), in1.end(), sin1[j]).first);
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


#if !defined(_CORPUS_)

#include <Eigen/Core>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iterator>
#include <iostream>

bool partial_compare(const std::string& x0, const std::string& x1);
bool partial_compare_equal(const std::string& x0, const std::string& x1);

template <typename T, typename U> class corpus {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>   Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1>                Vec;
  typedef Eigen::Matrix<Vec, Eigen::Dynamic, Eigen::Dynamic> Tensor;
  
  corpus();
  ~corpus();
  
  int  isin(const std::string& key);
  void init(const U* words, const int& nthresh, const int& Nthresh);
  corpus<T, U>&                   operator = (const corpus<T, U>& other);
  const void                      compute(const U* input);
  const std::vector<std::string>& getWords() const;
  const Tensor&                   getCorpus() const;
private:
  std::vector<std::string>       words0;
  std::vector<std::string>       words;
  std::vector<std::vector<int> > ptrs0;
  std::vector<std::vector<int> > ptrs;
  typedef std::vector<std::string>::iterator       vsitr;
  typedef std::vector<int>::iterator               viitr;
  typedef std::vector<std::vector<int> >::iterator vviitr;
  Tensor                         corpust;
  int                            nthresh;
  int                            Nthresh;
  T                              Midx;
   
  void getWordPtrs(const U* input);
  void corpusEach();
};

template <typename T, typename U> corpus<T,U>::corpus() {
  ;
}

template <typename T, typename U> corpus<T,U>::~corpus() {
  // auto called destructors for string.
  ;
}

template <typename T, typename U> void corpus<T,U>::init(const U* words, const int& nthresh, const int& Nthresh) {
  std::string buf;
  bool        flag  = false;
  for(int i = 0; words[i]; i ++) {
    if(words[i] == ',' || words[i] == '\n') {
      if(buf.size()) {
        this->words0.push_back(buf);
        this->ptrs0.push_back(std::vector<int>());
      }
      buf = std::string();
      if(words[i] == ',')
        flag = true;
      else
        flag = false;
      continue;
    } else if(words[i] == ' ' || words[i] == '\t')
      continue;
    if(!flag)
      buf += words[i];
  }
  this->nthresh = nthresh;
  this->Nthresh = Nthresh;
  this->Midx    = 1;
  return;
}

template <typename T, typename U> int corpus<T,U>::isin(const std::string& key) {
  return 0;
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

template <typename T, typename U> const void corpus<T,U>::compute(const U* input) {
  std::cerr << "Getting word pointers." << std::endl;
  getWordPtrs(input);
  std::cerr << "Corpus..." << std::endl;
  corpusEach();
  return;
}

template <typename T, typename U> const std::vector<std::string>& corpus<T,U>::getWords() const {
  return words;
}

template <typename T, typename U> const Eigen::Matrix<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic>& corpus<T,U>::getCorpus() const {
  return corpust;
}

template <typename T, typename U> void corpus<T,U>::getWordPtrs(const U* input) {
  std::sort(words0.begin(), words0.end());
  std::string work;
  for(int i = 0; input[i]; i ++) {
    work += input[i];
    if(std::binary_search(words0.begin(), words0.end(), work, partial_compare)) {
      std::pair<vsitr, vsitr> p = std::equal_range(words0.begin(), words0.end(), work, partial_compare);
      bool match = false;
      for(vsitr pp = p.first; pp != p.second; ++ pp) {
        if(work == *pp) {
          ptrs0[std::distance(words0.begin(), pp)].push_back(i);
          match = false;
          Midx  = T(i);
          break;
        } else if(partial_compare_equal(work, *pp))
          match = true;
      }
      if(match)
        continue;
    }
    i -= work.size() - 1;
    work = std::string();
  }
  for(vsitr itr = words0.begin(); itr != words0.end(); ++ itr) {
    const int idx = std::distance(words0.begin(), itr);
    if(ptrs0[idx].size()) {
      words.push_back(*itr);
      ptrs.push_back(ptrs0[idx]);
    }
  }
  std::cerr << words.size() << " words used." << std::endl;
  return;
}

template <typename T, typename U> void corpus<T,U>::corpusEach() {
  corpust = Tensor(words.size(), words.size());
  for(int i = 0; i < words.size(); i ++) {
    std::cerr << "Corpushing row : " << i << std::endl;
    if(!ptrs[i].size())
      continue;
    for(int j = 0; j < words.size(); j ++) {
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
        for(viitr itr = ptrs[k].begin(); itr != ptrs[k].end(); ++ itr) {
          while(ctru < ptrs[i].size() && ptrs[i][ctru] < *itr) ctru ++;
          if(ctru >= ptrs[i].size())
            continue;
          while(ctrv < ptrs[j].size() && ptrs[j][ctrv] < *itr) ctrv ++;
          if(ptrs[j][ctrv] < *itr)
            break;
          T work(0);
          const T buf0(log(std::abs(*itr - ptrs[i][ctru])));
          const T buf1(log(std::abs(*itr - ptrs[j][ctrv])));
          if(nthresh < buf0 && buf0 < Nthresh &&
             nthresh < buf1 && buf1 < Nthresh)
            work += buf0 * buf0 + buf1 * buf1;
          // XXX fixme:
          corpust(i, j)[k] += sqrt(work) / Midx;
        }
      }
    }
  }
  return;
}

#define _CORPUS_
#endif


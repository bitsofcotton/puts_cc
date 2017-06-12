#if !defined(_LWORD_)

#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>

using namespace std;

template <typename T> class word_t {
public:
  T*  str;
  int count;
};

template <typename T> class gram_t {
public:
  T*          str;
  vector<int> ptr0;
  vector<int> ptr1;
  gram_t() {
    this->str = 0;
  }
  gram_t(const gram_t& x) {
    this->str  = x.str;
    this->ptr0 = x.ptr0;
    this->ptr1 = x.ptr1;
  }
};

template <typename T> static int cmpwrap(const gram_t<T>& x0, const gram_t<T>& x1);
template <typename T> static int wordcmp(const word_t<T>& x0, const word_t<T>& x1);

template <typename T> class lword {
public:
  lword();
  ~lword();
  void init(const int& loop, const int& mthresh, const int& Mthresh);
  
  const std::vector<word_t<T> >& compute(const T* input);
private:
  std::vector<T>                        dict0;
  std::vector<std::vector<gram_t<T> > > dicts;
  std::vector<word_t<T> >               words;
  
  int loop;
  int mthresh;
  int Mthresh;
  
  void makeBigram(const T* input);
  void constructNwords();
  void bondLast();
  void makeWords();
  
  bool       isin(const T* key, const vector<gram_t<T> >& dict);
  gram_t<T>& find(const T* key, vector<gram_t<T> >& dict);
  void       assign(const T* key, vector<gram_t<T> >& dict, const gram_t<T>& val);
};

template <typename T> lword<T>::lword() {
  init(120, 2, 2);
}

template <typename T> lword<T>::~lword() {
  for(auto itr = dicts.begin(); itr != dicts.end(); ++ itr)
    for(auto itr2 = itr->begin(); itr2 != itr->end(); ++ itr2)
      if(itr2->str) {
        delete[] itr2->str;
        itr2->str = 0;
      }
  // Already freed.
/*
  for(auto itr = words.begin(); itr != words.end(); ++ itr)
    if(itr->str) {
      delete[] itr->str;
      itr->str = 0;
    }
*/
}

template <typename T> void lword<T>::init(const int& loop, const int& mthresh, const int& Mthresh) {
  this->loop    = loop;
  this->mthresh = mthresh;
  this->Mthresh = Mthresh;
}

template <typename T> int cmpwrap(const gram_t<T>& x0, const gram_t<T>& x1) {
  return std::strcmp(x0.str, x1.str) < 0;
}

template <typename T> int wordcmp(const word_t<T>& x0, const word_t<T>& x1) {
  return x0.count < x1.count || std::strcmp(x0.str, x1.str) < 0;
}

template <typename T> bool lword<T>::isin(const T* key, const vector<gram_t<T> >& dict) {
  gram_t<T> key0;
  key0.str = const_cast<T*>(key);
  return std::binary_search(dict.begin(), dict.end(), key0, cmpwrap<T>);
}

template <typename T> gram_t<T>& lword<T>::find(const T* key, vector<gram_t<T> >& dict) {
  if(!isin(key, dict)) {
    static gram_t<T> dummy;
    return dummy;
  }
  gram_t<T> key0;
  key0.str = const_cast<T*>(key);
  auto p = std::equal_range(dict.begin(), dict.end(), key0, cmpwrap<T>);
  return *(p.first);
}

template <typename T> void lword<T>::assign(const T* key, vector<gram_t<T> >& dict, const gram_t<T>& val) {
  gram_t<T> key0;
  key0.str = const_cast<T*>(key);
  auto p = std::equal_range(dict.begin(), dict.end(), key0, cmpwrap<T>);
  p.first->ptr0 = val.ptr0;
  p.first->ptr1 = val.ptr1;
  return;
}

template <typename T> const std::vector<word_t<T> >& lword<T>::compute(const T* input) {
  cerr << "bi-gramming..." << endl;
  makeBigram(input);
  cerr << "constructing longest words..." << endl;
  constructNwords();
  cerr << "gathering unknown words..." << endl;
  bondLast();
  cerr << "making word tables..." << endl;
  makeWords();
  return words;
}

template <typename T> void lword<T>::makeBigram(const T* input) {
  std::map<T, T> map0;
  std::map<std::string, vector<int> > map;
  for(int i = 1; input[i]; i ++) {
    T work[3];
    work[0] = input[i - 1];
    work[1] = input[i];
    work[2] = T(0);
    map[std::string(work)].push_back(i);
    map0[input[i]] = input[i];
  }
  for(auto itr = map0.begin(); itr != map0.end(); ++ itr)
    dict0.push_back(itr->first);
  std::sort(dict0.begin(), dict0.end());
  std::vector<gram_t<T> > buf;
  for(auto itr = map.begin(); itr != map.end(); ++ itr) {
    gram_t<T> work;
    work.str = new T[3];
    for(int i = 0; i < 2; i ++)
      work.str[i] = itr->first[i];
    work.str[2] = '\0';
    for(auto itr2 = itr->second.begin(); itr2 != itr->second.end(); ++ itr2) {
      work.ptr0.push_back(*itr2 - 1);
      work.ptr1.push_back(*itr2);
    }
    buf.push_back(work);
  }
  dicts.push_back(buf);
  std::sort(dicts[0].begin(), dicts[0].end(), cmpwrap<T>);
}

template <typename T> void lword<T>::constructNwords() {
  int i;
  for(i = 1; i < loop; i ++) {
    std::cerr << "constructing " << i + 2 << " table..." << std::endl;
    std::map<std::string, vector<int> > map0;
    std::map<std::string, vector<int> > map1;
    std::map<std::string, vector<int> > dmap;
    for(auto itr = dicts[i - 1].begin(); itr != dicts[i - 1].end(); ++ itr) {
      const T* key = itr->str;
      if(!isin(key, dicts[i - 1])) {
        std::cerr << "Something odd has occurred." << std::endl;
        continue;
      }
      gram_t<T> idxkey(find(key, dicts[i - 1]));
      for(auto itr2 = dict0.begin(); itr2 != dict0.end(); ++ itr2) {
        T key2[i + 2];
        for(int j = 1; j < i + 1; j ++)
          key2[j - 1] = key[j];
        key2[i]       = *itr2;
        key2[i + 1]   = '\0';
        if(!isin(key2, dicts[i - 1]))
          continue;
        T workkey[i + 3];
        for(int j = 0; j < i + 1; j ++)
          workkey[j] = key[j];
        workkey[i + 1] = *itr2;
        workkey[i + 2] = '\0';
        vector<int> idxwork[4];
        gram_t<T>   idxkey2(find(key2, dicts[i - 1]));
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
        const int diff(std::min(idxkey.ptr0.size(), idxkey2.ptr0.size()) - idxwork[0].size());
        const int diffM(std::max(idxkey.ptr0.size(), idxkey2.ptr0.size()) - idxwork[0].size());
        if(0 < diff && diffM < mthresh)
          continue;
        for(int i = 0; i < idxwork[0].size(); i ++) {
          map0[workkey].push_back(idxwork[0][i]);
          map1[workkey].push_back(idxwork[1][i]);
          dmap[key].push_back(idxwork[2][i]);
          dmap[key2].push_back(idxwork[3][i]);
        }
      }
    }
    // delete this stage garbage.
    for(auto itr = dmap.begin(); itr != dmap.end(); ++ itr) {
      std::sort(itr->second.begin(), itr->second.end());
      if(!isin(itr->first.c_str(), dicts[i - 1])) {
        std::cerr << "Something odd has occured." << std::endl;
        continue;
      }
      gram_t<T> before(find(itr->first.c_str(), dicts[i - 1]));
      gram_t<T> after;
      after.str = before.str;
      for(int j = 0; j < before.ptr0.size(); j ++) {
        bool flag = false;
        for(auto itr2 = dmap[itr->first].begin(); itr2 != dmap[itr->first].end(); ++ itr2)
          if(*itr2 == j) {
            flag = true;
            break;
          }
        if(flag)
          continue;
        after.ptr0.push_back(before.ptr0[j]);
        after.ptr1.push_back(before.ptr1[j]);
      }
      assign(itr->first.c_str(), dicts[i - 1], after);
    }
    // construct next stage.
    vector<gram_t<T> > workv;
    for(auto itr = map0.begin(); itr != map0.end(); ++ itr) {
      gram_t<T>   work;
      work.str = new T[i + 3];
      for(int j = 0; j < i + 2; j ++)
        work.str[j] = itr->first[j];
      work.str[i + 2] = '\0';
      const vector<int>& ptr0orig(itr->second);
      const vector<int>& ptr1orig(map1[itr->first]);
      for(int j = 0; j < ptr0orig.size(); j ++) {
        work.ptr0.push_back(ptr0orig[j]);
        work.ptr1.push_back(ptr1orig[j]);
      }
      workv.push_back(work);
    }
    if(workv.size() < 1)
      break;
    std::sort(workv.begin(), workv.end(), cmpwrap<T>);
    dicts.push_back(workv);
  }
  return;
}

template <typename T> void lword<T>::bondLast() {
  bool loopf = true;
  int  count = 0;
  while(loopf) {
    loopf = false;
    count ++;
    for(int i = 0; i < dicts.size(); i ++) {
      std::cerr << "bondLast: " << i << "/" << dicts.size() << " @ " << count << std::endl;
      for(auto itr = dicts[i].begin(); itr != dicts[i].end(); ++ itr)
        for(int j = 0; j < dicts.size(); j ++) {
          if(i + j + 1 >= loop)
            continue;
          for(auto itr2 = dicts[j].begin();  itr2 != dicts[j].end(); ++ itr2) {
            std::vector<int> ptr0;
            std::vector<int> ptr1;
            if(!itr->str || itr2->str)
              continue;
            if(itr->str[i + 1] == itr2->str[0]) {
              for(int ii = 0; ii < itr->ptr0.size(); ii ++)
                for(int jj = 0; jj < itr2->ptr0.size(); jj ++)
                  if(itr->ptr1[ii] == itr2->ptr0[jj]) {
                    ptr0.push_back(itr->ptr0[ii]);
                    ptr1.push_back(itr2->ptr1[jj]);
                  }
              const int diff(std::min(itr->ptr0.size(), itr2->ptr0.size()) - ptr0.size());
              const int diffM(std::max(itr->ptr0.size(), itr2->ptr0.size()) - ptr0.size());
              if(0 < diff && diffM < mthresh)
                continue;
              if(ptr0.size() < Mthresh)
                continue;
              // insert longer word.
              T key2[i + j + 3];
              for(int ii = 0; ii <= i; ii ++)
                key2[ii] = itr->str[ii];
              for(int jj = 0; jj <= j; jj ++)
                key2[i + 1 + jj] = itr2->str[jj];
              key2[i + j + 2] = '\0';
              if(!isin(key2, dicts[i + j])) {
                gram_t<T> work;
                work.str = new T[i + j + 1];
                for(int jj = 0; jj < i + j + 1; jj ++)
                  work.str[jj] = key2[jj];
                work.ptr0 = ptr0;
                work.ptr1 = ptr1;
                dicts[i + j].push_back(work);
              } else {
                gram_t<T> work(find(key2, dicts[i + j + 1]));
                int k = 0;
                for(int jj = 0; jj < work.ptr0.size(); jj ++)
                  for(; k < ptr0.size() && (jj == work.ptr0.size() - 1 || ptr0[k] <= work.ptr0[jj]); k ++)
                    if(work.ptr0[jj] <= ptr0[k] && (jj == work.ptr0.size() - 1 || ptr0[k] <= work.ptr0[jj + 1])) {
                      work.ptr0.insert(work.ptr0.begin() + jj, ptr0[k]);
                      work.ptr1.insert(work.ptr1.begin() + jj, ptr1[k]);
                    }
                assign(key2, dicts[i + j + 1], work);
              }
              // delete.
              std::vector<int> i0ptr0;
              std::vector<int> i0ptr1;
              std::vector<int> i1ptr0;
              std::vector<int> i1ptr1;
              for(int ii = 0; ii < itr->ptr0.size(); ii ++) {
                bool flag = false;
                for(auto itr3 = ptr0.begin(); itr3 != ptr0.end(); ++ itr3)
                  if(itr->ptr0[ii] == *itr3) {
                    flag = true;
                    break;
                  }
                if(!flag) {
                  i0ptr0.push_back(itr->ptr0[ii]);
                  i0ptr1.push_back(itr->ptr1[ii]);
                }
              }
              itr->ptr0 = i0ptr0;
              itr->ptr1 = i0ptr1;
              for(int ii = 0; ii < itr2->ptr0.size(); ii ++) {
                bool flag = false;
                for(auto itr3 = ptr1.begin(); itr3 != ptr1.end(); ++ itr3)
                  if(itr2->ptr1[ii] == *itr3) {
                    flag = true;
                    break;
                  }
                if(!flag) {
                  i1ptr0.push_back(itr2->ptr0[ii]);
                  i1ptr1.push_back(itr2->ptr1[ii]);
                }
              }
              itr2->ptr0 = i1ptr0;
              itr2->ptr1 = i1ptr1;
              loopf = true;
            }
          }
        }
    }
  }
  return;
}

template <typename T> void lword<T>::makeWords() {
  // for(auto itr = dicts.begin(); itr != dicts.end(); ++ itr) {
  for(int idx = 1; idx < dicts.size(); idx ++) {
    const vector<gram_t<T> >* itr(&dicts[idx]);
    for(auto itr2 = itr->begin(); itr2 != itr->end(); ++ itr2) {
      word_t<T> work;
      work.str   = itr2->str;
      work.count = itr2->ptr0.size();
      if(!work.count)
        continue;
      words.push_back(work);
    }
  }
  std::sort(words.begin(), words.end(), wordcmp<T>);
  return;
}

#define _LWORD_
#endif


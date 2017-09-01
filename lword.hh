#if !defined(_LWORD_)

#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <cstring>

using std::strcmp;
using std::vector;
using std::map;
using std::string;
using std::pair;
using std::binary_search;
using std::equal_range;
using std::sort;
using std::cerr;
using std::endl;
using std::max;
using std::min;

template <typename T> class word_t {
public:
  T*  str;
  int count;
  word_t() {
    str   = NULL;
    count = 0;
  }
  word_t(const word_t& other) {
    this->str   = other.str;
    this->count = other.count;
  }
  ~word_t() {
    // frees in another place.
    ;
  }
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

template <typename T> int cmpwrap(const gram_t<T>& x0, const gram_t<T>& x1) {
  return strcmp(x0.str, x1.str) < 0;
}

template <typename T> int cmpwrap2(const word_t<T>& x0, const word_t<T>& x1) {
  // gcc4.9 is bugly...
  return (x0.count < x1.count) || (x0.count == x1.count && strcmp(x0.str, x1.str) < 0);
  // return strcmp(x0.str, x1.str) < 0;
}


template <typename T> class lword {
public:
  lword();
  ~lword();
  void init(const int& loop, const int& mthresh, const int& Mthresh);
  
  const vector<word_t<T> >& compute(const T* input);
private:
  vector<T>                   dict0;
  vector<vector<gram_t<T> > > dicts;
  vector<word_t<T> >          words;
  typedef typename vector<T>::iterator          vTitr;
  typedef typename vector<gram_t<T> >::iterator vgTitr;
  typedef typename vector<word_t<T> >::iterator vwTitr;
  typedef vector<int>::iterator                 viitr;
  typedef typename vector<vector<gram_t<T> > >::iterator vvgTitr;
  typedef typename map<T, T>::iterator   mTTitr;
  typedef typename map<T, int>::iterator mTiitr;
  typedef map<string, vector<int> >::iterator msviitr;
  
  int loop;
  int mthresh;
  int Mthresh;
  
  void makeBigram(const T* input);
  void constructNwords();
  void bondLast();
  void balanceBonded();
  void makeWords();
  
  bool       isin(const T* key, const vector<gram_t<T> >& dict);
  gram_t<T>& find(const T* key, vector<gram_t<T> >& dict);
  void       assign(const T* key, vector<gram_t<T> >& dict, const gram_t<T>& val);
};

template <typename T> lword<T>::lword() {
  init(120, 2, 2);
}

template <typename T> lword<T>::~lword() {
  for(vvgTitr itr = dicts.begin(); itr != dicts.end(); ++ itr)
    for(vgTitr itr2 = itr->begin(); itr2 != itr->end(); ++ itr2)
      if(itr2->str) {
        delete[] itr2->str;
        itr2->str = 0;
      }
  // Already freed.
/*
  for(vsitr itr = words.begin(); itr != words.end(); ++ itr)
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

template <typename T> bool lword<T>::isin(const T* key, const vector<gram_t<T> >& dict) {
  gram_t<T> key0;
  key0.str = const_cast<T*>(key);
  return binary_search(dict.begin(), dict.end(), key0, cmpwrap<T>);
}

template <typename T> gram_t<T>& lword<T>::find(const T* key, vector<gram_t<T> >& dict) {
  if(!isin(key, dict)) {
    static gram_t<T> dummy;
    return dummy;
  }
  gram_t<T> key0;
  key0.str = const_cast<T*>(key);
  pair<vgTitr, vgTitr> p = equal_range(dict.begin(), dict.end(), key0, cmpwrap<T>);
  return *(p.first);
}

template <typename T> void lword<T>::assign(const T* key, vector<gram_t<T> >& dict, const gram_t<T>& val) {
  gram_t<T> key0;
  key0.str = const_cast<T*>(key);
  pair<vgTitr, vgTitr> p = equal_range(dict.begin(), dict.end(), key0, cmpwrap<T>);
  p.first->ptr0 = val.ptr0;
  p.first->ptr1 = val.ptr1;
  return;
}

template <typename T> const vector<word_t<T> >& lword<T>::compute(const T* input) {
  cerr << "bi-gramming..." << endl;
  makeBigram(input);
  cerr << "constructing longest words..." << endl;
  constructNwords();
  cerr << "gathering unknown words..." << endl;
  bondLast();
  cerr << "balancing bonded word tables..." << endl;
  balanceBonded();
  cerr << "making word tables..." << endl;
  makeWords();
  return words;
}

template <typename T> void lword<T>::makeBigram(const T* input) {
  map<T, T> map0;
  map<string, vector<int> > map;
  for(int i = 1; input[i]; i ++) {
    T work[3];
    work[0] = input[i - 1];
    work[1] = input[i];
    work[2] = T(0);
    map[string(work)].push_back(i);
    map0[input[i]] = input[i];
  }
  for(auto itr = map0.begin(); itr != map0.end(); ++ itr)
    dict0.push_back(itr->first);
  sort(dict0.begin(), dict0.end());
  vector<gram_t<T> > buf;
  for(auto itr = map.begin(); itr != map.end(); ++ itr) {
    gram_t<T> work;
    work.str = new T[3];
    for(int i = 0; i < 2; i ++)
      work.str[i] = itr->first[i];
    work.str[2] = '\0';
    for(viitr itr2 = itr->second.begin(); itr2 != itr->second.end(); ++ itr2) {
      work.ptr0.push_back(*itr2 - 1);
      work.ptr1.push_back(*itr2);
    }
    buf.push_back(work);
  }
  dicts.push_back(buf);
  sort(dicts[0].begin(), dicts[0].end(), cmpwrap<T>);
}

template <typename T> void lword<T>::constructNwords() {
  int i;
  for(i = 1; i < loop; i ++) {
    cerr << "constructing " << i + 2 << " table..." << endl;
    map<string, vector<int> > map0;
    map<string, vector<int> > map1;
    map<string, vector<int> > dmap;
    for(vgTitr itr = dicts[i - 1].begin(); itr != dicts[i - 1].end(); ++ itr) {
      const T* key = itr->str;
      if(!isin(key, dicts[i - 1])) {
        cerr << "Something odd has occurred." << endl;
        continue;
      }
      gram_t<T> idxkey(find(key, dicts[i - 1]));
      for(vTitr itr2 = dict0.begin(); itr2 != dict0.end(); ++ itr2) {
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
        const int diff(min(idxkey.ptr0.size(), idxkey2.ptr0.size()) - idxwork[0].size());
        const int diffM(max(idxkey.ptr0.size(), idxkey2.ptr0.size()) - idxwork[0].size());
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
    for(msviitr itr = dmap.begin(); itr != dmap.end(); ++ itr) {
      sort(itr->second.begin(), itr->second.end());
      if(!isin(itr->first.c_str(), dicts[i - 1])) {
        cerr << "Something odd has occured." << endl;
        continue;
      }
      gram_t<T> before(find(itr->first.c_str(), dicts[i - 1]));
      gram_t<T> after;
      after.str = before.str;
      for(int j = 0; j < before.ptr0.size(); j ++) {
        bool flag = false;
        for(viitr itr2 = dmap[itr->first].begin(); itr2 != dmap[itr->first].end(); ++ itr2)
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
    for(msviitr itr = map0.begin(); itr != map0.end(); ++ itr) {
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
    sort(workv.begin(), workv.end(), cmpwrap<T>);
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
      cerr << "bondLast: " << i << "/" << dicts.size() << " @ " << count << endl;
      for(vgTitr itr = dicts[i].begin(); itr != dicts[i].end(); ++ itr)
        for(int j = 0; j < dicts.size(); j ++) {
          if(i + j + 1 >= loop)
            continue;
          for(vgTitr itr2 = dicts[j].begin();  itr2 != dicts[j].end(); ++ itr2) {
            vector<int> ptr0;
            vector<int> ptr1;
            if(!itr->str || itr2->str)
              continue;
            if(itr->str[i + 1] == itr2->str[0]) {
              for(int ii = 0; ii < itr->ptr0.size(); ii ++)
                for(int jj = 0; jj < itr2->ptr0.size(); jj ++)
                  if(itr->ptr1[ii] == itr2->ptr0[jj]) {
                    ptr0.push_back(itr->ptr0[ii]);
                    ptr1.push_back(itr2->ptr1[jj]);
                  }
              const int diff(min(itr->ptr0.size(), itr2->ptr0.size()) - ptr0.size());
              const int diffM(max(itr->ptr0.size(), itr2->ptr0.size()) - ptr0.size());
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
              vector<int> i0ptr0;
              vector<int> i0ptr1;
              vector<int> i1ptr0;
              vector<int> i1ptr1;
              for(int ii = 0; ii < itr->ptr0.size(); ii ++) {
                bool flag = false;
                for(viitr itr3 = ptr0.begin(); itr3 != ptr0.end(); ++ itr3)
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
                for(viitr itr3 = ptr1.begin(); itr3 != ptr1.end(); ++ itr3)
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

template <typename T> void lword<T>::balanceBonded() {
  int  count = 0;
  bool loopf = true;
  while(loopf) {
    loopf = false;
    count ++;
    for(int i = 0; i < dicts.size(); i ++) {
      cerr << "balance: " << i << "/" << dicts.size() << " @ " << count << endl;
      map<T, int> wstat;
      map<T, int> wstat2;
      for(vgTitr itr = dicts[i].begin(); itr != dicts[i].end(); ++ itr) {
        for(int j = 0; j < dicts.size(); j ++)
          for(vgTitr itr2 = dicts[j].begin(); itr2 != dicts[j].end(); ++ itr2) {
            if(itr->str[i + 1] == itr2->str[0]) {
              wstat[itr2->str[1]] = 0;
              for(viitr litr = itr->ptr1.begin(); litr != itr->ptr1.end(); ++ litr)
                for(viitr litr2 = itr2->ptr0.begin(); litr2 != itr2->ptr0.end(); ++ litr2)
                  if(*litr == *litr2)
                    wstat[itr2->str[1]] ++;
            }
            if(itr->str[0] == itr2->str[j + 1]) {
              wstat2[itr->str[1]] = 0;
              for(viitr litr = itr->ptr0.begin(); litr != itr->ptr0.end(); ++ litr)
                for(viitr litr2 = itr2->ptr1.begin(); litr2 != itr2->ptr1.end(); ++ litr2)
                  if(*litr == *litr2)
                    wstat2[itr2->str[j]] ++;
            }
          }
        int stat(0), stat2(0);
        T   str(0),  str2(0);
        for(mTiitr itr2 = wstat.begin(); itr2 != wstat.end(); ++ itr2)
          if(stat < itr2->second) {
            str  = itr2->first;
            stat = itr2->second;
          }
        for(mTiitr itr2 = wstat2.begin(); itr2 != wstat2.end(); ++ itr2)
          if(stat < itr2->second) {
            str2  = itr2->first;
            stat2 = itr2->second;
          }
        for(int j = 0; j < dicts.size(); j ++)
          for(vgTitr itr2 = dicts[j].begin(); itr2 != dicts[j].end(); ++ itr2) {
            if(j >= 1 && itr->ptr0.size() <= stat && itr->str[i + 1] == itr2->str[0] && itr2->str[1] == str) {
              T key[j + 2];
              for(int k = 0; k < j + 1; k ++)
                key[k] = itr2->str[k + 1];
              key[j + 1] = '\0';
              gram_t<T> work;
              if(isin(key, dicts[j - 1]))
                work = find(key, dicts[j - 1]);
              else {
                work.str = new T[j + 2];
                for(int k = 0; k < j + 2; k ++)
                  work.str[k] = key[k];
              }
              if(itr->ptr1.size() > 0) {
                vector<int> rptr;
                for(viitr litr = itr->ptr1.begin(); litr != itr->ptr1.end(); ++ litr)
                  for(int iitr2 = 0; iitr2 < itr2->ptr0.size(); iitr2 ++)
                    if(*litr == itr2->ptr0[iitr2]) {
                      rptr.push_back(*litr);
                      work.ptr0.push_back(itr2->ptr0[iitr2] + 1);
                      work.ptr1.push_back(itr2->ptr1[iitr2]);
                      itr2->ptr0.erase(itr2->ptr0.begin() + iitr2);
                      itr2->ptr1.erase(itr2->ptr1.begin() + iitr2);
                      loopf = true;
                      break;
                    }
                for(int k = 0, l = 0; k < rptr.size(); k ++)
                  if(itr->ptr1[l ++] == rptr[k])
                    itr->ptr1[l - 1] ++;
              }
              assign(work.str, dicts[j - 1], work);
            }
            if(j >= 1 && itr->ptr0.size() <= stat2 && itr->str[0] == itr2->str[j + 1] && itr2->str[j] == str2) {
              T key[j + 2];
              for(int k = 0; k < j + 1; k ++)
                key[k] = itr->str[k];
              key[j + 1] = '\0';
              gram_t<T> work;
              if(isin(key, dicts[j - 1]))
                work = find(key, dicts[j - 1]);
              else {
                work.str = new T[j + 2];
                for(int k = 0; k < j + 2; k ++)
                  work.str[k] = key[k];
              }
              if(itr2->ptr1.size() > 0) {
                vector<int> rptr;
                for(viitr litr = itr2->ptr1.begin(); litr != itr2->ptr1.end(); ++ litr)
                  for(int iitr2 = 0; iitr2 < itr->ptr0.size(); iitr2 ++)
                    if(*litr == itr->ptr0[iitr2]) {
                      rptr.push_back(*litr);
                      work.ptr0.push_back(itr->ptr0[iitr2] + 1);
                      work.ptr1.push_back(itr->ptr1[iitr2]);
                      itr->ptr0.erase(itr->ptr0.begin() + iitr2);
                      itr->ptr1.erase(itr->ptr1.begin() + iitr2);
                      loopf = true;
                      break;
                    }
                for(int k = 0, l = 0; k < rptr.size(); k ++)
                  if(itr2->ptr1[l ++] == rptr[k])
                    itr->ptr0[l - 1] --;
              }
              assign(work.str, dicts[j - 1], work);
            }
          }
        if(itr->ptr0.size() <= stat && itr->ptr0.size() <= stat2 && i + 2 < dicts.size()) {
          T bstr[i + 4];
          bstr[0] = str2;
          for(int k = 0; k < i + 1; k ++)
            bstr[k + 1] = itr->str[k];
          bstr[i + 2] = str;
          bstr[i + 3] = '\0';
          if(isin(bstr, dicts[i + 2])) {
            gram_t<T>& work(find(bstr, dicts[i + 2]));
            for(int i = 0; i < itr->ptr0.size(); i ++) {
              work.ptr0.push_back(itr->ptr0[i]);
              work.ptr1.push_back(itr->ptr1[i]);
            }
          } else {
            gram_t<T> work(*itr);
            work.str = new T[i + 4];
            for(int k = 0; k < i + 4; k ++)
              work.str[k] = bstr[k];
            dicts[i + 2].push_back(work);
          }
          delete[] itr->str;
          dicts[i].erase(itr);
          break;
        } else if((itr->ptr0.size() <= stat || itr->ptr0.size() <= stat2) && i + 1 < dicts.size()) {
          T bstr[i + 3];
          if(itr->ptr0.size() <= stat) {
            for(int k = 0; k < i + 1; k ++)
              bstr[k] = itr->str[k];
            bstr[i + 1] = str;
            bstr[i + 2] = '\0';
          } else {
            bstr[0] = str2;
            for(int k = 0; k < i + 1; k ++)
              bstr[k + 1] = itr->str[k];
            bstr[i + 2] = '\0';
          }
          if(isin(bstr, dicts[i + 1])) {
            gram_t<T>& work(find(bstr, dicts[i + 1]));
            for(int i = 0; i < itr->ptr0.size(); i ++) {
              work.ptr0.push_back(itr->ptr0[i]);
              work.ptr1.push_back(itr->ptr1[i]);
            }
          } else {
            gram_t<T> work(*itr);
            work.str = new T[i + 3];
            for(int k = 0; k < i + 3; k ++)
              work.str[k] = bstr[k];
            dicts[i + 1].push_back(work);
          }
          delete[] itr->str;
          dicts[i].erase(itr);
          break;
        }
      }
    }
  }
  return;
}

template <typename T> void lword<T>::makeWords() {
  for(vvgTitr itr = dicts.begin(); itr != dicts.end(); ++ itr) {
    for(vgTitr itr2 = itr->begin(); itr2 != itr->end(); ++ itr2) {
      word_t<T> work;
      work.str   = itr2->str;
      work.count = itr2->ptr0.size();
      if(!work.count)
        continue;
      words.push_back(work);
    }
  }
  sort(words.begin(), words.end(), cmpwrap2<T>);
  return;
}

#define _LWORD_
#endif


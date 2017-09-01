#if !defined(_FILE2EIGEN_)

#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <sstream>

using std::istream;
using std::ostream;
using std::vector;
using std::endl;

template <typename T> ostream& operator << (ostream& os, const vector<T>& p) {
  os << p.size() << endl;
  for(int i = 0; i < p.size(); i ++)
    os << p[i] << endl;
  return os;
}

template <typename T> ostream& operator << (ostream& os, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& p) {
  os << p.rows() << endl << p.cols() << endl;
  for(int i = 0; i < p.rows(); i ++) {
    for(int j = 0; j < p.cols(); j ++)
      os << p(i, j) << " ";
    os << endl;
  }
  return os;
}

template <typename T> ostream& operator << (ostream& os, const Eigen::Matrix<T, Eigen::Dynamic, 1>& p) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> work(p.size(), 1);
  work.col(0) = p;
  os << work;
  return os;
}

template <typename T> istream& operator >> (istream& is, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& work) {
  int rows, cols;
  is >> rows;
  is >> cols;
  work = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(rows, cols);
  for(int i = 0; i < work.rows(); i ++)
    for(int j = 0; j < work.cols(); j ++)
      is >> work(i, j);
  return is;
}

template <typename T> istream& operator >> (istream& is, Eigen::Matrix<T, Eigen::Dynamic, 1>& work) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> buf;
  is >> buf;
  work = buf.col(0);
  return is;
}

#define _FILE2EIGEN_
#endif


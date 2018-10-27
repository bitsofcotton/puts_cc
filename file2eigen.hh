/* BSD 3-Clause License:
 * Copyright (c) 2018, bitsofcotton.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(_FILE2EIGEN_)

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


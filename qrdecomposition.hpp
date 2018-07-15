#ifndef AUTODIFF_QRDECOMPOSITION_HPP_
#define AUTODIFF_QRDECOMPOSITION_HPP_


#include "mat33.hpp"


namespace autodiff {


template <typename T>
struct QRDecomposition {
  Mat33<T> q;
  Mat33<T> r;
};


template <typename T>
std::ostream& operator<<(std::ostream &stream,const QRDecomposition<T> &qr)
{
  stream << "q:\n";
  stream << qr.q << "\n";
  stream << "r:\n";
  stream << qr.r << "\n";
  return stream;
}


}


#endif /* AUTODIFF_QRDECOMPOSITION_HPP_ */

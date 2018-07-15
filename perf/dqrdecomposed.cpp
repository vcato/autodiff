#include "dqrdecomposed.hpp"

#include "../dualmat33.hpp"
#include "../qrdecomposedexpr.hpp"
#include "../evalandaddderiv.hpp"



using autodiff::Mat33;
using autodiff::QRDecomposition;
using autodiff::dual;


void
  dqrDecomposed(
    const Mat33<float> &a,Mat33<float> &da,
    const QRDecomposition<float> &dresult
  )
{
  evalAndAddDeriv(qrDecomposed(dual(a,da)),dresult);
}

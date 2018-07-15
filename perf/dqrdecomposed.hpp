#include "../mat33.hpp"
#include "../qrdecomposition.hpp"


extern void
  dqrDecomposed(
    const autodiff::Mat33<float> &a,autodiff::Mat33<float> &da,
    const autodiff::QRDecomposition<float> &dresult
  );

#include "../random.hpp"
#include "../mat33.hpp"
#include "../randommat33.hpp"
#include "../qrdecomposition.hpp"
#include "dqrdecomposed.hpp"


using autodiff::Mat33;
using autodiff::zeroMat33;
using autodiff::QRDecomposition;


static void benchmark()
{
  RandomEngine random_engine(/*seed*/1);

  Mat33<float> a = randomMat33<float>(random_engine);
  Mat33<float> dq = randomMat33<float>(random_engine);
  Mat33<float> dr = randomMat33<float>(random_engine);
  
  Mat33<float> da = zeroMat33<float>();

  for (int i=0; i!=100000; ++i) {
    dqrDecomposed(a,da,QRDecomposition<float>{dq,dr});
  }
}



int main()
{
  benchmark();
}

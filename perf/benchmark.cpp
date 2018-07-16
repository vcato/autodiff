#include "../random.hpp"
#include "../mat33.hpp"
#include "../randommat33.hpp"
#include "../qrdecomposition.hpp"
#include "dqrdecomposed.hpp"


using std::cerr;
using std::string;
using autodiff::Mat33;
using autodiff::zeroMat33;
using autodiff::QRDecomposition;


static QRDecomposition<float> qrDecomposed(const Mat33<float> &a)
{
  using T = float;
  // Make q*r == a
  // q is an orthogonal matrix
  // r is an upper triangular matrix

  // a = q*r
  // [a11,a12,a13]   [q11,q12,q13]   [r11,r12,r13]
  // [a21,a22,a23] = [q21,q22,q23] * [  0,r22,r23]
  // [a31,a32,a33]   [q31,a32,a33]   [  0,  0,r33]

  // q1 = <q11,q21,q31>
  // q2 = <q12,q22,q32>
  // q3 = <q13,q23,q33>

  // a1 = <a11,a21,a31>
  // a2 = <a12,a22,a32>
  // a3 = <a13,a23,a33>

  auto a1 = vec3(col(a,0));
  auto a2 = vec3(col(a,1));
  auto a3 = vec3(col(a,2));

  // a1 = q1*r11
  // a2 = q1*r12 + q2*r22
  // a3 = q1*r13 + q2*r23 + q3*r33

  // each ui is orthogonal to uj when i!=j
  //   ui = ai - sum{j if j!=i}(dot(qj,ai)*qj)
  // qi is the normalized version of ui
  //   qi = ui/mag(ui)

  // u1 = a1
  auto u1 = a1;
  // q1 = u1/mag(u1)
  // q1*mag(u1) = u1
  // u1 = q1*mag(u1)
  // a1 = q1*r11
  // r11 = mag(u1)
  auto r11 = mag(u1);

  // q1 = u1/mag(u1)
  // q1 = u1/r11
  auto q1 = u1/r11;

  // u2 = a2 - q1*dot(q1,a2)
  // u2 + q1*dot(q1,a2) = a2
  // a2 = u2 + q1*dot(q1,a2)
  // a2 = q1*dot(q1,a2) + u2
  // q2 = u2/mag(u2)
  // q2*mag(u2) = u2
  // u2 = q2*mag(u2)
  // a2 = u2 + q1*dot(q1,a2)
  // a2 = q1*dot(q1,a2) + q2*mag(u2)
  // a2 = q1*r12 + q2*r22
  // r12 = dot(q1,a2)
  auto r12 = dot(a2,q1);
  // u2 = a2 - q1*dot(q1,a2)
  // u2 = a2 - q1*r12
  auto u2 = a2 - q1*r12;
  // r22 = mag(u2)
  auto r22 = mag(u2);
  // q2 = u2/mag(u2)
  // q2 = u2/r22
  auto q2 = u2/r22;

  // u3 = a3 - q1*dot(q1,a3) - q2*dot(q2,a3)
  // u3 + q1*dot(q1,a3) + q2*dot(q2,a3) = a3
  // a3 = u3 + q1*dot(q1,a3) + q2*dot(q2,a3)
  // a3 = q1*dot(q1,a3) + q2*dot(q2,a3) + u3
  // a3 = q1*dot(q1,a3) + q2*dot(q2,a3) + q3*mag(u3)
  // a3 = q1*r13 + q2*r23 + q3*r33
  // r13 = dot(q1,a3)
  auto r13 = dot(q1,a3);
  // r23 = dot(q2,a3)
  auto r23 = dot(q2,a3);
  // r33 = mag(u3)
  // u3 = a3 - q1*dot(q1,a3) - q2*dot(q2,a3)
  // u3 = a3 - q1*r13 - q2*r23
  auto u3 = a3 - q1*r13 - q2*r23;
  auto r33 = mag(u3);
  // q3 = u3/mag(u3)
  // q3 = u3/r33
  auto q3 = u3/r33;

  T r_values[3][3] = {
    {r11,r12,r13},
    {  0,r22,r23},
    {  0,  0,r33}
  };

  auto q = columns(q1,q2,q3);
  Mat33<T> r(r_values);

  return {q,r};
}


static void verify()
{
  RandomEngine random_engine(/*seed*/1);

  Mat33<float> a = randomMat33<float>(random_engine);
  Mat33<float> dq = randomMat33<float>(random_engine);
  Mat33<float> dr = randomMat33<float>(random_engine);

  Mat33<float> da = zeroMat33<float>();
  Mat33<float> fda = zeroMat33<float>();

  dqrDecomposed(a,da,QRDecomposition<float>{dq,dr});

  float h = 1e-3;

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      double old = a[i][j];
      a[i][j] = old-h;
      QRDecomposition<float> qr1 = qrDecomposed(a);
      a[i][j] = old+h;
      QRDecomposition<float> qr2 = qrDecomposed(a);
      a[i][j] = old;

      double v = 0;

      for (int i2=0; i2!=3; ++i2) {
        for (int j2=0; j2!=3; ++j2) {
          v += (qr2.q[i2][j2] - qr1.q[i2][j2])/(2*h) * dq[i2][j2];
          v += (qr2.r[i2][j2] - qr1.r[i2][j2])/(2*h) * dr[i2][j2];
        }
      }

      fda[i][j] = v;
    }
  }

  // Determine how much q[i][j] changes due to a change in a, multiply
  // this by the gradient of q[i][j] to get the

  cerr << "da: " << da << "\n";
  cerr << "fda: " << fda << "\n";
}


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



int main(int argc,char **argv)
{
  if (argc!=2) {
    cerr << "Usage: benchmark verify\n";
    cerr << "Usage: benchmark run\n";
    return EXIT_FAILURE;
  }

  string op = argv[1];

  if (op=="verify") {
    verify();
    return EXIT_SUCCESS;
  }
  else if (op=="run") {
    benchmark();
    return EXIT_SUCCESS;
  }
  else {
    cerr << "Unknown operation: " << op << "\n";
    return EXIT_FAILURE;
  }
}

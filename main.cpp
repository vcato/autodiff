#include <iostream>
#include <cmath>
#include <cfloat>
#include "random.hpp"
#include "vec3.hpp"
#include "dualvec3.hpp"
#include "mat33.hpp"
#include "dualmat33.hpp"
#include "scalarexpr.hpp"
#include "evaluator.hpp"
#include "scalarexpr.hpp"
#include "vec3expr.hpp"
#include "mat33expr.hpp"
#include "evalandaddderiv.hpp"
#include "qrdecomposition.hpp"
#include "qrdecomposedexpr.hpp"

#define assertNear(actual,expected,tolerance) \
  (assertNearHelper(actual,expected,tolerance,__FILE__,__LINE__))


using std::cerr;
using std::ostream;
using std::max;
using autodiff::QRDecomposition;
using autodiff::Mat33;
using autodiff::FloatMat33;
using autodiff::FloatVec3;
using autodiff::mat33Identity;
using autodiff::zeroMat33;
using autodiff::evalAndAddDeriv;
using autodiff::dual;
using autodiff::ScalarAdd;
using autodiff::DualFloat;
using autodiff::DualVec3;
using autodiff::DualMat33;
using autodiff::ColRef;
using autodiff::vec3;
using autodiff::rotX;
using autodiff::Vec3ExprVar;
using autodiff::Mat33ExprVar;
using autodiff::Evaluator;




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


static float weightedSum(const FloatVec3 &v,const FloatVec3 &w)
{
  float x = v.x()*w.x();
  float y = v.y()*w.y();
  float z = v.z()*w.z();

  return x + y + z;
}


static float weightedSum(const FloatMat33 &m,const FloatMat33 &w)
{
  float result = 0;

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      result += m[i][j]*w[i][j];
    }
  }

  return result;
}


static float max(float a,float b,float c)
{
  return std::max(std::max(a,b),c);
}


static float differenceBetween(float a,float b)
{
  return fabsf(a-b);
}


static float differenceBetween(const FloatVec3 &a,const FloatVec3 &b)
{
  float dx = differenceBetween(a.x(), b.x());
  float dy = differenceBetween(a.y(), b.y());
  float dz = differenceBetween(a.z(), b.z());

  return max(dx,dy,dz);
}


static float differenceBetween(const FloatMat33 &a,const FloatMat33 &b)
{
  float max_d = -FLT_MAX;

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      float d = differenceBetween(a.values[i][j],b.values[i][j]);
      max_d = max(max_d,d);
    }
  }

  return max_d;
}


template <typename T>
static T
  differenceBetween(const QRDecomposition<T> &a,const QRDecomposition<T> &b)
{
  return max(differenceBetween(a.q,b.q),differenceBetween(a.r,b.r));
}



static FloatVec3 randomVec3(RandomEngine &random_engine)
{
  float x = randomFloat(-1,1,random_engine);
  float y = randomFloat(-1,1,random_engine);
  float z = randomFloat(-1,1,random_engine);

  return FloatVec3{x,y,z};
}


static FloatMat33 randomMat33(RandomEngine &random_engine)
{
  FloatVec3 x = randomVec3(random_engine);
  FloatVec3 y = randomVec3(random_engine);
  FloatVec3 z = randomVec3(random_engine);

  float values[3][3] = {
    {x.x(), x.y(), x.z()},
    {y.x(), y.y(), y.z()},
    {z.x(), z.y(), z.z()},
  };

  return FloatMat33(values);
}


template <typename Function>
static float finiteDeriv(Function f,float &v)
{
  float h = 1e-3;
  float old_value = v;
  v = old_value - h;
  float value1 = f();
  v = old_value + h;
  float value2 = f();
  v = old_value;
  return (value2-value1)/(2*h);
}


template <typename Function>
static FloatVec3 finiteDeriv(Function f,FloatVec3 &v)
{
  float x = finiteDeriv(f,v.x());
  float y = finiteDeriv(f,v.y());
  float z = finiteDeriv(f,v.z());

  return {x,y,z};
}


template <typename Function>
static FloatMat33 finiteDeriv(Function f,FloatMat33 &m)
{
  FloatMat33 result = zeroMat33();

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      result.values[i][j] = finiteDeriv(f,m.values[i][j]);
    }
  }

  return result;
}



template <typename T>
static void
  assertNearHelper(
    const T& actual,
    const T& expected,
    float tolerance,
    const char *file,
    int line
  )
{
  float delta = differenceBetween(actual,expected);

  if (delta<=tolerance) {
    return;
  }

  cerr << "file: " << file << "\n";
  cerr << "line: " << line << "\n";
  cerr << "actual: " << actual << "\n";
  cerr << "expected: " << expected << "\n";
  cerr << "delta: " << delta << "\n";
  cerr << "tolerance: " << tolerance << "\n";
  assert(false);
}


namespace tests {


static void testRowAndCol()
{
  float values[3][3] = {
    {1,2,3},
    {4,5,6},
    {7,8,9}
  };

  FloatMat33 a(values);

  assert(row(a,0)[0]==1);
  assert(row(a,0)[1]==2);
  assert(row(a,0)[2]==3);

  assert(row(a,2)[0]==7);
  assert(row(a,2)[1]==8);
  assert(row(a,2)[2]==9);

  assert(col(a,0)[0]==1);
  assert(col(a,0)[1]==4);
  assert(col(a,0)[2]==7);

  assert(col(a,2)[0]==3);
  assert(col(a,2)[1]==6);
  assert(col(a,2)[2]==9);
}


static void testMat33Inv()
{
  {
    FloatMat33 a = mat33Identity();
    FloatMat33 a_inv = mat33Inv(a);
    assertNear(a_inv,mat33Identity(),1e-4);
  }
  {
    RandomEngine random_engine(/*seed*/1);
    FloatMat33 a = randomMat33(random_engine);
    FloatMat33 a_inv = mat33Inv(a);
    FloatMat33 a_times_a_inv = a*a_inv;
    assertNear(a_times_a_inv,mat33Identity(),1e-4);
  }
}


static void testQRDecomposition()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);
  QRDecomposition<float> qr = qrDecomposed(a);
  const FloatMat33 &q = qr.q;
  const FloatMat33 &r = qr.r;
  assertNear(q*r,a,0);
}


static void testScalarConstantEvaluator()
{
  float dresult = 1;
  float result = evalAndAddDeriv(5,dresult);
  assert(result==5);
}


static void testScalarAddEvaluator()
{
  {
    float a = 1;
    float da = 0;
    float b = 2;
    float db = 0;
    float dresult = 1;
    float result = evalAndAddDeriv(dual(a,da) + dual(b,db),dresult);

    assert(da==1);
    assert(db==1);
    assert(result==3);
  }
  {
    float a = 1;
    float da = 0;
    float dresult = 1;
    float result = evalAndAddDeriv(dual(a,da) + 2,dresult);

    assert(da==1);
    assert(result==3);
  }
}


static void testScalarSubEvaluator()
{
  float a = 1;
  float da = 0;
  float b = 2;
  float db = 0;
  float dresult = 1;
  float result = evalAndAddDeriv(dual(a,da) - dual(b,db),dresult);

  assert(da== 1);
  assert(db==-1);
  assert(result==-1);
}


static void testScalarMulEvaluator()
{
  float a = 1;
  float da = 0;
  float b = 2;
  float db = 0;
  float dresult = 1;
  float result = evalAndAddDeriv(dual(a,da) * dual(b,db),dresult);

  assert(da==2);
  assert(db==1);
  assert(result==2);
}


static void testScalarDivEvaluator()
{
  float a = 1;
  float da = 0;
  float b = 2;
  float db = 0;
  float dresult = 1;
  float result = evalAndAddDeriv(dual(a,da) / dual(b,db),dresult);
  auto f = [&]{ return a/b; };

  assertNear(da,finiteDeriv(f,a),1e-4);
  assertNear(db,finiteDeriv(f,b),1e-4);
  assert(result==a/b);
}


static void testScalarMagEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 v = randomVec3(random_engine);
  FloatVec3 dv{0,0,0};
  float dresult = 1;
  float result = evalAndAddDeriv(mag(dual(v,dv)),dresult);
  auto f = [&]{ return mag(v); };
  assert(result==f());
  assertNear(dv,finiteDeriv(f,v),1e-4);
}


static void testSqrtEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  float x = randomFloat(0,1,random_engine);
  float dx = 0;
  float dresult = 1;
  float result = evalAndAddDeriv(sqrt(dual(x,dx)),dresult);
  assert(result==sqrtf(x));
  auto f = [&]{ return sqrtf(x); };
  assertNear(dx,finiteDeriv(f,x),1e-5);
}


static void testCosEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  float a = randomFloat(-M_PI,M_PI,random_engine);
  float da = 0;
  float dresult = randomFloat(-1,1,random_engine);
  float result = evalAndAddDeriv(cos(dual(a,da)),dresult);
  assert(result==cosf(a));
  auto f = [&]{ return cosf(a); };
  assertNear(da,finiteDeriv(f,a),.005);
}


static void testSinEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  float a = randomFloat(-M_PI,M_PI,random_engine);
  float da = 0;
  float dresult = randomFloat(-1,1,random_engine);
  float result = evalAndAddDeriv(sin(dual(a,da)),dresult);
  assert(result==sinf(a));
  auto f = [&]{ return sinf(a); };
  assertNear(da,finiteDeriv(f,a),.005);
}


static void testDualVec3Evaluator()
{
  FloatVec3 a{1,2,3};
  FloatVec3 da{0,0,0};
  FloatVec3 dresult{1,0,0};
  FloatVec3 result = evalAndAddDeriv(dual(a,da),dresult);
  assert(result==a);
  auto f = [&]{ return a.x(); };
  assertNear(da,finiteDeriv(f,a),1e-4);
}


static void testCofactorEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      FloatMat33 dmat = zeroMat33();
      float dresult = 1;
      float result = evalAndAddDeriv(cofactor(dual(mat,dmat),i,j),dresult);

      auto f = [&]{ return cofactor(mat,i,j); };
      assert(result==f());
      assertNear(dmat,finiteDeriv(f,mat),1e-4);
    }
  }
}


static void testMat33Evaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = zeroMat33();
  auto m = dual(mat,dmat);
  using Element = ScalarAdd<DualFloat,float>;
  Element values[3][3] = {
    {(m[0][0] + 1).expr,
     (m[0][1] + 1).expr,
     (m[0][2] + 1).expr},
    {(m[1][0] + 1).expr,
     (m[1][1] + 1).expr,
     (m[1][2] + 1).expr},
    {(m[2][0] + 1).expr,
     (m[2][1] + 1).expr,
     (m[2][2] + 1).expr},
  };

  auto e = Mat33<Element>(values);
  FloatMat33 dresult = randomMat33(random_engine);

  FloatMat33 result = evalAndAddDeriv(e,dresult);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      assert(result[i][j]==mat[i][j]+1);
    }
  }

  assertNear(dmat,dresult,0);
}


static void testColEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = zeroMat33();
  float dresult = 1;
  float result = evalAndAddDeriv(col(dual(mat,dmat),0)[0],dresult);
  assert(result==mat[0][0]);
  auto f = [&]{ return mat[0][0]; };
  assertNear(dmat,finiteDeriv(f,mat),1e-5);

}


static void testVec3Evaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = zeroMat33();
  FloatVec3 dresult = randomVec3(random_engine);
  FloatVec3 result = evalAndAddDeriv(vec3(col(dual(mat,dmat),0)),dresult);
  ColRef<FloatMat33> c{mat,0};
  assert(result==vec3(c));
  auto f = [&]{ return weightedSum(vec3(col(mat,0)),dresult); };
  assertNear(dmat,finiteDeriv(f,mat),1e-4);
}


static void testCofactorMatrixEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = zeroMat33();
  FloatMat33 dresult = randomMat33(random_engine);
  FloatMat33 result = evalAndAddDeriv(cofactorMatrix(dual(mat,dmat)),dresult);
  assert(result==cofactorMatrix(mat));
  auto f = [&]{ return weightedSum(cofactorMatrix(mat),dresult); };
  assertNear(dmat,finiteDeriv(f,mat),2e-4);
}


static void testDeterminantEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = zeroMat33();
  float dresult = 1;
  float result = evalAndAddDeriv(determinant(dual(mat,dmat)),dresult);
  auto f = [&]{ return determinant(mat); };
  assert(result==f());
  assertNear(dmat,finiteDeriv(f,mat),1e-4);
}


static void testTransposeEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = zeroMat33();
  FloatMat33 dresult = randomMat33(random_engine);
  FloatMat33 result = evalAndAddDeriv(transpose(dual(mat,dmat)),dresult);
  assertNear(result,transpose(mat),0);
  auto f = [&]{ return weightedSum(transpose(mat),dresult); };
  assertNear(dmat,finiteDeriv(f,mat),2e-4);
}


static void testMat33MulEvaluator()
{
  {
    RandomEngine random_engine(/*seed*/1);
    FloatMat33 a = randomMat33(random_engine);
    FloatMat33 b = randomMat33(random_engine);
    FloatMat33 da = zeroMat33();
    FloatMat33 db = zeroMat33();
    FloatMat33 dresult = randomMat33(random_engine);
    FloatMat33 result = evalAndAddDeriv(dual(a,da)*dual(b,db),dresult);
    assert(result==a*b);
    auto f = [&]{ return weightedSum(a*b,dresult); };
    assertNear(da,finiteDeriv(f,a),1e-4);
    assertNear(db,finiteDeriv(f,b),1e-4);
  }
  {
    RandomEngine random_engine(/*seed*/1);
    FloatMat33 a = randomMat33(random_engine);
    FloatMat33 b = randomMat33(random_engine);
    FloatMat33 c = randomMat33(random_engine);
    FloatMat33 da = zeroMat33();
    FloatMat33 db = zeroMat33();
    FloatMat33 dc = zeroMat33();
    FloatMat33 dresult = randomMat33(random_engine);
    FloatMat33 result =
      evalAndAddDeriv(dual(a,da)*dual(b,db)*dual(c,dc),dresult);
    assert(result==a*b*c);
    auto f = [&]{ return weightedSum(a*b*c,dresult); };
    assertNear(da,finiteDeriv(f,a),1e-4);
    assertNear(db,finiteDeriv(f,b),1e-4);
    assertNear(dc,finiteDeriv(f,c),1e-4);
  }
}


static void testMat33DivEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  float divisor = randomFloat(-1,1,random_engine);
  float ddivisor = 0;

  FloatMat33 dmat = zeroMat33();
  FloatMat33 dresult = randomMat33(random_engine);
  auto e = dual(mat,dmat)/dual(divisor,ddivisor);
  FloatMat33 result = evalAndAddDeriv(e,dresult);
  auto f = [&]{ return weightedSum(mat/divisor,dresult); };
  assertNear(result,mat/divisor,0);
  assertNear(dmat,finiteDeriv(f,mat),1e-4);
  assertNear(ddivisor,finiteDeriv(f,divisor),1e-4);
}


static void testMat33InvEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);

  FloatMat33 dmat = zeroMat33();
  FloatMat33 dresult = randomMat33(random_engine);
  FloatMat33 result =
    evalAndAddDeriv(mat33Inv(dual(mat,dmat)),dresult);
  auto f = [&]{ return weightedSum(mat33Inv(mat),dresult); };
  assertNear(result,mat33Inv(mat),0);
  assertNear(dmat,finiteDeriv(f,mat),1e-3);
}


static void testDotEvaluator()
{
  {
    RandomEngine random_engine(/*seed*/1);
    FloatVec3 a = randomVec3(random_engine);
    FloatVec3 b = randomVec3(random_engine);
    FloatVec3 da{0,0,0};
    FloatVec3 db{0,0,0};

    auto e = dot(dual(a,da),dual(b,db));
    float dresult = 1;
    float result = evalAndAddDeriv(e,dresult);
    auto f = [&]{ return dot(a,b); };
    assert(result==f());
    assertNear(da,finiteDeriv(f,a),1e-4);
    assertNear(db,finiteDeriv(f,b),1e-4);
  }
  {
    float dv = 0;
    {
      auto a = vec3(2,dual(3,dv),4);
      auto b = vec3(5,6,7);
      float result = evalAndAddDeriv(dot(a,b),1);
      assert(result==dot(vec3(2,3,4),vec3(5,6,7)));
    }
    assertNear(dv,6.0f,1e-4);
  }
}


static void testVec3AddEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 a = randomVec3(random_engine);
  FloatVec3 b = randomVec3(random_engine);
  FloatVec3 da{0,0,0};
  FloatVec3 db{0,0,0};
  FloatVec3 dresult = randomVec3(random_engine);
  FloatVec3 result = evalAndAddDeriv(dual(a,da) + dual(b,db),dresult);
  auto f = [&]{ return weightedSum(a+b,dresult); };
  assertNear(result,a+b,0);
  assertNear(finiteDeriv(f,a),da,1e-4);
  assertNear(finiteDeriv(f,b),db,1e-4);
}


static void testVec3SubEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 a = randomVec3(random_engine);
  FloatVec3 b = randomVec3(random_engine);
  FloatVec3 da{0,0,0};
  FloatVec3 db{0,0,0};
  FloatVec3 dresult = randomVec3(random_engine);
  FloatVec3 result = evalAndAddDeriv(dual(a,da) - dual(b,db),dresult);
  auto f = [&]{ return weightedSum(a-b,dresult); };
  assertNear(result,a-b,0);
  assertNear(finiteDeriv(f,a),da,1e-4);
  assertNear(finiteDeriv(f,b),db,1e-4);
}


static void testVec3ScalarMulEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 a = randomVec3(random_engine);
  float b = randomFloat(-1,1,random_engine);
  FloatVec3 da{0,0,0};
  float db = 0;
  FloatVec3 dresult = randomVec3(random_engine);
  FloatVec3 result = evalAndAddDeriv(dual(a,da) * dual(b,db),dresult);
  auto f = [&]{ return weightedSum(a*b,dresult); };
  assertNear(result,a*b,0);
  assertNear(finiteDeriv(f,a),da,1e-4);
  assertNear(finiteDeriv(f,b),db,1e-4);
}


static void testScalarVec3MulEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 a = randomVec3(random_engine);
  float b = randomFloat(-1,1,random_engine);
  FloatVec3 da{0,0,0};
  float db = 0;
  FloatVec3 dresult = randomVec3(random_engine);
  FloatVec3 result = evalAndAddDeriv(dual(b,db) * dual(a,da),dresult);
  auto f = [&]{ return weightedSum(b*a,dresult); };
  assertNear(result,b*a,0);
  assertNear(finiteDeriv(f,a),da,1e-4);
  assertNear(finiteDeriv(f,b),db,1e-4);
}


static void testVec3DivEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 a = randomVec3(random_engine);
  float b = randomFloat(-1,1,random_engine);
  FloatVec3 da{0,0,0};
  float db = 0;
  FloatVec3 dresult = randomVec3(random_engine);
  FloatVec3 result = evalAndAddDeriv(dual(a,da) / dual(b,db),dresult);
  auto f = [&]{ return weightedSum(a/b,dresult); };
  assertNear(result,a/b,0);
  assertNear(finiteDeriv(f,a),da,1e-4);
  assertNear(finiteDeriv(f,b),db,1e-4);
}


static void testRotXEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  float a = randomFloat(-1,1,random_engine);
  float da = 0;
  FloatMat33 dresult = randomMat33(random_engine);
  FloatMat33 result = evalAndAddDeriv(rotX(dual(a,da)),dresult);
  assert(result==rotX(a));
  auto f = [&]{ return weightedSum(rotX(a),dresult); };
  assertNear(da,finiteDeriv(f,a),1e-4);
}


static void testColumnsEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 c1 = randomVec3(random_engine);
  FloatVec3 dc1{0,0,0};
  FloatVec3 c2 = randomVec3(random_engine);
  FloatVec3 dc2{0,0,0};
  FloatVec3 c3 = randomVec3(random_engine);
  FloatVec3 dc3{0,0,0};
  FloatMat33 dresult = randomMat33(random_engine);
  FloatMat33 result =
    evalAndAddDeriv(columns(dual(c1,dc1),dual(c2,dc2),dual(c3,dc3)),dresult);
  assert(result==columns(c1,c2,c3));
  auto f = [&]{ return weightedSum(columns(c1,c2,c3),dresult); };
  assertNear(dc1,finiteDeriv(f,c1),2e-4);
  assertNear(dc2,finiteDeriv(f,c2),2e-4);
  assertNear(dc3,finiteDeriv(f,c3),2e-4);
}


static void testExprVar1()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);
  FloatMat33 da = zeroMat33();
  Mat33ExprVar<DualMat33> a_var(dual(a,da));
  FloatMat33 dresult = randomMat33(random_engine);
  FloatMat33 result = evalAndAddDeriv(a_var,dresult);
  assertNear(result,a,0);
  auto f = [&]{ return weightedSum(a,dresult); };
  assertNear(da,finiteDeriv(f,a),1e-4);
}


static void testExprVar2()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 a = randomVec3(random_engine);
  FloatVec3 da{0,0,0};
  FloatVec3 dresult = randomVec3(random_engine);
  {
    Vec3ExprVar<DualVec3> a_var(dual(a,da));
    Evaluator<DualVec3> eval(a_var.dual());
    eval.addDeriv(dresult);
  }
  auto f = [&]{ return weightedSum(a,dresult); };
  assertNear(da,finiteDeriv(f,a),1e-4);
}


static void testExprVar3()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);
  FloatMat33 da = zeroMat33();
  FloatVec3 dresult = randomVec3(random_engine);

  {
    Mat33ExprVar<DualMat33> a_var(dual(a,da));
    Vec3ExprVar<decltype(vec3(col(a_var,0)))> a1(vec3(col(a_var,0)));
    FloatVec3 result = evalAndAddDeriv(a1,dresult);
    assertNear(result,vec3(col(a,0)),0);
  }
  auto f = [&]{ return weightedSum(vec3(col(a,0)),dresult); };
  assertNear(da,finiteDeriv(f,a),1e-4);
}



static void testQRDecompositionEvaluator()
{
  RandomEngine random_engine(/*seed*/1);

  FloatMat33 a = randomMat33(random_engine);
  FloatMat33 da = zeroMat33();
  FloatMat33 dq = randomMat33(random_engine);
  FloatMat33 dr = randomMat33(random_engine);
  QRDecomposition<float> dresult{dq,dr};
  QRDecomposition<float> result =
    evalAndAddDeriv(qrDecomposed(dual(a,da)),dresult);

  auto f = [&]{
    auto qr = qrDecomposed(a);
    auto &q = qr.q;
    auto &r = qr.r;
    return weightedSum(q,dq) + weightedSum(r,dr);
  };

  assertNear(result,qrDecomposed(a),1e-4);
  assertNear(da,finiteDeriv(f,a),2e-4);
}

}


int main()
{
  tests::testRowAndCol();
  tests::testMat33Inv();
  tests::testQRDecomposition();

  tests::testScalarConstantEvaluator();
  tests::testScalarAddEvaluator();
  tests::testScalarSubEvaluator();
  tests::testScalarMulEvaluator();
  tests::testScalarDivEvaluator();
  tests::testSqrtEvaluator();
  tests::testCosEvaluator();
  tests::testSinEvaluator();

  tests::testDualVec3Evaluator();
  tests::testDotEvaluator();
  tests::testVec3AddEvaluator();
  tests::testVec3SubEvaluator();
  tests::testVec3ScalarMulEvaluator();
  tests::testScalarVec3MulEvaluator();
  tests::testVec3DivEvaluator();
  tests::testScalarMagEvaluator();

  tests::testMat33Evaluator();
  tests::testColEvaluator();
  tests::testVec3Evaluator();
  tests::testCofactorEvaluator();
  tests::testCofactorMatrixEvaluator();
  tests::testDeterminantEvaluator();
  tests::testTransposeEvaluator();
  tests::testMat33MulEvaluator();
  tests::testMat33DivEvaluator();
  tests::testMat33InvEvaluator();
  tests::testRotXEvaluator();
  tests::testColumnsEvaluator();

  tests::testExprVar1();
  tests::testExprVar2();
  tests::testExprVar3();
  tests::testQRDecompositionEvaluator();
}

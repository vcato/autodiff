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


using std::cerr;
using std::ostream;
using std::max;


#define assertNear(actual,expected,tolerance) \
  (assertNearHelper(actual,expected,tolerance,__FILE__,__LINE__))



static float differenceBetween(float a,float b)
{
  return fabsf(a-b);
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


static float max(float a,float b,float c)
{
  return std::max(std::max(a,b),c);
}


static float differenceBetween(const FloatVec3 &a,const FloatVec3 &b)
{
  float dx = differenceBetween(a.x(), b.x());
  float dy = differenceBetween(a.y(), b.y());
  float dz = differenceBetween(a.z(), b.z());

  return max(dx,dy,dz);
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


namespace tests {

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


static void testScalarConstantEvaluator()
{
  float dresult = 1;
  float result = evalAndAddDeriv(expr(5),dresult);
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
    float result = evalAndAddDeriv(expr(a,da) + expr(b,db),dresult);

    assert(da==1);
    assert(db==1);
    assert(result==3);
  }
  {
    float a = 1;
    float da = 0;
    float dresult = 1;
    float result = evalAndAddDeriv(expr(a,da) + expr(2),dresult);

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
  float result = evalAndAddDeriv(expr(a,da) - expr(b,db),dresult);

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
  float result = evalAndAddDeriv(expr(a,da) * expr(b,db),dresult);

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
  float result = evalAndAddDeriv(expr(a,da) / expr(b,db),dresult);
  auto f = [&]{ return a/b; };

  assertNear(da,finiteDeriv(f,a),1e-4);
  assertNear(db,finiteDeriv(f,b),1e-4);
  assert(result==a/b);
}


static void testCosEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  float a = randomFloat(-M_PI,M_PI,random_engine);
  float da = 0;
  float dresult = randomFloat(-1,1,random_engine);
  float result = evalAndAddDeriv(cos(expr(a,da)),dresult);
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
  float result = evalAndAddDeriv(sin(expr(a,da)),dresult);
  assert(result==sinf(a));
  auto f = [&]{ return sinf(a); };
  assertNear(da,finiteDeriv(f,a),.005);
}


static void testDualVec3Evaluator()
{
  FloatVec3 a{1,2,3};
  FloatVec3 da{0,0,0};
  FloatVec3 dresult{1,0,0};
  evalAndAddDeriv(expr(a,da),dresult);
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
      float result = evalAndAddDeriv(cofactor(expr(mat,dmat),i,j),dresult);

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
  auto m = expr(mat,dmat);
  using Element = ScalarAdd<DualFloat,float>;
  Element values[3][3] = {
    {(m[0][0] + expr(1)).expr,
     (m[0][1] + expr(1)).expr,
     (m[0][2] + expr(1)).expr},
    {(m[1][0] + expr(1)).expr,
     (m[1][1] + expr(1)).expr,
     (m[1][2] + expr(1)).expr},
    {(m[2][0] + expr(1)).expr,
     (m[2][1] + expr(1)).expr,
     (m[2][2] + expr(1)).expr},
  };

  auto e = expr(Mat33<Element>(values));
  FloatMat33 dresult = randomMat33(random_engine);

  evalAndAddDeriv(e,dresult);

  assertNear(dmat,dresult,0);
}


static void testCofactorMatrixEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = zeroMat33();
  FloatMat33 dresult = randomMat33(random_engine);
  evalAndAddDeriv(cofactorMatrix(expr(mat,dmat)),dresult);
  auto f = [&]{ return weightedSum(cofactorMatrix(mat),dresult); };
  assertNear(dmat,finiteDeriv(f,mat),2e-4);
}


static void testDeterminantEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = zeroMat33();
  float dresult = 1;
  float result = evalAndAddDeriv(determinant(expr(mat,dmat)),dresult);
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
  FloatMat33 result = evalAndAddDeriv(transpose(expr(mat,dmat)),dresult);
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
    evalAndAddDeriv(expr(a,da)*expr(b,db),dresult);
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
    evalAndAddDeriv(expr(a,da)*expr(b,db)*expr(c,dc),dresult);
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
  auto e = expr(mat,dmat)/expr(divisor,ddivisor);
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
    evalAndAddDeriv(mat33Inv(expr(mat,dmat)),dresult);
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

    auto e = dot(expr(a,da),expr(b,db));
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
      auto a = vec3(expr(2),expr(3,dv),expr(4));
      auto b = vec3(expr(5),expr(6),   expr(7));
      evalAndAddDeriv(dot(a,b),1);
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
  FloatVec3 result = evalAndAddDeriv(expr(a,da) + expr(b,db),dresult);
  auto f = [&]{ return weightedSum(a+b,dresult); };
  assertNear(result,a+b,0);
  assertNear(finiteDeriv(f,a),da,1e-4);
  assertNear(finiteDeriv(f,b),db,1e-4);
}


static void testRotXEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  float a = randomFloat(-1,1,random_engine);
  float da = 0;
  FloatMat33 dresult = randomMat33(random_engine);
  evalAndAddDeriv(rotX(expr(a,da)),dresult);
  auto f = [&]{ return weightedSum(rotX(a),dresult); };
  assertNear(da,finiteDeriv(f,a),1e-4);
}


}


int main()
{
  tests::testMat33Inv();

  tests::testScalarConstantEvaluator();
  tests::testScalarAddEvaluator();
  tests::testScalarSubEvaluator();
  tests::testScalarMulEvaluator();
  tests::testScalarDivEvaluator();
  tests::testCosEvaluator();
  tests::testSinEvaluator();

  tests::testDualVec3Evaluator();
  tests::testDotEvaluator();
  tests::testVec3AddEvaluator();

  tests::testCofactorEvaluator();
  tests::testMat33Evaluator();
  tests::testCofactorMatrixEvaluator();
  tests::testDeterminantEvaluator();
  tests::testTransposeEvaluator();
  tests::testMat33MulEvaluator();
  tests::testMat33DivEvaluator();
  tests::testMat33InvEvaluator();
  tests::testRotXEvaluator();
}

#include <iostream>
#include "mat33.hpp"
#include "vec3.hpp"
#include "random.hpp"
#include "dualfloat.hpp"
#include "evalandaddderiv.hpp"
#include "scalarexpr.hpp"
#include "vec3expr.hpp"
#include "mat33expr.hpp"
#include "finitederivfloat.hpp"
#include "assertnear.hpp"

using autodiff::dual;
using autodiff::vec3;
using autodiff::Vec3ExprType;
using autodiff::ScalarExprType;
using std::sin;
using std::cos;
using autodiff::ScalarExpr;
using autodiff::Mat33Expr;



namespace {
template <typename Angle>
struct Rot {
  Angle angle;
};
}


template <typename Angle>
struct autodiff::Evaluator<Rot<Angle>> {
  ScalarExprVar<Angle> angle;
  SCALAR_VAR(          cos_angle,  cos(angle));
  SCALAR_VAR(          sin_angle,  sin(angle));

  Evaluator(Rot<Angle> arg)
  : angle(arg.angle)
  {
  }

  FloatMat33 value() const
  {
    float c = cos_angle.value();
    float s = sin_angle.value();

    float values[3][3] = {
      { c,s,0},
      {-s,c,0},
      { 0,0,1}
    };

    return {values};
  }

  void addDeriv(FloatMat33 dresult)
  {
    cos_angle.addDeriv( dresult[0][0]);
    cos_angle.addDeriv( dresult[1][1]);
    sin_angle.addDeriv( dresult[0][1]);
    sin_angle.addDeriv(-dresult[1][0]);
  }
};


namespace {
template <typename AngleExpr,typename Angle=ScalarExprType<AngleExpr>>
static Mat33Expr<Rot<Angle>> rot(const AngleExpr &angle)
{
  return {{internal(angle)}};
}
}


namespace {
template <typename Offset>
struct TransX {
  Offset offset;
};
}


template <typename Offset>
struct autodiff::Evaluator<TransX<Offset>> {
  ScalarExprVar<Offset> offset;

  Evaluator(TransX<Offset> expr)
  : offset(expr.offset)
  {
  }

  FloatMat33 value() const
  {
    float v = offset.value();

    float values[3][3] = {
      {1,0,v},
      {0,1,0},
      {0,0,1}
    };

    return FloatMat33(values);
  }

  void addDeriv(const FloatMat33 &deriv)
  {
    offset.addDeriv(deriv[0][2]);
  }
};


namespace {
template <typename OffsetExpr,typename Offset=ScalarExprType<OffsetExpr>>
static Mat33Expr<TransX<Offset>> transX(const OffsetExpr &offset)
{
  return {{autodiff::internal(offset)}};
}
}


namespace {
template <typename V>
struct Mag2 {
  V v;
};
}


namespace autodiff {
template <typename V>
struct autodiff::Evaluator<Mag2<V>> {
  Vec3ExprVar<V> v;
  SCALAR_VAR(result,dot(v,v));

  Evaluator(Mag2<V> arg)
  : v(arg.v)
  {
  }

  float value() const
  {
    return result.value();
  }

  void addDeriv(float dresult)
  {
    result.addDeriv(dresult);
  }
};
}


namespace {
template <typename VExpr,typename V = Vec3ExprType<VExpr>>
static ScalarExpr<Mag2<V>> mag2(const VExpr &v)
{
  return {{internal(v)}};
}
}


namespace {
template <typename Rot1,typename Rot2>
struct Error {
  Rot1 rot1;
  Rot2 rot2;
};
}


namespace autodiff {
template <typename Rot1,typename Rot2>
struct autodiff::Evaluator<Error<Rot1,Rot2>> {
  ScalarExprVar<Rot1> rot1;
  ScalarExprVar<Rot2> rot2;
  MAT33_VAR(          body1_transform, rot(rot1)                            );
  MAT33_VAR(          body2_transform, body1_transform*transX(10)*rot(rot2) );
  VEC3_VAR(           local          , vec3(20,0,1)                         );
  VEC3_VAR(           global         , vec3(25,0,1)                         );
  VEC3_VAR(           predicted      , body2_transform*local                );
  SCALAR_VAR(         result         , mag2(predicted-global                ));

  Evaluator(Error<Rot1,Rot2> arg)
  : rot1(arg.rot1),
    rot2(arg.rot2)
  {
  }

  float value() const
  {
    return result.value();
  }

  void addDeriv(float dresult)
  {
    result.addDeriv(dresult);
  }
};
}


namespace {
template <typename Rot1,typename Rot2>
static ScalarExpr<Error<Rot1,Rot2>> error(Rot1 rot1,Rot2 rot2)
{
  return {{autodiff::internal(rot1),autodiff::internal(rot2)}};
}
}


static float error(float rot1,float rot2)
{
  return autodiff::Evaluator<Error<float,float>>({rot1,rot2}).value();
}


static void test1()
{
  float rot1 = 0;
  float rot2 = 0;

  float e = error(rot1,rot2);
  assertNear(e,25.0f,0);
}


static void test2()
{
  float rot1 = 1.2;
  float rot2 = 7.5;
  float drot1 = 0;
  float drot2 = 0;
  float dresult = 1;

  float result =
    evalAndAddDeriv(error(dual(rot1,drot1),dual(rot2,drot2)),dresult);

  assert(result == error(rot1,rot2));

  auto f = [&]{ return error(rot1,rot2); };

  float tolerance = 0.1;
  assertNear(finiteDeriv(f,rot1),drot1,tolerance);
  assertNear(finiteDeriv(f,rot2),drot2,tolerance);
}


int main()
{
  test1();
  test2();
}

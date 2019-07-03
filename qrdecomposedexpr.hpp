#include "mat33.hpp"
#include "mat33expr.hpp"


namespace autodiff {


template <typename A>
struct QRDecomposed {
  A a;
};


template <typename T>
T differenceBetween(const QRDecomposition<T> &a,const QRDecomposition<T> &b)
{
  return std::max(differenceBetween(a.q,b.q),differenceBetween(a.r,b.r));
}



template <typename AExpr,typename A=Mat33ExprType<AExpr>>
Mat33Expr<QRDecomposed<A>> qrDecomposed(const AExpr &a)
{
  return {{internal(a)}};
}

template <typename A>
struct Evaluator<QRDecomposed<A>> {
  Mat33ExprVar<A> a;

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

  // a1 = q1*r11
  // a2 = q1*r12 + q2*r22
  // a3 = q1*r13 + q2*r23 + q3*r33

  // each ui is orthogonal to uj when i!=j
  //   ui = ai - sum{j if j!=i}(dot(qj,ai)*qj)
  // qi is the normalized version of ui
  //   qi = ui/mag(ui)

  // a1 = <a11,a21,a31>
  // a2 = <a12,a22,a32>
  // a3 = <a13,a23,a33>
  VEC3_VAR(a1, vec3(col(a,0)));
  VEC3_VAR(a2, vec3(col(a,1)));
  VEC3_VAR(a3, vec3(col(a,2)));

  // We should be able to use Vec3ExprVar<decltype(a1)> here.
  decltype(a1) &u1 = a1;

  // u1 = a1
  // q1 = u1/mag(u1)
  // q1*mag(u1) = u1
  // u1 = q1*mag(u1)
  // a1 = q1*r11
  // r11 = mag(u1)
  SCALAR_VAR(r11, mag(u1));

  // q1 = u1/mag(u1)
  // q1 = u1/r11
  VEC3_VAR(q1, u1/r11);

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
  SCALAR_VAR(r12, dot(a2,q1));

  // u2 = a2 - q1*dot(q1,a2)
  // u2 = a2 - q1*r12
  VEC3_VAR(u2, a2 - q1*r12);

  // r22 = mag(u2)
  SCALAR_VAR(r22, mag(u2));

  // q2 = u2/mag(u2)
  // q2 = u2/r22
  VEC3_VAR(q2, u2/r22);

  // u3 = a3 - q1*dot(q1,a3) - q2*dot(q2,a3)
  // u3 + q1*dot(q1,a3) + q2*dot(q2,a3) = a3
  // a3 = u3 + q1*dot(q1,a3) + q2*dot(q2,a3)
  // a3 = q1*dot(q1,a3) + q2*dot(q2,a3) + u3
  // a3 = q1*dot(q1,a3) + q2*dot(q2,a3) + q3*mag(u3)
  // a3 = q1*r13 + q2*r23 + q3*r33
  // r13 = dot(q1,a3)
  SCALAR_VAR(r13 ,dot(q1,a3));

  // r23 = dot(q2,a3)
  SCALAR_VAR(r23 ,dot(q2,a3));

  // r33 = mag(u3)
  // u3 = a3 - q1*dot(q1,a3) - q2*dot(q2,a3)
  // u3 = a3 - q1*r13 - q2*r23
  VEC3_VAR(u3, a3 - q1*r13 - q2*r23);

  // r33 = mag(u3)
  SCALAR_VAR(r33, mag(u3));

  // q3 = u3/mag(u3)
  // q3 = u3/r33
  VEC3_VAR(q3, u3/r33);

  MAT33_VAR(q, columns(q1,q2,q3));

  Evaluator(const QRDecomposed<A> &expr)
  : a({expr.a})
  {
  }

  QRDecomposition<float> value() const
  {
    float r_values[3][3] = {
      {r11.value(),r12.value(),r13.value()},
      {          0,r22.value(),r23.value()},
      {          0,          0,r33.value()}
    };

    return {q.value(),Mat33<float>(r_values)};
  }

  void addDeriv(const QRDecomposition<float> &deriv)
  {
    q.addDeriv(deriv.q);
    r11.addDeriv(deriv.r[0][0]);
    r12.addDeriv(deriv.r[0][1]);
    r13.addDeriv(deriv.r[0][2]);
    r22.addDeriv(deriv.r[1][1]);
    r23.addDeriv(deriv.r[1][2]);
    r33.addDeriv(deriv.r[2][2]);
  }
};

}

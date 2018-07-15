#include "dualmat33.hpp"
#include "scalarexpr.hpp"
#include "evaluator.hpp"


template <typename Expr>
struct Mat33Expr {
  Expr expr;
};


template <typename Expr>
struct ExprVar<Mat33Expr<Expr>> {
  Evaluator<Expr> eval;
  FloatMat33 value_member = eval.value();
  FloatMat33 deriv_member = zeroMat33();

  ExprVar(const Mat33Expr<Expr> &expr)
  : eval(expr.expr)
  {
  }

  FloatMat33 value() const { return value_member; }

  void addDeriv(const FloatMat33 &dvalue)
  {
    deriv_member += dvalue;
  }

  ~ExprVar()
  {
    eval.addDeriv(deriv_member);
  }
};


template <>
struct Mat33Expr<DualMat33> {
  DualMat33 expr;

  template <typename Self>
  static auto element(Self &arg,int i,int j)
  {
    assert(i>=0);
    assert(i<3);
    assert(j>=0);
    assert(j<3);
    return ScalarExpr<DualFloat>{arg.expr[i][j]};
  }

  auto operator[](int i)             { return row(*this,i); }
  auto operator[](int i) const       { return row(*this,i); }
};


inline Mat33Expr<DualMat33> expr(const FloatMat33 &value,FloatMat33 &deriv)
{
  return {dual(value,deriv)};
}


template <typename T>
struct Evaluator<Mat33<T>> {
  Mat33<T> m_expr;
  Mat33<Evaluator<T>> m_eval;

  Evaluator(const Mat33<T> &arg)
  : m_expr(arg),
    m_eval(
      {
        {arg[0][0],arg[0][1],arg[0][2]},
        {arg[1][0],arg[1][1],arg[1][2]},
        {arg[2][0],arg[2][1],arg[2][2]}
      }
    )
  {
  }

  FloatMat33 value() const
  {
    float values[3][3];

    for (int i=0; i!=3; ++i) {
      for (int j=0; j!=3; ++j) {
        values[i][j] = m_eval[i][j].value();
      }
    }

    return {values};
  }

  void addDeriv(const FloatMat33 &deriv)
  {
    for (int i=0; i!=3; ++i) {
      for (int j=0; j!=3; ++j) {
        m_eval[i][j].addDeriv(deriv[i][j]);
      }
    }
  }
};


template <typename A,typename B>
struct Mat33Mul {
  A a;
  B b;
};


template <typename A,typename B>
Mat33Expr<Mat33Mul<A,B>> operator*(Mat33Expr<A> a,Mat33Expr<B> b)
{
  return {{a.expr,b.expr}};
}


template <typename A,typename B>
struct Evaluator<Mat33Mul<A,B>> {
  Evaluator<A> a_eval;
  FloatMat33 a;
  Evaluator<B> b_eval;
  FloatMat33 b;

  Evaluator(const Mat33Mul<A,B> &expr)
  : a_eval(expr.a),
    a(a_eval.value()),
    b_eval(expr.b),
    b(b_eval.value())
  {
  }

  FloatMat33 value() const { return a*b; }

  void addDeriv(const FloatMat33 &dvalue)
  {
    a_eval.addDeriv(dvalue*transpose(b));
    b_eval.addDeriv(transpose(a)*dvalue);
  }
};




template <typename ExprWrapper,typename Value>
Value evalAndAddDeriv(const ExprWrapper &mat33_expr,const Value& dresult)
{
  auto e = mat33_expr.expr;
  Evaluator<decltype(e)> eval(e);
  Value result = eval.value();
  eval.addDeriv(dresult);
  return result;
}


inline void addDeriv(const DualFloat &dual,float dresult)
{
  dual.deriv += dresult;
}


template <typename T>
inline void addDeriv(const Mat33<T> &m,const FloatMat33 &dm)
{
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      addDeriv(m[i][j],dm[i][j]);
    }
  }
}


template <typename M>
struct Mat33ExprVar {
  Evaluator<M> eval;
  FloatMat33 _value = eval.value();
  FloatMat33 deriv = zeroMat33();

  FloatMat33 value() const { return _value; }
  void addDeriv(const FloatMat33 &dvalue) { deriv += dvalue; }
  DualMat33 dual() { return ::dual(_value,deriv); }

  Mat33ExprVar(const M &m)
  : eval(m)
  {
  }

  ~Mat33ExprVar()
  {
    eval.addDeriv(deriv);
  }
};


namespace {
template <typename T>
static Mat33<float> evaluate(const Mat33<T> &m)
{
  float values[3][3] = {
    {evaluate(m[0][0]),evaluate(m[0][1]),evaluate(m[0][2])},
    {evaluate(m[1][0]),evaluate(m[1][1]),evaluate(m[1][2])},
    {evaluate(m[2][0]),evaluate(m[2][1]),evaluate(m[2][2])},
  };

  return Mat33<float>(values);
}
}


// Specialization in the case of a DualMat33 expression, where we don't
// need to keep the value separately, and we don't need to add the
// derivative in the destructor.
template <>
struct Mat33ExprVar<DualMat33> {
  DualMat33 expr;

  Mat33ExprVar(const DualMat33 &e)
  : expr(e)
  {
  }

  FloatMat33 value() const { return evaluate(expr); }

  void addDeriv(const FloatMat33 &deriv) { ::addDeriv(expr,deriv); }

  DualMat33 dual() { return expr; }
};


namespace {
template <typename A,typename B>
struct Mat33Div {
  A a;
  B b;
};
}


template <typename A,typename B>
inline Mat33Expr<Mat33Div<A,B>> operator/(Mat33Expr<A> a,ScalarExpr<B> b)
{
  return {{a.expr,b.expr}};
}


template <typename A,typename B>
struct Evaluator<Mat33Div<A,B>> {
  Mat33ExprVar<A> a;
  ScalarExprVar<B> b;
  Evaluator<decltype(genMat33Div(expr(a),expr(b)).expr)>
    result{ genMat33Div(expr(a),expr(b)).expr };

  Evaluator(const Mat33Div<A,B> &expr)
  : a(expr.a),
    b(expr.b)
  {
  }

  FloatMat33 value() const
  {
    return result.value();
  }

  void addDeriv(const FloatMat33 &deriv)
  {
    result.addDeriv(deriv);
  }
};


template <typename M>
struct Cofactor {
  M m;
  int i;
  int j;
};


template <typename M>
ScalarExpr<Cofactor<M>> cofactor(const Mat33Expr<M> &arg,int i,int j)
{
  return {{arg.expr,i,j}};
}


template <typename M>
struct Evaluator<Cofactor<M>> {
  Mat33ExprVar<M> m;
  Evaluator<decltype(genCofactor(expr(m),0,0).expr)> result;

  Evaluator(const Cofactor<M> &c)
  : m{c.m},
    result{genCofactor(expr(m),c.i,c.j).expr}
  {
  }

  float value() const { return result.value(); }

  void addDeriv(float dresult) { result.addDeriv(dresult); }
};


template <typename M>
struct CofactorMatrixFunc {
  M m;
};


template <typename M>
Mat33Expr<CofactorMatrixFunc<M>> cofactorMatrix(const Mat33Expr<M> &m)
{
  return {{m.expr}};
}


template <typename M>
struct Evaluator<CofactorMatrixFunc<M>> {
  Mat33ExprVar<M> m;
  Evaluator<decltype(genCofactorMatrix(expr(m)).expr)> result_eval =
    genCofactorMatrix(expr(m)).expr;

  Evaluator(const CofactorMatrixFunc<M> &expr)
  : m(expr.m)
  {
  }

  FloatMat33 value() const
  {
    return result_eval.value();
  }

  void addDeriv(const FloatMat33 &dresult)
  {
    result_eval.addDeriv(dresult);
  }
};


template <typename M>
Mat33Expr<DualMat33> expr(Mat33ExprVar<M> &v)
{
  return {v.dual()};
}


template <typename T>
Mat33Expr<Mat33<T>> expr(const Mat33<T> &expr)
{
  return {expr};
}


template <typename T>
Mat33Expr<Mat33<T>> mat33(const ScalarExpr<T> (&arg)[3][3])
{
  T values[3][3] = {
    {arg[0][0].expr,arg[0][1].expr,arg[0][2].expr},
    {arg[1][0].expr,arg[1][1].expr,arg[1][2].expr},
    {arg[2][0].expr,arg[2][1].expr,arg[2][2].expr},
  };

  return {{values}};
}


template <typename M>
struct Transpose {
  M m;
};


template <typename M>
Mat33Expr<Transpose<M>> transpose(const Mat33Expr<M> &m)
{
  return {{m.expr}};
}


template <typename M>
struct Evaluator<Transpose<M>> {
  Mat33ExprVar<M> m;

  Evaluator(const Transpose<M> &expr)
  : m(expr.m)
  {
  }

  FloatMat33 value() const
  {
    return transpose(m.value());
  }

  void addDeriv(const FloatMat33 &dresult)
  {
    m.addDeriv(transpose(dresult));
  }
};


template <typename M>
struct Determinant {
  M m;
};


template <typename M>
ScalarExpr<Determinant<M>> determinant(const Mat33Expr<M> &m)
{
  const M &expr = m.expr;
  return {{expr}};
}


template <typename M>
struct Evaluator<Determinant<M>> {
  Mat33ExprVar<M> m;
  Evaluator<decltype(genDeterminant(expr(m)).expr)>
    result{genDeterminant(expr(m)).expr};

  Evaluator(const Determinant<M> &expr)
  : m(expr.m)
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


template <typename M>
struct Inv {
  M m;
};


template <typename T>
Mat33Expr<Inv<T>> mat33Inv(const Mat33Expr<T> &arg)
{
  return {{arg.expr}};
}


template <typename M>
struct Evaluator<Inv<M>>
{
  Mat33ExprVar<M> m;
  Evaluator<
    decltype(genMat33Inv(expr(m)).expr)
  > result_eval{genMat33Inv(expr(m)).expr};

  Evaluator(const Inv<M> &expr_arg)
  : m(expr_arg.m)
  {
  }

  FloatMat33 value() const
  {
    return result_eval.value();
  }

  void addDeriv(const FloatMat33 &dresult)
  {
    result_eval.addDeriv(dresult);
  }
};


template <typename Angle>
struct RotX {
  Angle angle;
};


template <typename Angle>
Mat33Expr<RotX<Angle>> rotX(ScalarExpr<Angle> angle)
{
  return {{angle.expr}};
}


template <typename Angle>
struct Evaluator<RotX<Angle>> {
  ScalarExprVar<Angle> angle;
  ScalarExprVar<Cos<Angle>> cosa{cos(expr(angle)).expr};
  ScalarExprVar<Sin<Angle>> sina{sin(expr(angle)).expr};

  Evaluator(const RotX<Angle> &expr)
  : angle(expr.angle)
  {
  }

  FloatMat33 value() const
  {
    return rotX(sina.value(),cosa.value());
  }

  void addDeriv(const FloatMat33 &dvalue)
  {
    cosa.addDeriv( dvalue[1][1]);
    sina.addDeriv( dvalue[1][2]);
    sina.addDeriv(-dvalue[2][1]);
    cosa.addDeriv( dvalue[2][2]);
  }
};


template <typename C1,typename C2,typename C3>
struct ColumnsFunc
{
  C1 c1;
  C2 c2;
  C3 c3;
};


template <typename C1,typename C2,typename C3>
struct Evaluator<ColumnsFunc<C1,C2,C3>> {
  Evaluator<C1> c1_eval;
  Evaluator<C2> c2_eval;
  Evaluator<C3> c3_eval;

  Evaluator(const ColumnsFunc<C1,C2,C3> &expr)
  : c1_eval(expr.c1),
    c2_eval(expr.c2),
    c3_eval(expr.c3)
  {
  }

  FloatMat33 value() const
  {
    FloatVec3 c1 = c1_eval.value();
    FloatVec3 c2 = c2_eval.value();
    FloatVec3 c3 = c3_eval.value();
    return columns(c1,c2,c3);
  }

  void addDeriv(const FloatMat33 &deriv)
  {
    c1_eval.addDeriv(vec3(col(deriv,0)));
    c2_eval.addDeriv(vec3(col(deriv,1)));
    c3_eval.addDeriv(vec3(col(deriv,2)));
  }
};


template <typename C1,typename C2,typename C3>
Mat33Expr<ColumnsFunc<C1,C2,C3>>
  columns(const Vec3Expr<C1> &c1,const Vec3Expr<C2> &c2,const Vec3Expr<C3> &c3)
{
  return {{c1.expr,c2.expr,c3.expr}};
}

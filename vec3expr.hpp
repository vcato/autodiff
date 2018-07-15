#ifndef AUTODIFF_VEC3EXPR_HPP
#define AUTODIFF_VEC3EXPR_HPP



#include "dualvec3.hpp"
#include "scalarexpr.hpp"

namespace autodiff {

template <typename Expr>
struct Vec3Expr {
  Expr expr;
};


template <typename T> struct Vec3ExprTypeHelper;

template <typename T>
struct Vec3ExprTypeHelper<Vec3Expr<T>> {
  using type = T;
};


template <>
struct Vec3ExprTypeHelper<FloatVec3> {
  using type = FloatVec3;
};


template <>
struct Vec3ExprTypeHelper<DualVec3> {
  using type = DualVec3;
};


template <typename T>
using Vec3ExprType = typename Vec3ExprTypeHelper<T>::type;


template <>
struct Vec3Expr<DualVec3> {
  DualVec3 expr;

  ScalarExpr<DualFloat> x() const { return {expr.x()}; }
  ScalarExpr<DualFloat> y() const { return {expr.y()}; }
  ScalarExpr<DualFloat> z() const { return {expr.z()}; }
};


template <typename T> auto internal(Vec3Expr<T> e)
{
  return e.expr;
}


inline DualVec3 internal(const DualVec3& e)
{
  return e;
}


inline FloatVec3 internal(const FloatVec3& e)
{
  return e;
}


template <typename X,typename Y,typename Z>
struct Vec3Func {
  X x;
  Y y;
  Z z;
};


template <
  typename AExpr,
  typename BExpr,
  typename CExpr,
  typename A=ScalarExprType<AExpr>,
  typename B=ScalarExprType<BExpr>,
  typename C=ScalarExprType<CExpr>
>
static Vec3Expr<Vec3Func<A,B,C>>
  vec3(
    const AExpr &a,
    const BExpr &b,
    const CExpr &c
  )
{
  auto f = Vec3Func<A,B,C>{internal(a),internal(b),internal(c)};
  return Vec3Expr<Vec3Func<A,B,C>>{f};
}



template <typename T>
struct Evaluator<Vec3<T>> {
  Vec3<T> expr;

  Evaluator(const Vec3<T> &expr)
  : expr(expr)
  {
  }

  FloatVec3 value() const { return expr; }

  void addDeriv(const FloatVec3 &) { }
};


template <typename X,typename Y,typename Z>
struct Evaluator<Vec3Func<X,Y,Z>> {
  Evaluator<X> x_eval;
  Evaluator<Y> y_eval;
  Evaluator<Z> z_eval;

  Evaluator(const Vec3Func<X,Y,Z> &expr)
  : x_eval(expr.x),
    y_eval(expr.y),
    z_eval(expr.z)
  {
  }

  FloatVec3 value() const
  {
    return FloatVec3(x_eval.value(),y_eval.value(),z_eval.value());
  }

  void addDeriv(const FloatVec3 &deriv)
  {
    x_eval.addDeriv(deriv.x());
    y_eval.addDeriv(deriv.y());
    z_eval.addDeriv(deriv.z());
  }
};


inline FloatVec3 evaluate(const DualVec3 &dual)
{
  float x = evaluate(dual.x());
  float y = evaluate(dual.y());
  float z = evaluate(dual.z());
  return FloatVec3{x,y,z};
}


inline Vec3Expr<DualVec3> expr(FloatVec3 value,FloatVec3 &deriv)
{
  return {dual(value,deriv)};
}


template <typename A,typename B>
struct Dot {
  A a;
  B b;
};


template <
  typename AExpr,
  typename BExpr,
  typename A=Vec3ExprType<AExpr>,
  typename B=Vec3ExprType<BExpr>
>
ScalarExpr<Dot<A,B>> dot(const AExpr &a,const BExpr &b)
{
  return {{internal(a),internal(b)}};
}


template <>
struct Evaluator<DualVec3> {
  DualVec3 expr;

  Evaluator(const DualVec3 &arg)
  : expr(arg)
  {
  }

  FloatVec3 value() const { return evaluate(expr); }

  void addDeriv(const FloatVec3& deriv) const
  {
    expr.x().deriv += deriv.x();
    expr.y().deriv += deriv.y();
    expr.z().deriv += deriv.z();
  }
};


template <typename M>
struct Vec3ExprVar {
  Evaluator<M> eval;
  FloatVec3 _value = eval.value();
  mutable FloatVec3 deriv{0,0,0};

  FloatVec3 value() const { return _value; }
  DualVec3 dual() const { return autodiff::dual(_value,deriv); }

  Vec3ExprVar(const M &m)
  : eval(m)
  {
  }

  ~Vec3ExprVar()
  {
    eval.addDeriv(deriv);
  }
};


template <typename M>
struct Vec3ExprVar<Vec3Expr<M>> : Vec3ExprVar<M> {
  Vec3ExprVar(const Vec3Expr<M> &arg) : Vec3ExprVar<M>(arg.expr) { }
};


template <typename M>
struct Vec3ExprTypeHelper<Vec3ExprVar<M>> {
  using type = DualVec3;
};


template <typename M>
DualVec3 internal(const Vec3ExprVar<Vec3Expr<M>> &v)
{
  return v.dual();
}

template <typename M>
Vec3Expr<DualVec3> expr(Vec3ExprVar<M> &v)
{
  return {v.dual()};
}


template <typename A,typename B>
struct Evaluator<Dot<A,B>> {
  Vec3ExprVar<A> a;
  Vec3ExprVar<B> b;

  Evaluator(const Dot<A,B> &expr)
  : a(expr.a),
    b(expr.b)
  {
  }

  float value() const { return genDot(a.value(),b.value()); }

  void addDeriv(float dresult)
  {
    evalAndAddDeriv(genDot(expr(a),expr(b)),dresult);
  }
};


template <typename A,typename B>
struct Vec3Add {
  A a;
  B b;
};


template <
  typename AExpr,
  typename BExpr,
  typename A=Vec3ExprType<AExpr>,
  typename B=Vec3ExprType<BExpr>
>
Vec3Expr<Vec3Add<A,B>> operator+(const AExpr &a,const BExpr &b)
{
  return {{internal(a),internal(b)}};
}


template <typename A,typename B>
struct Evaluator<Vec3Add<A,B>> {
  Evaluator<A> a_eval;
  Evaluator<B> b_eval;

  Evaluator(const Vec3Add<A,B> &expr)
  : a_eval(expr.a),
    b_eval(expr.b)
  {
  }

  FloatVec3 value() const
  {
    return a_eval.value() + b_eval.value();
  }

  void addDeriv(const FloatVec3 &deriv)
  {
    a_eval.addDeriv(deriv);
    b_eval.addDeriv(deriv);
  }
};


template <typename A,typename B>
struct Vec3Sub {
  A a;
  B b;
};


template <
  typename AExpr,
  typename BExpr,
  typename A=Vec3ExprType<AExpr>,
  typename B=Vec3ExprType<BExpr>
>
Vec3Expr<Vec3Sub<A,B>> operator-(const AExpr &a,const BExpr &b)
{
  return {{internal(a),internal(b)}};
}


template <typename A,typename B>
struct Evaluator<Vec3Sub<A,B>> {
  Evaluator<A> a_eval;
  Evaluator<B> b_eval;

  Evaluator(const Vec3Sub<A,B> &expr)
  : a_eval(expr.a),
    b_eval(expr.b)
  {
  }

  FloatVec3 value() const
  {
    return a_eval.value() - b_eval.value();
  }

  void addDeriv(const FloatVec3 &deriv)
  {
    a_eval.addDeriv(deriv);
    b_eval.addDeriv(-deriv);
  }
};


template <typename V>
struct Vec3Mag {
  V v;
};


template <typename V>
struct Evaluator<Vec3Mag<V>> {
  Vec3ExprVar<V> v;
  Evaluator<decltype(sqrt(dot(expr(v),expr(v))).expr)> result_eval =
    sqrt(dot(expr(v),expr(v))).expr;

  Evaluator(const Vec3Mag<V> &expr)
  : v(expr.v)
  {
  }

  float value() const { return result_eval.value(); }

  void addDeriv(float dvalue)
  {
    result_eval.addDeriv(dvalue);
  }
};


template <typename VExpr,typename V=Vec3ExprType<VExpr>>
ScalarExpr<Vec3Mag<V>> mag(const VExpr &v)
{
  return {{internal(v)}};
}


inline ScalarExpr<Vec3Mag<DualVec3>> mag(const DualVec3 &v)
{
  return {{internal(v)}};
}


template <typename A,typename B>
struct Vec3Mul {
  A a;
  B b;
};


template <typename A,typename B>
struct Evaluator<Vec3Mul<A,B>> {
  Vec3ExprVar<A> a;
  ScalarExprVar<B> b;

  static auto result(Vec3ExprVar<A> &a,ScalarExprVar<B> &b) {
    return
      vec3(
        expr(a).x()*expr(b),
        expr(a).y()*expr(b),
        expr(a).z()*expr(b)
      ).expr;
  }

  Evaluator<decltype(result(a,b))> result_eval{result(a,b)};

  Evaluator(const Vec3Mul<A,B> &expr)
  : a(expr.a),
    b(expr.b)
  {
  }

  FloatVec3 value() const { return result_eval.value(); }

  void addDeriv(const FloatVec3 &deriv)
  {
    result_eval.addDeriv(deriv);
  }
};


template <
  typename AExpr,
  typename BExpr,
  typename A=Vec3ExprType<AExpr>,
  typename B=ScalarExprType<BExpr>
>
Vec3Expr<Vec3Mul<A,B>> operator*(const AExpr& a,const BExpr &b)
{
  return {{internal(a),internal(b)}};
}


template <
  typename AExpr,
  typename BExpr,
  typename A=Vec3ExprType<AExpr>,
  typename B=ScalarExprType<BExpr>
>
Vec3Expr<Vec3Mul<A,B>> operator*(const BExpr &b,const AExpr &a)
{
  return {{internal(a),internal(b)}};
}


template <typename A,typename B>
struct Vec3Div {
  A a;
  B b;
};


template <typename A,typename B>
struct Evaluator<Vec3Div<A,B>> {
  Vec3ExprVar<A> a;
  ScalarExprVar<B> b;

  static auto result(Vec3ExprVar<A> &a,ScalarExprVar<B> &b) {
    return
      vec3(
        expr(a).x()/expr(b),
        expr(a).y()/expr(b),
        expr(a).z()/expr(b)
      ).expr;
  }

  Evaluator<decltype(result(a,b))> result_eval{result(a,b)};

  Evaluator(const Vec3Div<A,B> &expr)
  : a(expr.a),
    b(expr.b)
  {
  }

  FloatVec3 value() const { return result_eval.value(); }

  void addDeriv(const FloatVec3 &deriv)
  {
    result_eval.addDeriv(deriv);
  }
};


template <
  typename AExpr,
  typename BExpr,
  typename A=Vec3ExprType<AExpr>,
  typename B=ScalarExprType<BExpr>
>
Vec3Expr<Vec3Div<A,B>> operator/(const AExpr &a,const BExpr& b)
{
  return {{internal(a),internal(b)}};
}


}


#endif /* AUTODIFF_VEC3EXPR_HPP */

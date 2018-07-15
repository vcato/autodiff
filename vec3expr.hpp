#include "dualvec3.hpp"
#include "scalarexpr.hpp"
#include "exprvar.hpp"

template <typename Expr>
struct Vec3Expr {
  Expr expr;
};


template <>
struct Vec3Expr<DualVec3> {
  DualVec3 expr;

  ScalarExpr<DualFloat> x() const { return {expr.x()}; }
  ScalarExpr<DualFloat> y() const { return {expr.y()}; }
  ScalarExpr<DualFloat> z() const { return {expr.z()}; }
};


template <typename E>
Vec3Expr<DualVec3> expr(ExprVar<Vec3Expr<E>> &v)
{
  return expr(v.value_member,v.deriv_member);
}


template <typename Expr>
struct ExprVar<Vec3Expr<Expr>> {
  Evaluator<Expr> eval;
  FloatVec3 value_member;
  FloatVec3 deriv_member;

  ExprVar(const Vec3Expr<Expr> &expr)
  : eval(expr.expr),
    value_member(eval.value()),
    deriv_member{0,0,0}
  {
  }

  FloatVec3 value() const { return value_member; }

  ~ExprVar()
  {
    eval.addDeriv(deriv_member);
  }
};


template <typename X,typename Y,typename Z>
struct Vec3Func {
  X x;
  Y y;
  Z z;
};


template <typename A,typename B,typename C>
static Vec3Expr<Vec3Func<A,B,C>>
  vec3(
    ScalarExpr<A> a,
    ScalarExpr<B> b,
    ScalarExpr<C> c
  )
{
  auto f = Vec3Func<A,B,C>{a.expr,b.expr,c.expr};
  return Vec3Expr<Vec3Func<A,B,C>>{f};
}


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


inline float dot(const FloatVec3 &a,const FloatVec3 &b)
{
  return genDot(a,b);
}


template <typename A,typename B>
struct Dot {
  A a;
  B b;
};


template <typename A,typename B>
ScalarExpr<Dot<A,B>> dot(Vec3Expr<A> a,Vec3Expr<B> b)
{
  return {{a.expr,b.expr}};
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
  FloatVec3 deriv{0,0,0};

  FloatVec3 value() const { return _value; }
  DualVec3 dual() { return ::dual(_value,deriv); }

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


template <typename A,typename B>
Vec3Expr<Vec3Add<A,B>> operator+(const Vec3Expr<A> &a,const Vec3Expr<B> &b)
{
  return {{a.expr,b.expr}};
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


template <typename A,typename B>
Vec3Expr<Vec3Sub<A,B>> operator-(const Vec3Expr<A> &a,const Vec3Expr<B> &b)
{
  return {{a.expr,b.expr}};
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


#if 1
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
#endif


template <typename V>
ScalarExpr<Vec3Mag<V>> mag(const Vec3Expr<V> &v)
{
  return {{v.expr}};
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


template <typename A,typename B>
Vec3Expr<Vec3Mul<A,B>> operator*(const Vec3Expr<A> &a,ScalarExpr<B> b)
{
  return {{a.expr,b.expr}};
}


template <typename A,typename B>
Vec3Expr<Vec3Mul<A,B>> operator*(ScalarExpr<B> b,const Vec3Expr<A> &a)
{
  return {{a.expr,b.expr}};
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


template <typename A,typename B>
Vec3Expr<Vec3Div<A,B>> operator/(const Vec3Expr<A> &a,ScalarExpr<B> b)
{
  return {{a.expr,b.expr}};
}

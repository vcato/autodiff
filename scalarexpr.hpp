#ifndef SCALAREXPR_HPP
#define SCALAREXPR_HPP


#include "evaluator.hpp"


template <typename T> struct ScalarExprTypeHelper;



template <typename Expr>
struct ScalarExpr {
  Expr expr;
};


inline ScalarExpr<DualFloat> expr(float value,float &deriv)
{
  return {{value,deriv}};
}


template <>
struct Evaluator<DualFloat> {
  DualFloat expr;

  Evaluator(const DualFloat &arg)
  : expr(arg)
  {
  }

  float value() const { return expr.value; }

  void addDeriv(float deriv) const
  {
    expr.deriv += deriv;
  }
};


template <typename M>
struct ScalarExprVar {
  Evaluator<M> eval;
  float _value = eval.value();
  float deriv = 0;

  float value() const { return _value; }
  void addDeriv(float dvalue) { deriv += dvalue; }
  DualFloat dual() { return ::dual(_value,deriv); }

  ScalarExprVar(const M &m)
  : eval(m)
  {
  }

  ~ScalarExprVar()
  {
    eval.addDeriv(deriv);
  }
};


template <typename M>
struct ScalarExprVar<ScalarExpr<M>> : ScalarExprVar<M> {
  ScalarExprVar(const ScalarExpr<M> &expr) : ScalarExprVar<M>(expr.expr) { }
};


// Specialization when the expression is just a DualFloat.  We don't
// need to store the value and the dual separately and we don't need
// to add the derivative in the destructor.
template <>
struct ScalarExprVar<DualFloat> {
  DualFloat expr;

  ScalarExprVar(const DualFloat &m)
  : expr(m)
  {
  }

  float value() const { return expr.value; }

  void addDeriv(float dvalue) { expr.deriv += dvalue; }

  DualFloat dual() { return expr; }
};


template <typename T>
struct ScalarExprTypeHelper<ScalarExprVar<T>&> {
  using type = DualFloat;
};


template <typename T>
DualFloat internal(ScalarExprVar<T> &v)
{
  return v.dual();
}


inline float evaluate(const DualFloat &dual)
{
  return dual.value;
}


template <typename M>
ScalarExpr<DualFloat> expr(ScalarExprVar<M> &v)
{
  return {v.dual()};
}


template <typename A,typename B>
struct ScalarAdd {
  A a;
  B b;
};


template <typename T> auto internal(ScalarExpr<T> e)
{
  return e.expr;
}


inline float internal(float e)
{
  return e;
}


inline DualFloat internal(const DualFloat &e)
{
  return e;
}


template <typename T>
struct ScalarExprTypeHelper<ScalarExpr<T>> {
  using type = T;
};


template <typename T>
struct ScalarExprTypeHelper<ScalarExpr<T>&> {
  using type = T;
};


template <typename T>
struct ScalarExprTypeHelper<const ScalarExpr<T> &> {
  using type = T;
};


template <>
struct ScalarExprTypeHelper<float> {
  using type = float;
};


template <>
struct ScalarExprTypeHelper<DualFloat> {
  using type = DualFloat;
};


template <>
struct ScalarExprTypeHelper<DualFloat&> {
  using type = DualFloat;
};


template <>
struct ScalarExprTypeHelper<int> {
  using type = float;
};


template <typename T>
using ScalarExprType = typename ScalarExprTypeHelper<T>::type;

template <
  typename AExpr,
  typename BExpr,
  typename A = ScalarExprType<AExpr>,
  typename B = ScalarExprType<BExpr>
>
ScalarExpr<ScalarAdd<A,B>>
  operator+(AExpr &&a,BExpr &&b)
{
  return {{internal(a),internal(b)}};
}


template <typename A,typename B>
struct ScalarSub {
  A a;
  B b;
};


template <
  typename AExpr,
  typename BExpr,
  typename A=ScalarExprType<AExpr>,
  typename B=ScalarExprType<BExpr>
>
ScalarExpr<ScalarSub<A,B>>
  operator-(AExpr &&a,BExpr &&b)
{
  return {{internal(a),internal(b)}};
}


template <typename A,typename B>
struct ScalarMul {
  A a;
  B b;
};


template <
  typename AExpr,
  typename BExpr,
  typename A=ScalarExprType<AExpr>,
  typename B=ScalarExprType<BExpr>
>
ScalarExpr<ScalarMul<A,B>>
  operator*(AExpr &&a,BExpr &&b)
{
  return {{internal(a),internal(b)}};
}


template <typename A,typename B>
struct ScalarDiv {
  A a;
  B b;
};


template <
  typename AExpr,
  typename BExpr,
  typename A=ScalarExprType<AExpr>,
  typename B=ScalarExprType<BExpr>
>
ScalarExpr<ScalarDiv<A,B>> operator/(AExpr &&a,BExpr &&b)
{
  return {{internal(a),internal(b)}};
}


template <typename A>
struct Cos {
  A a;
};


template <typename AExpr,typename A=ScalarExprType<AExpr>>
inline ScalarExpr<Cos<A>> cos(AExpr &&a)
{
  return {{internal(a)}};
}


template <typename A>
struct Evaluator<Cos<A>> {
  ScalarExprVar<A> a;

  Evaluator(const Cos<A> &expr)
  : a(expr.a)
  {
  }

  float value() const { return cosf(a.value()); }

  void addDeriv(float dvalue)
  {
    a.addDeriv(dvalue*-sinf(a.value()));
  }
};


template <typename A>
struct Sin {
  A a;
};


template <typename AExpr,typename A=ScalarExprType<AExpr>>
ScalarExpr<Sin<A>> sin(AExpr &&a)
{
  return {{internal(a)}};
}


template <typename A>
struct Evaluator<Sin<A>> {
  ScalarExprVar<A> a;

  Evaluator(const Sin<A> &expr)
  : a(expr.a)
  {
  }

  float value() const { return sinf(a.value()); }

  void addDeriv(float dvalue)
  {
    a.addDeriv(dvalue*cosf(a.value()));
  }
};


template <>
struct Evaluator<float> {
  float expr;

  Evaluator(float arg)
  : expr(arg)
  {
  }

  float value() const { return expr; }

  void addDeriv(float) const
  {
  }
};


template <typename A,typename B>
struct Evaluator<ScalarAdd<A,B>> {
  Evaluator<A> a_eval;
  Evaluator<B> b_eval;

  Evaluator(ScalarAdd<A,B> expr)
  : a_eval(expr.a), b_eval(expr.b)
  {
  }

  float value() const
  {
    return a_eval.value() + b_eval.value();
  }

  void addDeriv(float deriv)
  {
    a_eval.addDeriv(deriv);
    b_eval.addDeriv(deriv);
  }
};


template <typename A,typename B>
struct Evaluator<ScalarSub<A,B>> {
  Evaluator<A> a;
  Evaluator<B> b;

  Evaluator(ScalarSub<A,B> expr)
  : a(expr.a), b(expr.b)
  {
  }

  float value() const
  {
    return a.value() - b.value();
  }

  void addDeriv(float deriv)
  {
    a.addDeriv(deriv);
    b.addDeriv(-deriv);
  }
};


template <typename A,typename B>
struct Evaluator<ScalarMul<A,B>> {
  Evaluator<A> a_eval;
  float a = a_eval.value();
  Evaluator<B> b_eval;
  float b = b_eval.value();

  Evaluator(ScalarMul<A,B> expr)
  : a_eval(expr.a),
    b_eval(expr.b)
  {
  }

  float value() const
  {
    return a*b;
  }

  void addDeriv(float deriv)
  {
    a_eval.addDeriv(deriv*b);
    b_eval.addDeriv(deriv*a);
  }
};


template <typename A,typename B>
struct Evaluator<ScalarDiv<A,B>> {
  Evaluator<A> a_eval;
  float a = a_eval.value();
  Evaluator<B> b_eval;
  float b = b_eval.value();

  Evaluator(ScalarDiv<A,B> expr)
  : a_eval(expr.a),
    b_eval(expr.b)
  {
  }

  float value() const
  {
    return a/b;
  }

  void addDeriv(float deriv)
  {
    a_eval.addDeriv(deriv* b/(b*b));
    b_eval.addDeriv(deriv*-a/(b*b));
  }
};


template <typename X>
struct Sqrt {
  X x;
};


template <typename X>
struct Evaluator<Sqrt<X>> {
  Evaluator<X> x_eval;
  float sqrt_x = sqrtf(x_eval.value());

  Evaluator(const Sqrt<X> &expr)
  : x_eval(expr.x)
  {
  }

  float value() const { return sqrt_x; }

  void addDeriv(float dvalue)
  {
    x_eval.addDeriv(dvalue*0.5/sqrt_x);
  }
};


template <typename XExpr,typename X=ScalarExprType<XExpr>>
ScalarExpr<Sqrt<X>> sqrt(XExpr &&x)
{
  return {{internal(x)}};
}



inline ScalarExpr<float> expr(float x)
{
  return {x};
}


#endif /* SCALAREXPR_HPP */

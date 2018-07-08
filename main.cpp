#include <iostream>
#include <cmath>
#include <cfloat>
#include "random.hpp"


using std::cerr;
using std::ostream;
using std::max;


#define assertNear(actual,expected,tolerance) \
  (assertNearHelper(actual,expected,tolerance,__FILE__,__LINE__))



namespace {
template <typename X,typename Y,typename Z>
struct Vec3Func {
  X x;
  Y y;
  Z z;
};
}


namespace {
struct DualFloat {
  float value;
  float &deriv;
};
}


namespace {
static float evaluate(const DualFloat &dual)
{
  return dual.value;
}
}


namespace {
static void addDeriv(const DualFloat &dual,float dresult)
{
  dual.deriv += dresult;
}
}


namespace {
template <typename Expr>
struct ScalarExpr {
  Expr expr;
};
}


namespace {
template <typename Expr>
static float evaluate(const ScalarExpr<Expr> &expr)
{
  return evaluate(expr.expr);
}
}


namespace {
template <typename Self>
struct RowRef {
  Self &self;
  int i;

  decltype(auto) operator[](int j)
  {
    return self.element(self,i,j);
  }
};
}


namespace {
template <typename Self>
static auto row(Self &arg,int i)
{
  return RowRef<Self>{arg,i};
}
}

namespace {
template <typename T>
struct Mat33 {
  Mat33(const T (&arg)[3][3])
  : values{
      {arg[0][0],arg[0][1],arg[0][2]},
      {arg[1][0],arg[1][1],arg[1][2]},
      {arg[2][0],arg[2][1],arg[2][2]},
    }
  {
  }

  Mat33 &operator+=(Mat33 arg)
  {
    for (int i=0; i!=3; ++i) {
      for (int j=0; j!=3; ++j) {
        values[i][j] += arg.values[i][j];
      }
    }

    return *this;
  }

  friend ostream& operator<<(ostream &stream,const Mat33 &arg)
  {
    stream << "\n";

    for (int i=0; i!=3; ++i) {
      for (int j=0; j!=3; ++j) {
        stream << "  ";
        stream << arg.values[i][j];
      }
      stream << "\n";
    }

    return stream;
  }

  template <typename Self>
  static decltype(auto) element(Self &arg,int i,int j)
  {
    assert(i>=0);
    assert(i<3);
    assert(j>=0);
    assert(j<3);
    return arg.values[i][j];
  }

  auto operator[](int i)             { return row(*this,i); }
  auto operator[](int i) const       { return row(*this,i); }

  T values[3][3];
};
}


static float differenceBetween(float a,float b)
{
  return fabsf(a-b);
}


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


using FloatMat33 = Mat33<float>;
using DualMat33 = Mat33<DualFloat>;


namespace {
template <typename T>
static void addDeriv(const Mat33<T> &m,const FloatMat33 &dm)
{
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      addDeriv(m[i][j],dm[i][j]);
    }
  }
}
}


namespace {
template <typename A,typename B>
struct ScalarAdd {
  A a;
  B b;
};
}


namespace {
template <typename A,typename B>
static float evaluate(const ScalarAdd<A,B> &expr)
{
  return evaluate(expr.a) + evaluate(expr.b);
}
}


template <typename A,typename B>
static ScalarExpr<ScalarAdd<A,B>>
  operator+(const ScalarExpr<A> &a,const ScalarExpr<B> &b)
{
  return {{a.expr,b.expr}};
}


namespace {
template <typename A,typename B>
struct ScalarSub {
  A a;
  B b;
};
}


namespace {
template <typename A,typename B>
static float evaluate(const ScalarSub<A,B> &expr)
{
  return evaluate(expr.a) - evaluate(expr.b);
}
}


template <typename A,typename B>
static ScalarExpr<ScalarSub<A,B>>
  operator-(const ScalarExpr<A> &a,const ScalarExpr<B> &b)
{
  return {{a.expr,b.expr}};
}


namespace {
template <typename A,typename B>
struct ScalarMul {
  A a;
  B b;
};
}


namespace {
template <typename A,typename B>
static float evaluate(const ScalarMul<A,B> &expr)
{
  return evaluate(expr.a)*evaluate(expr.b);
}
}


namespace {
template <typename A,typename B>
static void addDeriv(const ScalarMul<A,B> &expr,float dresult)
{
  addDeriv(expr.a, dresult*evaluate(expr.b));
  addDeriv(expr.b, dresult*evaluate(expr.a));
}
}


template <typename A,typename B>
static ScalarExpr<ScalarMul<A,B>>
  operator*(const ScalarExpr<A> &a,const ScalarExpr<B> &b)
{
  return {{a.expr,b.expr}};
}


static FloatMat33 mat33All(float arg)
{
  float values[3][3] = {
    {arg,arg,arg},
    {arg,arg,arg},
    {arg,arg,arg},
  };

  return FloatMat33(values);
}


template <typename T> static T zero();

template <>
Mat33<float> zero<Mat33<float>>()
{
  return mat33All(0);
}


template <>
float zero<float>()
{
  return 0;
}


static void ddiv(DualFloat a_arg,DualFloat b_arg,float dresult)
{
  // d(a/b) = (b*da - a*db)/(b*b);
  float b = evaluate(b_arg);
  float a = evaluate(a_arg);
  addDeriv(a_arg,dresult* b/(b*b));
  addDeriv(b_arg,dresult*-a/(b*b));
}


static DualFloat dual(float av,float &da)
{
  return {av,da};
}


namespace {
template <typename T>
struct Vec3 {
  T _x, _y, _z;

  T x() const { return _x; }
  T y() const { return _y; }
  T z() const { return _z; }

  T& x() { return _x; }
  T& y() { return _y; }
  T& z() { return _z; }

  friend ostream& operator<<(ostream &stream,const Vec3 &self)
  {
    stream << "[" << self.x() << "," << self.y() << "," << self.z() << "]";
    return stream;
  }
};
}


using DualVec3 = Vec3<DualFloat>;
using FloatVec3 = Vec3<float>;


namespace {
static FloatVec3 evaluate(const DualVec3 &dual)
{
  float x = evaluate(dual.x());
  float y = evaluate(dual.y());
  float z = evaluate(dual.z());
  return FloatVec3{x,y,z};
}
}


namespace {
template <typename X,typename Y,typename Z>
static void addDeriv(Vec3Func<X,Y,Z> f,const FloatVec3 &deriv)
{
  addDeriv(f.x,deriv.x);
  addDeriv(f.y,deriv.y);
  addDeriv(f.z,deriv.z);
}
}


namespace {
template <typename X,typename Y,typename Z>
static FloatVec3 evaluate(Vec3Func<X,Y,Z> f)
{
  float x = evaluate(f.x);
  float y = evaluate(f.y);
  float z = evaluate(f.z);
  return {x,y,z};
}
}


static DualVec3 dual(FloatVec3 v,FloatVec3 &dv)
{
  return {{v.x(),dv.x()},{v.y(),dv.y()},{v.z(),dv.z()}};
}


static DualMat33 dual(Mat33<float> v,Mat33<float> &dv)
{
  DualFloat values[3][3] = {
    {dual(v[0][0],dv[0][0]),dual(v[0][1],dv[0][1]),dual(v[0][2],dv[0][2])},
    {dual(v[1][0],dv[1][0]),dual(v[1][1],dv[1][1]),dual(v[1][2],dv[1][2])},
    {dual(v[2][0],dv[2][0]),dual(v[2][1],dv[2][1]),dual(v[2][2],dv[2][2])},
  };

  return DualMat33(values);
}


namespace {
template <typename Expr>
struct Mat33Expr {
  Expr expr;
};
}


namespace {
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
}


namespace {
template <typename A,typename B>
struct ScalarDiv {
  A a;
  B b;
};
}


namespace {
template <typename A,typename B>
struct Mat33Div {
  A a;
  B b;
};
}


namespace {
template <typename A,typename B>
static Mat33Expr<Mat33Div<A,B>> operator/(Mat33Expr<A> a,ScalarExpr<B> b)
{
  return {{a.expr,b.expr}};
}
}


namespace {
template <typename A,typename B>
static ScalarExpr<ScalarDiv<A,B>> operator/(ScalarExpr<A> a,ScalarExpr<B> b)
{
  return {{a.expr,b.expr}};
}
}


static ScalarExpr<DualFloat> expr(const DualFloat &expr)
{
  return {{expr}};
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
static FloatMat33 finiteDeriv(Function f,FloatMat33 &m)
{
  FloatMat33 result = mat33All(0);

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


template <typename Function>
static FloatVec3 finiteDeriv(Function f,FloatVec3 &v)
{
  float dx = finiteDeriv(f,v.x());
  float dy = finiteDeriv(f,v.y());
  float dz = finiteDeriv(f,v.z());
  return FloatVec3{dx,dy,dz};
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


namespace {
template <typename Expr>
struct Vec3Expr {
  Expr expr;
};
}


namespace {
template <>
struct Vec3Expr<DualVec3> {
  DualVec3 expr;

  ScalarExpr<DualFloat> x() const { return {expr.x()}; }
  ScalarExpr<DualFloat> y() const { return {expr.y()}; }
  ScalarExpr<DualFloat> z() const { return {expr.z()}; }
};
}


static Vec3Expr<DualVec3> expr(const DualVec3 &expr)
{
  return {{expr}};
}


template <typename A,typename B>
static auto genDot(const A a,const B b)
{
  return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
}


static float dot(const FloatVec3 &a,const FloatVec3 &b)
{
  return genDot(a,b);
}


namespace {
template <typename A,typename B>
struct Dot {
  A a;
  B b;
};
}


template <typename A,typename B>
static ScalarExpr<Dot<A,B>> dot(Vec3Expr<A> a,Vec3Expr<B> b)
{
  return {{a.expr,b.expr}};
}


namespace {
template <typename Expr> struct Evaluator;
}


namespace {
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
}


namespace {
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
}


namespace {
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
    ::addDeriv(expr.x(),deriv.x());
    ::addDeriv(expr.y(),deriv.y());
    ::addDeriv(expr.z(),deriv.z());
  }
};
}


namespace {
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
}


namespace {
template <typename M>
struct ScalarExprVar {
  Evaluator<M> eval;
  float _value = eval.value();
  float value() const { return _value; }
  float deriv = 0;
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
}


namespace {
template <typename A,typename B>
struct Evaluator<ScalarAdd<A,B>> {
  ScalarExprVar<A> a;
  ScalarExprVar<B> b;

  Evaluator(ScalarAdd<A,B> expr)
  : a(expr.a), b(expr.b)
  {
  }

  float value() const
  {
    return a.value() + b.value();
  }

  void addDeriv(float deriv)
  {
    ::addDeriv(a.dual(),deriv);
    ::addDeriv(b.dual(),deriv);
  }
};
}


namespace {
template <typename A,typename B>
struct Evaluator<ScalarSub<A,B>> {
  ScalarExprVar<A> a;
  ScalarExprVar<B> b;

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
    ::addDeriv(a.dual(),deriv);
    ::addDeriv(b.dual(),-deriv);
  }
};
}


namespace {
template <typename A,typename B>
struct Evaluator<ScalarMul<A,B>> {
  ScalarExprVar<A> a;
  ScalarExprVar<B> b;

  Evaluator(ScalarMul<A,B> expr)
  : a(expr.a),
    b(expr.b)
  {
  }

  float value() const
  {
    return a.value() * b.value();
  }

  void addDeriv(float deriv)
  {
    ::addDeriv(a.dual(),deriv*b.value());
    ::addDeriv(b.dual(),deriv*a.value());
  }
};
}


namespace {
template <typename A,typename B>
struct Evaluator<ScalarDiv<A,B>> {
  Evaluator<A> a_eval;
  Evaluator<B> b_eval;
  float a = a_eval.value();
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
    float da = 0;
    float db = 0;
    ddiv(dual(a,da),dual(b,db),deriv);
    a_eval.addDeriv(da);
    b_eval.addDeriv(db);
  }
};
}


namespace {
template <typename M>
struct Cofactor {
  M m;
  int i;
  int j;
};
}


template <typename Mat>
static auto genCofactor(const Mat &arg,int i,int j)
{
  auto a11 = arg[(i+1)%3][(j+1)%3];
  auto a12 = arg[(i+1)%3][(j+2)%3];
  auto a21 = arg[(i+2)%3][(j+1)%3];
  auto a22 = arg[(i+2)%3][(j+2)%3];

  return a11*a22 - a12*a21;
}


static float cofactor(const FloatMat33 &arg,int i,int j)
{
  return genCofactor(arg,i,j);
}


// Specialization when the expression is just a DualFloat.  We don't
// need to store the value and the dual separately and we don't need
// to add the derivative in the destructor.
namespace {
template <>
struct ScalarExprVar<DualFloat> {
  DualFloat expr;

  ScalarExprVar(const DualFloat &m)
  : expr(m)
  {
  }

  float value() const { return evaluate(expr); }

  DualFloat dual() { return expr; }
};
}


namespace {
template <typename M>
struct Vec3ExprVar {
  Evaluator<M> eval;
  FloatVec3 _value = eval.value();
  FloatVec3 value() const { return _value; }
  FloatVec3 deriv{0,0,0};
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
}


namespace {
template <typename M>
struct Mat33ExprVar {
  Evaluator<M> eval;
  FloatMat33 _value = eval.value();
  FloatMat33 value() const { return _value; }
  FloatMat33 deriv = mat33All(0);
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
}


template <typename ExprWrapper,typename Value>
static Value
  evalAndAddDeriv(const ExprWrapper &mat33_expr,const Value& dresult)
{
  auto e = mat33_expr.expr;
  Evaluator<decltype(e)> eval(e);
  Value result = eval.value();
  eval.addDeriv(dresult);
  return result;
}


namespace {
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
    evalAndAddDeriv(genDot(expr(a.dual()),expr(b.dual())),dresult);
  }
};
}



// Specialization in the case of a DualMat33 expression, where we don't
// need to keep the value separately, and we don't need to add the
// derivative in the destructor.
namespace {
template <>
struct Mat33ExprVar<DualMat33> {
  DualMat33 expr;

  Mat33ExprVar(const DualMat33 &e)
  : expr(e)
  {
  }

  FloatMat33 value() const { return evaluate(expr); }

  DualMat33 dual() { return expr; }
};
}


template <typename T>
static Mat33Expr<Mat33<T>> expr(const Mat33<T> &expr)
{
  return {expr};
}


namespace {
template <typename M>
struct Evaluator<Cofactor<M>> {
  Mat33ExprVar<M> m;
  Evaluator<decltype(genCofactor(expr(m.dual()),0,0).expr)> result;

  Evaluator(const Cofactor<M> &c)
  : m{c.m},
    result{genCofactor(expr(m.dual()),c.i,c.j).expr}
  {
  }

  float value() const { return result.value(); }

  void addDeriv(float dresult) { result.addDeriv(dresult); }
};
}


namespace {
template <typename A,typename B>
struct Evaluator<Mat33Div<A,B>> {
  Mat33ExprVar<A> a;
  ScalarExprVar<B> b;

  Evaluator(const Mat33Div<A,B> &expr)
  : a(expr.a),
    b(expr.b)
  {
  }

  FloatMat33 value() const
  {
    return a.value() / b.value();
  }

  void addDeriv(const FloatMat33 &deriv)
  {
    for (int i=0; i!=3; ++i) {
      for (int j=0; j!=3; ++j) {
        ddiv(a.dual()[i][j],b.dual(),deriv[i][j]);
      }
    }
  }
};
}


template <typename T>
static auto transpose(const Mat33<T> &a)
{
  T values[3][3] = {
    {a[0][0],a[1][0],a[2][0]},
    {a[0][1],a[1][1],a[2][1]},
    {a[0][2],a[1][2],a[2][2]},
  };

  return Mat33<T>(values);
}


namespace {
template <typename M>
struct Transpose {
  M m;
};
}


static void dtranspose(const DualMat33 &a,const FloatMat33 &dresult)
{
  addDeriv(a,transpose(dresult));
}


namespace {
template <typename M>
static void addDeriv(Transpose<M> transpose,const FloatMat33 &deriv)
{
  dtranspose(transpose.m,deriv);
}
}


template <typename M>
static Mat33Expr<Transpose<M>> transpose(const Mat33Expr<M> &m)
{
  return {{m.expr}};
}

namespace {
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
    dtranspose(m.dual(),dresult);
  }
};
}


template <typename M>
static ScalarExpr<Cofactor<M>>
  cofactor(const Mat33Expr<M> &arg,int i,int j)
{
  return {{arg.expr,i,j}};
}


template <typename T>
static Mat33<T> mat33(const T (&values)[3][3])
{
  return Mat33<T>(values);
}


template <typename T>
static Mat33Expr<Mat33<T>> mat33(const ScalarExpr<T> (&arg)[3][3])
{
  T values[3][3] = {
    {arg[0][0].expr,arg[0][1].expr,arg[0][2].expr},
    {arg[1][0].expr,arg[1][1].expr,arg[1][2].expr},
    {arg[2][0].expr,arg[2][1].expr,arg[2][2].expr},
  };

  return {{values}};
}


template <typename Mat>
static auto genCofactorMatrix(const Mat &arg)
{
  using V = decltype(genCofactor(arg,0,0));

  V values[3][3] = {
    {genCofactor(arg,0,0),genCofactor(arg,0,1),genCofactor(arg,0,2)},
    {genCofactor(arg,1,0),genCofactor(arg,1,1),genCofactor(arg,1,2)},
    {genCofactor(arg,2,0),genCofactor(arg,2,1),genCofactor(arg,2,2)},
  };

  return mat33(values);
}


static auto cofactorMatrix(const FloatMat33 &arg)
{
  return genCofactorMatrix(arg);
}


namespace {
template <typename M>
struct CofactorMatrixFunc {
  M m;
};
}


template <typename M>
static Mat33Expr<CofactorMatrixFunc<M>> cofactorMatrix(const Mat33Expr<M> &m)
{
  return {{m.expr}};
}


namespace {
template <typename M>
struct Evaluator<CofactorMatrixFunc<M>> {
  Mat33ExprVar<M> m;
  Evaluator<decltype(genCofactorMatrix(expr(m.dual())).expr)> result_eval =
    genCofactorMatrix(expr(m.dual())).expr;

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
}


template <typename Mat>
static auto genDeterminant(const Mat &a)
{
  return
    a[0][0] * genCofactor(a,0,0) +
    a[0][1] * genCofactor(a,0,1) +
    a[0][2] * genCofactor(a,0,2);
}


static auto determinant(const FloatMat33 &a)
{
  return genDeterminant(a);
}



namespace {
template <typename M>
struct Determinant {
  M m;
};
}


template <typename M>
static ScalarExpr<Determinant<M>> determinant(const Mat33Expr<M> &m)
{
  const M &expr = m.expr;
  return {{expr}};
}


namespace {
template <typename M>
struct Evaluator<Determinant<M>> {
  Mat33ExprVar<M> m;
  Evaluator<decltype(genDeterminant(expr(m.dual())).expr)>
    result{genDeterminant(expr(m.dual())).expr};

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
}


namespace {
static FloatMat33 operator/(const FloatMat33 &a,float b)
{
  float values[3][3];

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      values[i][j] = a[i][j]/b;
    }
  }

  return FloatMat33(values);
}
}


static FloatMat33 operator*(const FloatMat33 &a,const FloatMat33 &b)
{
  float values[3][3];

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      float sum = 0;

      for (int k=0; k!=3; ++k) {
        sum += a[i][k]*b[k][j];
      }

      values[i][j] = sum;
    }
  }

  return FloatMat33(values);
}


template <typename T>
static auto genMat33Inv(const T &arg)
{
  return transpose(cofactorMatrix(arg))/determinant(arg);
}


template <typename T>
static auto mat33Inv(const Mat33<T> &arg)
{
  return genMat33Inv(arg);
}


namespace {
template <typename M>
struct Inv {
  M m;
};
}


template <typename T>
static Mat33Expr<Inv<T>> mat33Inv(const Mat33Expr<T> &arg)
{
  return {{arg.expr}};
}


namespace {
template <typename M>
struct Evaluator<Inv<M>>
{
  Mat33ExprVar<M> m;
  Evaluator<
    decltype(genMat33Inv(expr(m.dual())).expr)
  > result_eval{genMat33Inv(expr(m.dual())).expr};

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
}


static FloatMat33 mat33Identity()
{
  float values[3][3] = {
    {1,0,0},
    {0,1,0},
    {0,0,1},
  };

  return FloatMat33(values);
}


static ScalarExpr<float> expr(float x)
{
  return {x};
}


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

namespace tests {

static void testMat33Inv()
{
  {
    FloatMat33 a = mat33Identity();
    FloatMat33 a_inv = mat33Inv(a);
    FloatMat33 a_times_a_inv = a*a_inv;
    assertNear(a_times_a_inv,mat33Identity(),1e-4);
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
    float result = evalAndAddDeriv(expr(dual(a,da)) + expr(dual(b,db)),dresult);

    assert(da==1);
    assert(db==1);
    assert(result==3);
  }
  {
    float a = 1;
    float da = 0;
    float dresult = 1;
    float result = evalAndAddDeriv(expr(dual(a,da)) + expr(2),dresult);

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
  float result = evalAndAddDeriv(expr(dual(a,da)) - expr(dual(b,db)),dresult);

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
  float result = evalAndAddDeriv(expr(dual(a,da)) * expr(dual(b,db)),dresult);

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
  float result = evalAndAddDeriv(expr(dual(a,da)) / expr(dual(b,db)),dresult);
  auto f = [&]{ return a/b; };

  assertNear(da,finiteDeriv(f,a),1e-4);
  assertNear(db,finiteDeriv(f,b),1e-4);
  assert(result==a/b);
}


static void testDualVec3Evaluator()
{
  FloatVec3 a{1,2,3};
  FloatVec3 da{0,0,0};
  FloatVec3 dresult{1,0,0};
  evalAndAddDeriv(expr(dual(a,da)),dresult);
  auto f = [&]{ return a.x(); };
  assertNear(da,finiteDeriv(f,a),1e-4);
}


static void testCofactorEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      FloatMat33 dmat = mat33All(0);
      DualMat33 dual_mat = dual(mat,dmat);
      float dresult = 1;
      float result = evalAndAddDeriv(cofactor(expr(dual_mat),i,j),dresult);

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
  FloatMat33 dmat = mat33All(0);
  DualMat33 m = dual(mat,dmat);
  using Element = ScalarAdd<DualFloat,float>;
  Element values[3][3] = {
    {(expr(m[0][0]) + expr(1)).expr,
     (expr(m[0][1]) + expr(1)).expr,
     (expr(m[0][2]) + expr(1)).expr},
    {(expr(m[1][0]) + expr(1)).expr,
     (expr(m[1][1]) + expr(1)).expr,
     (expr(m[1][2]) + expr(1)).expr},
    {(expr(m[2][0]) + expr(1)).expr,
     (expr(m[2][1]) + expr(1)).expr,
     (expr(m[2][2]) + expr(1)).expr},
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

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      FloatMat33 dmat = mat33All(0);
      FloatMat33 dresult = mat33All(0);
      dresult[i][j] = 1;

      FloatMat33 result =
        evalAndAddDeriv(cofactorMatrix(expr(dual(mat,dmat))),dresult);
      auto f = [&]{ return cofactorMatrix(mat)[i][j]; };
      assert(result[i][j]==f());
      assertNear(dmat,finiteDeriv(f,mat),1e-4);
    }
  }
}


static void testDeterminantEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  FloatMat33 dmat = mat33All(0);
  float dresult = 1;
  float result = evalAndAddDeriv(determinant(expr(dual(mat,dmat))),dresult);
  auto f = [&]{ return determinant(mat); };
  assert(result==f());
  assertNear(dmat,finiteDeriv(f,mat),1e-4);
}


static void testTransposeEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      FloatMat33 dmat = mat33All(0);
      FloatMat33 dresult = mat33All(0);
      dresult[i][j] = 1;
      FloatMat33 result =
        evalAndAddDeriv(transpose(expr(dual(mat,dmat))),dresult);
      auto f = [&]{ return transpose(mat)[i][j]; };
      assert(result[i][j]==f());
      assertNear(dmat,finiteDeriv(f,mat),1e-4);
    }
  }
}


static void testMat33DivEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);
  float divisor = randomFloat(-1,1,random_engine);
  float ddivisor = 0;

  FloatMat33 dmat = mat33All(0);
  FloatMat33 dresult = mat33All(0);
  dresult[0][0] = 1;
  auto e = expr(dual(mat,dmat))/expr(dual(divisor,ddivisor));
  FloatMat33 result = evalAndAddDeriv(e,dresult);
  auto f = [&]{ return (mat/divisor)[0][0]; };
  assert(result[0][0]==f());
  assertNear(dmat,finiteDeriv(f,mat),1e-4);
}


static void testMat33InvEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 mat = randomMat33(random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      FloatMat33 dmat = mat33All(0);
      FloatMat33 dresult = mat33All(0);
      dresult[i][j] = 1;
      FloatMat33 result =
        evalAndAddDeriv(mat33Inv(expr(dual(mat,dmat))),dresult);
      auto f = [&]{ return mat33Inv(mat)[i][j]; };
      assert(result[i][j]==f());
      assertNear(dmat,finiteDeriv(f,mat),2e-4);
    }
  }
}


static void testDotEvaluator()
{
  RandomEngine random_engine(/*seed*/1);
  FloatVec3 a = randomVec3(random_engine);
  FloatVec3 b = randomVec3(random_engine);
  FloatVec3 da{0,0,0};
  FloatVec3 db{0,0,0};

  auto e = dot(expr(dual(a,da)),expr(dual(b,db)));
  float dresult = 1;
  float result = evalAndAddDeriv(e,dresult);
  auto f = [&]{ return dot(a,b); };
  assert(result==f());
  assertNear(da,finiteDeriv(f,a),1e-4);
  assertNear(db,finiteDeriv(f,b),1e-4);
}


#if 0
static void testDot2()
{
  float dv = 0;
  {
    DualFloat v{3,dv};
    auto a = var(vec3(expr(2),expr(v),expr(4)));
    auto b = var(vec3(expr(5),expr(6),expr(7)));
    auto e = var(dot(a,b));
    addDeriv(e,1);
  }
  assertNear(dv,6.0f,1e-4);
}
#endif


#if 0
static void testAddingVectors()
{
  float dv = 0;

  {
    DualFloat v{3,dv};
    auto a = vec3(1,2,v);
    auto b = vec3(4,5,6);
    auto e = var(a+b);
    addDeriv(e,1);
  }

  assertNear(dv,1.0f,1e-4);
}
#endif


}


int main()
{
  tests::testMat33Inv();
  tests::testScalarConstantEvaluator();
  tests::testScalarAddEvaluator();
  tests::testScalarSubEvaluator();
  tests::testScalarMulEvaluator();
  tests::testScalarDivEvaluator();
  tests::testDualVec3Evaluator();
  tests::testCofactorEvaluator();
  tests::testMat33Evaluator();
  tests::testCofactorMatrixEvaluator();
  tests::testDeterminantEvaluator();
  tests::testTransposeEvaluator();
  tests::testMat33DivEvaluator();
  tests::testMat33InvEvaluator();
  tests::testDotEvaluator();
  //tests::testDot2();
  // tests::testAddingVectors();
}

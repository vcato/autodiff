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
static float evaluate(float arg)
{
  return arg;
}
}


namespace {
static void addDeriv(const float &,float)
{
  // Nothing to do if we are adding a derivative to a constant.
}
}


namespace {
template <typename X,typename Y,typename Z>
struct VectorFunc {
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
template <typename T>
struct Mat33 {
  Mat33(T (&arg)[3][3])
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
  struct RowRef {
    Self &self;
    int i;

    auto &operator[](int j)
    {
      return self(i,j);
    }
  };

  template <typename Self>
  static auto row(Self &arg,int i)
  {
    return RowRef<Self>{arg,i};
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
  decltype(auto) operator()(int i,int j)       { return element(*this,i,j); }
  decltype(auto) operator()(int i,int j) const { return element(*this,i,j); }

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


namespace {
template <typename Expr>
static void addDeriv(const ScalarExpr<Expr> &scalar_expr,float deriv)
{
  addDeriv(scalar_expr.expr,deriv);
}
}


namespace {
template <typename T>
static void addDeriv(const Mat33<T> &m,const Mat33<float> &dm)
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
struct Add {
  A a;
  B b;
};
}


namespace {
template <typename A,typename B>
static float evaluate(const Add<A,B> &expr)
{
  return evaluate(expr.a) + evaluate(expr.b);
}
}


namespace {
template <typename A,typename B>
static void addDeriv(const Add<A,B> &expr,float dresult)
{
  addDeriv(expr.a,dresult);
  addDeriv(expr.b,dresult);
}
}


template <typename AExpr,typename BExpr>
static Add<AExpr,BExpr> operator+(const AExpr &a,const BExpr &b)
{
  return {a,b};
}


template <typename A,typename B>
static ScalarExpr<Add<A,B>>
  operator+(const ScalarExpr<A> &a,const ScalarExpr<B> &b)
{
  return {{a.expr,b.expr}};
}


namespace {
template <typename A,typename B>
struct Sub {
  A a;
  B b;
};
}


namespace {
template <typename A,typename B>
static float evaluate(const Sub<A,B> &expr)
{
  return evaluate(expr.a) - evaluate(expr.b);
}
}


namespace {
template <typename A,typename B>
static void addDeriv(const Sub<A,B> &expr,float dresult)
{
  addDeriv(expr.a, dresult);
  addDeriv(expr.b,-dresult);
}
}


template <typename AExpr,typename BExpr>
static Sub<AExpr,BExpr> operator-(const AExpr &a,const BExpr &b)
{
  return {a,b};
}


template <typename A,typename B>
static ScalarExpr<Sub<A,B>>
  operator-(const ScalarExpr<A> &a,const ScalarExpr<B> &b)
{
  return {{a.expr,b.expr}};
}


namespace {
template <typename A,typename B>
struct Mul {
  A a;
  B b;
};
}


namespace {
template <typename A,typename B>
static float evaluate(const Mul<A,B> &expr)
{
  return evaluate(expr.a)*evaluate(expr.b);
}
}


namespace {
template <typename A,typename B>
static void addDeriv(const Mul<A,B> &expr,float dresult)
{
  addDeriv(expr.a, dresult*evaluate(expr.b));
  addDeriv(expr.b, dresult*evaluate(expr.a));
}
}


template <typename AExpr,typename BExpr>
static Mul<AExpr,BExpr> operator*(const AExpr &a,const BExpr &b)
{
  return {a,b};
}


template <typename A,typename B>
static ScalarExpr<Mul<A,B>>
  operator*(const ScalarExpr<A> &a,const ScalarExpr<B> &b)
{
  return {{a.expr,b.expr}};
}


using FloatMat33 = Mat33<float>;


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


using DualMat33 = Mat33<DualFloat>;


static DualFloat dual(float av,float &da)
{
  return {av,da};
}


namespace {
template <typename T>
struct Vector {
  T x, y, z;

  friend ostream& operator<<(ostream &stream,const Vector &self)
  {
    stream << "[" << self.x << "," << self.y << "," << self.z << "]";
    return stream;
  }
};
}


using DualVector = Vector<DualFloat>;

using FloatVector = Vector<float>;


namespace {
template <typename X,typename Y,typename Z>
static void addDeriv(VectorFunc<X,Y,Z> f,const FloatVector &deriv)
{
  addDeriv(f.x,deriv.x);
  addDeriv(f.y,deriv.y);
  addDeriv(f.z,deriv.z);
}
}


namespace {
template <typename X,typename Y,typename Z>
static FloatVector evaluate(VectorFunc<X,Y,Z> f)
{
  float x = evaluate(f.x);
  float y = evaluate(f.y);
  float z = evaluate(f.z);
  return {x,y,z};
}
}


static DualVector dual(FloatVector v,FloatVector &dv)
{
  return {{v.x,dv.x},{v.y,dv.y},{v.z,dv.z}};
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
static void addDeriv(const ScalarDiv<A,B> &expr,float dresult)
{
  float aval = evaluate(expr.a);
  float bval = evaluate(expr.b);
  float da = 0;
  float db = 0;

  ddiv(dual(aval,da),dual(bval,db),dresult);
  addDeriv(expr.a, da);
  addDeriv(expr.b, db);
}
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


static float differenceBetween(const FloatVector &a,const FloatVector &b)
{
  float dx = differenceBetween(a.x, b.x);
  float dy = differenceBetween(a.y, b.y);
  float dz = differenceBetween(a.z, b.z);

  return max(dx,dy,dz);
}


template <typename Function>
static FloatVector finiteDeriv(Function f,FloatVector &v)
{
  float dx = finiteDeriv(f,v.x);
  float dy = finiteDeriv(f,v.y);
  float dz = finiteDeriv(f,v.z);
  return FloatVector{dx,dy,dz};
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


static FloatVector randomVector(RandomEngine &random_engine)
{
  float x = randomFloat(-1,1,random_engine);
  float y = randomFloat(-1,1,random_engine);
  float z = randomFloat(-1,1,random_engine);

  return FloatVector{x,y,z};
}


static FloatMat33 randomMat33(RandomEngine &random_engine)
{
  FloatVector x = randomVector(random_engine);
  FloatVector y = randomVector(random_engine);
  FloatVector z = randomVector(random_engine);

  float values[3][3] = {
    {x.x, x.y, x.z},
    {y.x, y.y, y.z},
    {z.x, z.y, z.z},
  };

  return FloatMat33(values);
}


namespace {
template <typename Expr>
struct VectorExpr {
  Expr expr;
};
}


#if 1
template <typename A,typename B>
static auto dot(const A a,const B b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
#else
namespace {
template <typename A,typename B>
struct Dot {
  A a;
  B b;
};
}


template <typename A,typename B>
static ScalarExpr<Dot<A,B>> dot(VectorExpr<A> a,VectorExpr<B> b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

static float dot(const FloatVector &a,const FloatVector &b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
#endif


namespace {
template <typename Expr>
struct ScalarExprVar {
  Expr expr;
};
}


namespace {
template <typename Expr> struct Evaluator;
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
struct Evaluator<Mat33<DualFloat>> {
  Mat33<DualFloat> m_expr;
  FloatMat33 m_value;

  Evaluator(const Mat33<DualFloat> &arg)
  : m_expr(arg), m_value(evaluate(arg))
  {
  }

  const FloatMat33& value() const { return m_value; }

  void addDeriv(const FloatMat33 &deriv) const
  {
    ::addDeriv(m_expr,deriv);
  }
};
}



namespace {
template <typename A,typename B>
struct Evaluator<Add<A,B>> {
  Add<A,B> expr;

  Evaluator(Add<A,B> expr)
  : expr(expr)
  {
  }

  float value() const
  {
    return evaluate(expr);
  }

  void addDeriv(float deriv) const
  {
    return ::addDeriv(expr,deriv);
  }
};
}


namespace {
template <typename A,typename B>
struct Evaluator<Mul<A,B>> {
  Mul<A,B> expr;
  float a;
  float b;

  Evaluator(Mul<A,B> expr)
  : expr(expr),
    a(evaluate(expr.a)),
    b(evaluate(expr.b))
  {
  }

  float value() const
  {
    return a*b;
  }

  void addDeriv(float deriv)
  {
    ::addDeriv(expr.b,deriv*a);
    ::addDeriv(expr.a,deriv*b);
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
static auto cofactor(const Mat &arg,int i,int j)
{
  auto a11 = arg[(i+1)%3][(j+1)%3];
  auto a12 = arg[(i+1)%3][(j+2)%3];
  auto a21 = arg[(i+2)%3][(j+1)%3];
  auto a22 = arg[(i+2)%3][(j+2)%3];

  return a11*a22 - a12*a21;
}


static void dcofactor(const DualMat33 &arg,int i,int j,float dresult)
{
  addDeriv(cofactor(arg,i,j),dresult);
}


namespace {
template <typename M>
struct Evaluator<Cofactor<M>> {
  Evaluator<M> m;
  int i;
  int j;

  Evaluator(const Cofactor<M> &c)
  : m{c.m},
    i(c.i),
    j(c.j)
  {
  }

  float value() const
  {
    return cofactor(m.value(),i,j);
  }

  void addDeriv(float dresult)
  {
    FloatMat33 dm = mat33All(0);
    DualMat33 dual_m = dual(m.value(),dm);
    dcofactor(dual_m,i,j,dresult);
    m.addDeriv(dm);
  }
};
}


namespace {
template <typename A,typename B>
struct Evaluator<Mat33Div<A,B>> {
  Evaluator<A> a_eval;
  Evaluator<B> b_eval;

  Evaluator(const Mat33Div<A,B> &expr)
  : a_eval(expr.a),
    b_eval(expr.b)
  {
  }

  FloatMat33 value() const
  {
    return a_eval.value() / b_eval.value();
  }

  void addDeriv(const FloatMat33 &deriv)
  {
    FloatMat33 da = mat33All(0);
    DualMat33 a = dual(a_eval.value(),da);
    float db = 0;
    DualFloat b = dual(b_eval.value(),db);

    for (int i=0; i!=3; ++i) {
      for (int j=0; j!=3; ++j) {
        ::addDeriv(expr(a[i][j])/expr(b),deriv[i][j]);
      }
    }

    a_eval.addDeriv(da);
    b_eval.addDeriv(db);
  }
};
}

namespace {
template <typename Expr>
struct VectorExprVar {
  Expr expr;

  FloatVector value = evaluate(expr);
  FloatVector dvalue{0,0,0};

  DualFloat x{value.x,dvalue.x};
  DualFloat y{value.y,dvalue.y};
  DualFloat z{value.z,dvalue.z};

  VectorExprVar(Expr expr)
  : expr(expr)
  {
  }

  ~VectorExprVar()
  {
    addDeriv(expr,dvalue);
  }
};
}


namespace {
template <typename Expr>
static float evaluate(const ScalarExprVar<Expr> &e)
{
  return evaluate(e.expr);
}
}


namespace {
template <typename Expr>
static void addDeriv(const ScalarExprVar<Expr> &e,float deriv)
{
  addDeriv(e.expr,deriv);
}
}


template <typename Expr>
static ScalarExprVar<Expr> var(const ScalarExpr<Expr> &expr)
{
  return {expr.expr};
}


template <typename Expr>
static VectorExprVar<Expr> var(const VectorExpr<Expr> &expr)
{
  return {expr.expr};
}


template <typename Expr>
static ScalarExprVar<Expr> var(const Expr &expr)
{
  return {expr};
}


namespace {
template <>
struct VectorExpr<DualVector> {
  ScalarExpr<DualFloat> x, y, z;
};
}


static VectorExpr<DualVector> expr(const DualVector &expr)
{
  return {{expr.x},{expr.y},{expr.z}};
}


static void ddot(DualVector a,DualVector b,float dresult)
{
  addDeriv(dot(expr(a),expr(b)),dresult);
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


template <typename M>
static void addDeriv(Transpose<M> transpose,const FloatMat33 &deriv)
{
  dtranspose(transpose.m,deriv);
}


template <typename M>
static Mat33Expr<Transpose<M>> transpose(const Mat33Expr<M> &m)
{
  return {{m.expr}};
}

namespace {
template <typename M>
struct Evaluator<Transpose<M>> {
  Evaluator<M> m_eval;

  Evaluator(const Transpose<M> &expr)
  : m_eval(expr.m)
  {
  }

  FloatMat33 value() const
  {
    return transpose(m_eval.value());
  }

  void addDeriv(const FloatMat33 &dresult)
  {
    const FloatMat33 &m = m_eval.value();
    FloatMat33 dm = mat33All(0);
    dtranspose(dual(m,dm),dresult);
    m_eval.addDeriv(dm);
  }
};
}


template <typename M>
static ScalarExpr<Cofactor<M>>
  cofactor(const Mat33Expr<M> &arg,int i,int j)
{
  return {{arg.expr,i,j}};
}


template <typename Mat>
static auto cofactorMatrix(const Mat &arg)
{
  using V = decltype(cofactor(arg,0,0));

  V values[3][3] = {
    {cofactor(arg,0,0),cofactor(arg,0,1),cofactor(arg,0,2)},
    {cofactor(arg,1,0),cofactor(arg,1,1),cofactor(arg,1,2)},
    {cofactor(arg,2,0),cofactor(arg,2,1),cofactor(arg,2,2)},
  };

  return Mat33<V>(values);
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


static void dcofactorMatrix(const DualMat33 &arg,const FloatMat33 &dresult)
{
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      dcofactor(arg,i,j,dresult[i][j]);
    }
  }
}


namespace {
template <typename M>
struct Evaluator<CofactorMatrixFunc<M>> {
  Evaluator<M> m_eval;
  FloatMat33 m;

  Evaluator(const CofactorMatrixFunc<M> &expr)
  : m_eval(expr.m),
    m(m_eval.value())
  {
  }

  FloatMat33 value() const
  {
    return cofactorMatrix(m);
  }

  void addDeriv(const FloatMat33 &dresult)
  {
    FloatMat33 dm = mat33All(0);

    dcofactorMatrix(dual(m,dm),dresult);
    m_eval.addDeriv(dm);
  }
};
}


template <typename Mat>
static auto determinant(const Mat &a)
{
  return
    a[0][0] * cofactor(a,0,0) +
    a[0][1] * cofactor(a,0,1) +
    a[0][2] * cofactor(a,0,2);
}


// We need a determinant function that works on a Mat33Expr

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


static void ddeterminant(const DualMat33 &a,float dresult)
{
  addDeriv(determinant(a),dresult);
}


namespace {
template <typename M>
struct Evaluator<Determinant<M>> {
  Evaluator<M> m_eval;

  Evaluator(const Determinant<M> &expr)
  : m_eval(expr.m)
  {
  }

  float value() const
  {
    return determinant(m_eval.value());
  }

  void addDeriv(float dresult)
  {
    const FloatMat33 &m = m_eval.value();
    FloatMat33 dm = mat33All(0);

    ddeterminant(dual(m,dm),dresult);

    m_eval.addDeriv(dm);
  }
};
}


template <typename M>
static void addDeriv(Determinant<M> expr,float deriv)
{
  ddeterminant(expr.m,deriv);
}


static Mat33Expr<DualMat33> expr(const DualMat33 &expr)
{
  return {expr};
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
  Evaluator<M> m_eval;
  FloatMat33 m = m_eval.value();
  FloatMat33 dm = mat33All(0);
  DualMat33 m_dual = dual(m,dm);
  Evaluator<
    decltype(genMat33Inv(expr(m_dual)).expr)
  > result_eval{genMat33Inv(expr(m_dual)).expr};

  Evaluator(const Inv<M> &expr_arg)
  : m_eval(expr_arg.m)
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

  ~Evaluator()
  {
    m_eval.addDeriv(dm);
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
static VectorExpr<VectorFunc<A,B,C>>
  vector(
    ScalarExpr<A> a,
    ScalarExpr<B> b,
    ScalarExpr<C> c
  )
{
  auto f = VectorFunc<A,B,C>{a.expr,b.expr,c.expr};
  return VectorExpr<VectorFunc<A,B,C>>{f};
}


template <typename Expr>
static float
  evalAndAddDeriv(const ScalarExpr<Expr> &scalar_expr,float dresult)
{
  auto e = scalar_expr.expr;
  Evaluator<decltype(e)> eval(e);
  float result = eval.value();
  eval.addDeriv(dresult);
  return result;
}


template <typename Expr>
static FloatMat33
  evalAndAddDeriv(const Mat33Expr<Expr> &mat33_expr,const FloatMat33& dresult)
{
  auto e = mat33_expr.expr;
  Evaluator<decltype(e)> eval(e);
  FloatMat33 result = eval.value();
  eval.addDeriv(dresult);
  return result;
}


namespace tests {
static void testScalarAddEvaluator()
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


#if 1
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
#endif


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


static void testAddOperatorDeriv()
{
  float av = 5;
  float bv = 6;
  float da = 0;
  float db = 0;
  auto a = dual(av,da);
  auto b = dual(bv,db);
  addDeriv(a+b,1);
  assertNear(da,1.0f,0);
  assertNear(db,1.0f,0);
}


static void testMultiplyDeriv()
{
  float a = 5;
  float b = 6;
  float da = 0;
  float db = 0;
  addDeriv(expr(dual(a,da))*expr(dual(b,db)),/*dresult*/1);
  auto f = [&](){ return a*b; };
  float fda = finiteDeriv(f,a);
  float fdb = finiteDeriv(f,b);
  assertNear(da,fda,1e-3);
  assertNear(db,fdb,1e-3);
}


static void testDot()
{
  RandomEngine random_engine{/*seed*/1};
  FloatVector v1 = randomVector(random_engine);
  FloatVector v2 = randomVector(random_engine);
  FloatVector dv1{0,0,0};
  FloatVector dv2{0,0,0};
  ddot(dual(v1,dv1),dual(v2,dv2),1);
  auto f = [&](){ return dot(v1,v2); };
  FloatVector fdv1 = finiteDeriv(f,v1);
  FloatVector fdv2 = finiteDeriv(f,v2);
  assertNear(dv1,fdv1,1e-4);
  assertNear(dv2,fdv2,1e-4);
}


static void testDot2()
{
  float dv = 0;
  {
    DualFloat v{3,dv};
    auto a = var(vector(expr(2),expr(v),expr(4)));
    auto b = var(vector(expr(5),expr(6),expr(7)));
    auto e = var(dot(a,b));
    addDeriv(e,1);
  }
  assertNear(dv,6.0f,1e-4);
}


#if 0
static void testAddingVectors()
{
  float dv = 0;

  {
    DualFloat v{3,dv};
    auto a = vector(1,2,v);
    auto b = vector(4,5,6);
    auto e = var(a+b);
    addDeriv(e,1);
  }

  assertNear(dv,1.0f,1e-4);
}
#endif


static void testMatTranspose()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      auto f = [&](){ return transpose(a)[i][j]; };
      FloatMat33 da = mat33All(0);
      FloatMat33 e = mat33All(0);
      e[i][j] = 1;
      dtranspose(dual(a,da),e);
      FloatMat33 fda = finiteDeriv(f,a);
      assertNear(da,fda,1e-4);
    }
  }
}


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


static void testCofactorDeriv()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      float dresult = 1;
      FloatMat33 da = mat33All(0);
      dcofactor(dual(a,da),i,j,dresult);
      auto f = [&](){ return cofactor(a,i,j); };
      FloatMat33 fda = finiteDeriv(f,a);
      assertNear(da,fda,1e-4);
    }
  }
}


static void testMat33CofactorMatrixDeriv()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      FloatMat33 dresult = mat33All(0);
      dresult[i][j] = 1;
      FloatMat33 da = mat33All(0);
      dcofactorMatrix(dual(a,da),dresult);
      auto f = [&](){ return cofactorMatrix(a)[i][j]; };
      float fdaij = finiteDeriv(f,a[i][j]);
      assertNear(da[i][j],fdaij,1e-4);
    }
  }
}


static void testDDeterminant()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);
  FloatMat33 da = mat33All(0);

  ddeterminant(dual(a,da),/*dresult*/1);
  auto f = [&](){ return determinant(a); };
  FloatMat33 fda = finiteDeriv(f,a);
  assertNear(da,fda,1e-4);
}


static void testDeterminantDeriv()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);
  FloatMat33 da = mat33All(0);

  addDeriv(determinant(expr(dual(a,da))).expr,/*dresult*/1);
  auto f = [&](){ return determinant(a); };
  FloatMat33 fda = finiteDeriv(f,a);
  assertNear(da,fda,1e-4);
}
}


int main()
{
  tests::testScalarAddEvaluator();
  tests::testScalarMulEvaluator();
  tests::testCofactorEvaluator();
  tests::testCofactorMatrixEvaluator();
  tests::testDeterminantEvaluator();
  tests::testTransposeEvaluator();
  tests::testMat33DivEvaluator();
  tests::testMat33InvEvaluator();
  tests::testAddOperatorDeriv();
  tests::testMultiplyDeriv();
  tests::testDot();
  tests::testDot2();
  // tests::testAddingVectors();
  tests::testMatTranspose();
  tests::testMat33Inv();
  tests::testCofactorDeriv();
  tests::testMat33CofactorMatrixDeriv();
  tests::testDDeterminant();
  tests::testDeterminantDeriv();
}

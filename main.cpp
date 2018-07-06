#include <iostream>
#include <cmath>
#include <cfloat>
#include "random.hpp"

using std::cerr;
using std::ostream;
using std::max;


#define assertNear(actual,expected,tolerance) \
  (assertNearHelper(actual,expected,tolerance,__FILE__,__LINE__))



static void addDeriv(const float &,float)
{
  // Nothing to do if we are adding a derivative to a constant.
}


static float evaluate(float arg)
{
  return arg;
}


namespace {
struct DualFloat {
  float value;
  float &deriv;

  friend float evaluate(const DualFloat &dual)
  {
    return dual.value;
  }

  friend void addDeriv(const DualFloat &dual,float dresult)
  {
    dual.deriv += dresult;
  }
};
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


template <typename T>
static void addDeriv(const Mat33<T> &m,const Mat33<float> &dm)
{
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      addDeriv(m[i][j],dm[i][j]);
    }
  }
}


namespace {
template <typename AExpr,typename BExpr>
struct Add {
  AExpr a;
  BExpr b;

  friend auto evaluate(const Add &expr)
  {
    return evaluate(expr.a) + evaluate(expr.b);
  }

  friend void addDeriv(const Add &expr,float dresult)
  {
    addDeriv(expr.a,dresult);
    addDeriv(expr.b,dresult);
  }
};
}


template <typename AExpr,typename BExpr>
static Add<AExpr,BExpr> operator+(const AExpr &a,const BExpr &b)
{
  return {a,b};
}


namespace {
template <typename AExpr,typename BExpr>
struct Sub {
  AExpr a;
  BExpr b;

  friend auto evaluate(const Sub<AExpr,BExpr> &expr)
  {
    return evaluate(expr.a) - evaluate(expr.b);
  }

  friend void addDeriv(const Sub &expr,float dresult)
  {
    addDeriv(expr.a, dresult);
    addDeriv(expr.b,-dresult);
  }
};
}


template <typename AExpr,typename BExpr>
static Sub<AExpr,BExpr> operator-(const AExpr &a,const BExpr &b)
{
  return {a,b};
}


namespace {
template <typename AExpr,typename BExpr>
struct Mul {
  AExpr a;
  BExpr b;

  friend auto evaluate(const Mul &expr)
  {
    return evaluate(expr.a)*evaluate(expr.b);
  }

  friend void addDeriv(const Mul &expr,float dresult)
  {
    addDeriv(expr.a, dresult*evaluate(expr.b));
    addDeriv(expr.b, dresult*evaluate(expr.a));
  }
};
}


template <typename AExpr,typename BExpr>
static Mul<AExpr,BExpr> operator*(const AExpr &a,const BExpr &b)
{
  return {a,b};
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


static void ddiv(DualMat33 a,DualFloat b,const FloatMat33 &dresult);


namespace {
template <typename AExpr,typename BExpr>
struct Div {
  AExpr a;
  BExpr b;
  using Result = decltype(evaluate(a)/evaluate(b));
  using A = decltype(evaluate(a));
  using B = decltype(evaluate(b));

  friend auto evaluate(const Div &expr)
  {
    return evaluate(expr.a) * evaluate(expr.b);
  }

  friend void addDeriv(const Div &expr,Result dresult)
  {
    A da = zero<A>();
    B db = zero<B>();
    auto aval = evaluate(expr.a);
    auto bval = evaluate(expr.b);

    ddiv(dual(aval,da),dual(bval,db),dresult);

    addDeriv(expr.a,da);
    addDeriv(expr.b,db);
  }
};
}


template <typename AExpr,typename BExpr>
static Div<AExpr,BExpr> operator/(AExpr a,BExpr b)
{
  return {a,b};
}


static void
  ddiv(
    DualMat33 a,
    DualFloat b,
    const FloatMat33 &dresult
  )
{
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      addDeriv(a[i][j]/b,dresult[i][j]);
    }
  }
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


template <typename A,typename B>
static auto dot(const A a,const B b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}


static void ddot(DualVector a,DualVector b,float dresult)
{
  addDeriv(dot(a,b),dresult);
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


static void dtranspose(const DualMat33 &a,const FloatMat33 &dresult)
{
  addDeriv(a,transpose(dresult));
}


template <typename T>
static auto cofactor(const Mat33<T> &arg,int i,int j)
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


template <typename T>
static auto cofactorMatrix(const Mat33<T> &arg)
{
  using V = decltype(cofactor(arg,0,0));

  V values[3][3] = {
    {cofactor(arg,0,0),cofactor(arg,0,1),cofactor(arg,0,2)},
    {cofactor(arg,1,0),cofactor(arg,1,1),cofactor(arg,1,2)},
    {cofactor(arg,2,0),cofactor(arg,2,1),cofactor(arg,2,2)},
  };

  return Mat33<V>(values);
}


static void dcofactorMatrix(const DualMat33 &arg,const FloatMat33 &dresult)
{
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      dcofactor(arg,i,j,dresult[i][j]);
    }
  }
}


template <typename T>
static auto determinant(const Mat33<T> &a)
{
  return
    a[0][0] * cofactor(a,0,0) +
    a[0][1] * cofactor(a,0,1) +
    a[0][2] * cofactor(a,0,2);
}


static void ddeterminant(const DualMat33 &a,float dresult)
{
  addDeriv(determinant(a),dresult);
}


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
static auto mat33Inv(const Mat33<T> &arg)
{
  return transpose(cofactorMatrix(arg))/determinant(arg);
}


static void
  dmat33Inv(const DualMat33 &arg,const FloatMat33 &dresult)
{
  addDeriv(mat33Inv(arg),dresult);
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


namespace {
template <typename X,typename Y,typename Z>
struct VectorExpr {
  X x;
  Y y;
  Z z;
};
}


template <typename A,typename B,typename C>
static VectorExpr<A,B,C> vector(const A& a,const B& b,const C& c)
{
  return {a,b,c};
}


namespace tests {
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
  addDeriv(dual(a,da)*dual(b,db),/*dresult*/1);
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


namespace {
template <typename Expr>
struct ExprVar {
  Expr expr;
};
}


template <typename Expr>
static ExprVar<Expr> var(const Expr &expr)
{
  return {expr};
}


template <typename Expr>
static void addDeriv(ExprVar<Expr> &e,float deriv)
{
  addDeriv(e.expr,deriv);
}


static void testDot2()
{
  float dv = 0;
  DualFloat v{3,dv};
  auto a = vector(2,v,4);
  auto b = vector(5,6,7);
  auto e = var(dot(a,b));
  addDeriv(e,1);
  assertNear(dv,6.0f,1e-4);
}


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


static void testMat33DivDeriv()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);
  float b = randomFloat(-1,1,random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      FloatMat33 dresult = mat33All(0);
      dresult[i][j] = 1;
      FloatMat33 da = mat33All(0);
      float db = 0;
      addDeriv(dual(a,da)/dual(b,db),dresult);
      auto f = [&](){ return (a/b)[i][j]; };
      float fdaij = finiteDeriv(f,a[i][j]);
      assertNear(fdaij,da[i][j],1e-4);
      float fdb = finiteDeriv(f,b);
      assertNear(fdb,db,1e-4);
    }
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


static void testDeterminantDeriv()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);
  FloatMat33 da = mat33All(0);

  ddeterminant(dual(a,da),/*dresult*/1);
  auto f = [&](){ return determinant(a); };
  FloatMat33 fda = finiteDeriv(f,a);
  assertNear(da,fda,1e-4);
}


static void testMatInvDeriv()
{
  RandomEngine random_engine(/*seed*/1);
  FloatMat33 a = randomMat33(random_engine);

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      auto f = [&](){ return mat33Inv(a)[i][j]; };
      FloatMat33 dresult = mat33All(0);
      dresult[i][j] = 1;
      FloatMat33 da = mat33All(0);
      dmat33Inv(dual(a,da),dresult);
      FloatMat33 fda = finiteDeriv(f,a);
      assertNear(da,fda,2e-4);
    }
  }
}
}


int main()
{
  tests::testAddOperatorDeriv();
  tests::testMultiplyDeriv();
  tests::testDot();
  tests::testDot2();
  tests::testMatTranspose();
  tests::testMat33Inv();
  tests::testMat33DivDeriv();
  tests::testCofactorDeriv();
  tests::testMat33CofactorMatrixDeriv();
  tests::testDeterminantDeriv();
  tests::testMatInvDeriv();
}

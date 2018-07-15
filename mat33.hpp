#ifndef MAT33_HPP
#define MAT33_HPP

#include <iostream>
#include <cassert>
#include <cmath>
#include "vec3.hpp"


template <typename Matrix> decltype(auto) element(Matrix &matrix,int i,int j)
{
  return matrix.element(matrix,i,j);
}



template <typename Matrix>
struct RowRef {
  Matrix &matrix;
  int i;

  decltype(auto) operator[](int j)
  {
    return element(matrix,i,j);
  }
};


template <typename Self>
auto row(Self &arg,int i)
{
  return RowRef<Self>{arg,i};
}


template <typename Matrix>
struct ColRef {
  Matrix &matrix;
  int j;

  decltype(auto) operator[](int i)
  {
    return element(matrix,i,j);
  }
};


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

  friend std::ostream& operator<<(std::ostream &stream,const Mat33 &arg)
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

  bool operator==(const Mat33 &arg) const
  {
    for (int i=0; i!=3; ++i) {
      for (int j=0; j!=3; ++j) {
        if (values[i][j] != arg.values[i][j]) {
          return false;
        }
      }
    }

    return true;
  }

  T values[3][3];
};


template <typename T>
auto col(Mat33<T> &arg,int j)
{
  return ColRef<Mat33<T>>{arg,j};
}


template <typename T>
auto col(const Mat33<T> &arg,int j)
{
  return ColRef<const Mat33<T>>{arg,j};
}


inline Mat33<float> operator*(const Mat33<float> &a,const Mat33<float> &b)
{
  using T = float;
  T values[3][3];

  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      T sum = 0;

      for (int k=0; k!=3; ++k) {
        sum += a[i][k]*b[k][j];
      }

      values[i][j] = sum;
    }
  }

  return Mat33<T>(values);
}


using FloatMat33 = Mat33<float>;


template <typename T>
Mat33<T> transpose(const Mat33<T> &a)
{
  T values[3][3] = {
    {a[0][0],a[1][0],a[2][0]},
    {a[0][1],a[1][1],a[2][1]},
    {a[0][2],a[1][2],a[2][2]},
  };

  return Mat33<T>(values);
}


template <typename T>
Mat33<T> mat33(const T (&values)[3][3])
{
  return Mat33<T>(values);
}


inline FloatMat33 mat33All(float arg)
{
  float values[3][3] = {
    {arg,arg,arg},
    {arg,arg,arg},
    {arg,arg,arg},
  };

  return FloatMat33(values);
}


inline FloatMat33 zeroMat33()
{
  return mat33All(0);
}


template <typename A,typename B>
auto genMat33Div(A&& a,B&& b)
{
  decltype(a[0][0]/b) values[3][3] = {
    { a[0][0]/b, a[0][1]/b, a[0][2]/b },
    { a[1][0]/b, a[1][1]/b, a[1][2]/b },
    { a[2][0]/b, a[2][1]/b, a[2][2]/b }
  };

  return mat33(values);
}


template <typename Mat>
auto genCofactor(const Mat &arg,int i,int j)
{
  auto a11 = arg[(i+1)%3][(j+1)%3];
  auto a12 = arg[(i+1)%3][(j+2)%3];
  auto a21 = arg[(i+2)%3][(j+1)%3];
  auto a22 = arg[(i+2)%3][(j+2)%3];

  return a11*a22 - a12*a21;
}


inline float cofactor(const FloatMat33 &arg,int i,int j)
{
  return genCofactor(arg,i,j);
}


template <typename Mat>
auto genCofactorMatrix(const Mat &arg)
{
  using V = decltype(genCofactor(arg,0,0));

  V values[3][3] = {
    {genCofactor(arg,0,0),genCofactor(arg,0,1),genCofactor(arg,0,2)},
    {genCofactor(arg,1,0),genCofactor(arg,1,1),genCofactor(arg,1,2)},
    {genCofactor(arg,2,0),genCofactor(arg,2,1),genCofactor(arg,2,2)},
  };

  return mat33(values);
}


inline auto cofactorMatrix(const FloatMat33 &arg)
{
  return genCofactorMatrix(arg);
}


template <typename Mat>
auto genDeterminant(const Mat &a)
{
  return
    a[0][0] * genCofactor(a,0,0) +
    a[0][1] * genCofactor(a,0,1) +
    a[0][2] * genCofactor(a,0,2);
}


inline auto determinant(const FloatMat33 &a)
{
  return genDeterminant(a);
}


inline FloatMat33 operator/(const FloatMat33 &a,float b)
{
  return genMat33Div(a,b);
}


template <typename T>
auto genMat33Inv(const T &arg)
{
  return transpose(cofactorMatrix(arg))/determinant(arg);
}


template <typename T>
auto mat33Inv(const Mat33<T> &arg)
{
  return genMat33Inv(arg);
}


inline FloatMat33 mat33Identity()
{
  float values[3][3] = {
    {1,0,0},
    {0,1,0},
    {0,0,1},
  };

  return FloatMat33(values);
}


inline FloatMat33 rotX(float sina,float cosa)
{
  float values[3][3] = {
    {1,     0,    0},
    {0,  cosa, sina},
    {0, -sina, cosa}
  };

  return mat33(values);
}

inline FloatMat33 rotX(float angle)
{
  return rotX(sinf(angle),cosf(angle));
}


template <typename T>
inline Vec3<T> vec3(ColRef<const Mat33<T>> c)
{
  return {c[0],c[1],c[2]};
}


template <typename T>
inline Vec3<T> vec3(ColRef<Mat33<T>> c)
{
  return {c[0],c[1],c[2]};
}


template <typename T>
inline Mat33<T> columns(const Vec3<T> &c1,const Vec3<T> &c2,const Vec3<T> &c3)
{
  T values[3][3] = {
    {c1.x(),c2.x(),c3.x()},
    {c1.y(),c2.y(),c3.y()},
    {c1.z(),c2.z(),c3.z()},
  };

  return {values};
}


#endif /* MAT33_HPP */

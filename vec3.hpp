#ifndef VEC3_HPP
#define VEC3_HPP


#include <iostream>
#include "differencebetween.hpp"
#include "max3.hpp"


namespace autodiff {


template <typename T>
class Vec3 {
  public:
    Vec3(T x,T y,T z) : _x(x), _y(y), _z(z) { }

    T x() const { return _x; }
    T y() const { return _y; }
    T z() const { return _z; }

    T& x() { return _x; }
    T& y() { return _y; }
    T& z() { return _z; }

    friend std::ostream& operator<<(std::ostream &stream,const Vec3 &self)
    {
      stream << "[" << self.x() << "," << self.y() << "," << self.z() << "]";
      return stream;
    }

    bool operator==(const Vec3 &arg) const
    {
      return _x==arg._x && _y==arg._y && _z==arg._z;
    }

  private:
    T _x, _y, _z;
};


template <typename A,typename B>
auto genDot(const A& a,const B& b)
{
  return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
}


template <typename T>
T dot(const Vec3<T> &a,const Vec3<T> &b)
{
  return genDot(a,b);
}


inline Vec3<float> operator+(const Vec3<float> &a,const Vec3<float> &b)
{
  return {a.x() + b.x(), a.y() + b.y(), a.z() + b.z()};
}


inline Vec3<float> vec3(float x,float y,float z)
{
  return {x,y,z};
}


inline Vec3<float> vec3(int x,int y,int z)
{
  return {float(x),float(y),float(z)};
}


inline Vec3<float> operator-(const Vec3<float> &a,const Vec3<float> &b)
{
  return {a.x() - b.x(), a.y() - b.y(), a.z() - b.z()};
}


template <typename T>
Vec3<T> operator-(const Vec3<T> &v)
{
  return {-v.x(),-v.y(),-v.z()};
}


inline Vec3<float> operator/(const Vec3<float> &a,float b)
{
  return {a.x()/b, a.y()/b, a.z()/b};
}


inline Vec3<float> operator*(const Vec3<float> &a,float b)
{
  return {a.x()*b, a.y()*b, a.z()*b};
}


inline Vec3<float> operator*(float b,const Vec3<float> &a)
{
  return {a.x()*b, a.y()*b, a.z()*b};
}


template <typename T>
T mag(const Vec3<T> &v)
{
  return sqrt(dot(v,v));
}


using FloatVec3 = Vec3<float>;


inline float differenceBetween(const FloatVec3 &a,const FloatVec3 &b)
{
  float dx = ::differenceBetween(a.x(), b.x());
  float dy = ::differenceBetween(a.y(), b.y());
  float dz = ::differenceBetween(a.z(), b.z());

  return max(dx,dy,dz);
}


}



#endif /* VEC3_HPP */

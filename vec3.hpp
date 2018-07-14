#ifndef VEC3_HPP
#define VEC3_HPP


#include <iostream>


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

  private:
    T _x, _y, _z;
};


template <typename A,typename B>
auto genDot(const A a,const B b)
{
  return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
}


template <typename T>
T dot(const Vec3<T> &a,const Vec3<T> &b)
{
  return genDot(a,b);
}


template <typename T>
Vec3<T> operator+(const Vec3<T> &a,const Vec3<T> &b)
{
  return {a.x() + b.x(), a.y() + b.y(), a.z() + b.z()};
}


template <typename T>
Vec3<T> operator-(const Vec3<T> &a,const Vec3<T> &b)
{
  return {a.x() - b.x(), a.y() - b.y(), a.z() - b.z()};
}



template <typename T>
Vec3<T> operator/(const Vec3<T> &a,T b)
{
  return {a.x()/b, a.y()/b, a.z()/b};
}


template <typename T>
Vec3<T> operator*(const Vec3<T> &a,T b)
{
  return {a.x()*b, a.y()*b, a.z()*b};
}


template <typename T>
T mag(const Vec3<T> &v)
{
  return sqrt(dot(v,v));
}


using FloatVec3 = Vec3<float>;



#endif /* VEC3_HPP */

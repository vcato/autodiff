#include "random.hpp"

template <typename T>
autodiff::Vec3<T> randomVec3(RandomEngine &random_engine)
{
  T x = random<T>(-1,1,random_engine);
  T y = random<T>(-1,1,random_engine);
  T z = random<T>(-1,1,random_engine);

  return autodiff::Vec3<T>{x,y,z};
}

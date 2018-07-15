#include "randomvec3.hpp"

template <typename T>
autodiff::Mat33<T> randomMat33(RandomEngine &random_engine)
{
  using autodiff::Vec3;

  Vec3<T> x = randomVec3<T>(random_engine);
  Vec3<T> y = randomVec3<T>(random_engine);
  Vec3<T> z = randomVec3<T>(random_engine);

  T values[3][3] = {
    {x.x(), x.y(), x.z()},
    {y.x(), y.y(), y.z()},
    {z.x(), z.y(), z.z()},
  };

  return autodiff::Mat33<T>(values);
}

#include "random.hpp"

#include <random>

struct RandomEngine::Impl {
  std::mt19937 std_engine;

  Impl(int seed)
  : std_engine(seed)
  {
  }
};


RandomEngine::RandomEngine(int seed)
: impl_ptr(new Impl(seed))
{
}


RandomEngine::~RandomEngine()
{
  delete impl_ptr;
}


template <typename T>
T random(T low,T high,RandomEngine &random_engine)
{
  std::mt19937 &std_engine = random_engine.impl().std_engine;

  return std::uniform_real_distribution<T>(low,high)(std_engine);
}


template float random<float>(float,float,RandomEngine &);

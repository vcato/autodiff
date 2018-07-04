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


float randomFloat(float low,float high,RandomEngine &random_engine)
{
  std::mt19937 &std_engine = random_engine.impl().std_engine;

  return std::uniform_real_distribution<float>(low,high)(std_engine);
}

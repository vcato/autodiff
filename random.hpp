#ifndef RANDOM_HPP
#define RANDOM_HPP


#include <cassert>


struct RandomEngine {
  RandomEngine(int seed);

  ~RandomEngine();

  struct Impl;

  Impl *impl_ptr;
  Impl &impl() { assert(impl_ptr); return *impl_ptr; }
};


template <typename T>
extern T random(T low,T high,RandomEngine &random_engine);


#endif /* RANDOM_HPP */

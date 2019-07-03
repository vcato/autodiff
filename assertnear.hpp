#include "differencebetween.hpp"


template <typename T>
void
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

  std::cerr << "file: " << file << "\n";
  std::cerr << "line: " << line << "\n";
  std::cerr << "actual: " << actual << "\n";
  std::cerr << "expected: " << expected << "\n";
  std::cerr << "delta: " << delta << "\n";
  std::cerr << "tolerance: " << tolerance << "\n";
  assert(false);
}


#define assertNear(actual,expected,tolerance) \
  (assertNearHelper(actual,expected,tolerance,__FILE__,__LINE__))

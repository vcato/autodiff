#include "../mat33expr.hpp"


#if 0
void
  testFunction(
    const FloatMat33 &a,FloatMat33 &da,
    const FloatMat33 &b,FloatMat33 &db,
    const FloatMat33 &c,FloatMat33 &dc,
    const FloatMat33 &dresult
  )
{
  evalAndAddDeriv(expr(a,da)*expr(b,db)*expr(c,dc),dresult);
}
#else
void
  testFunction(
    float a,float &da,
    float b,float &db,
    float c,float &dc,
    float dresult
  )
{
  evalAndAddDeriv(expr(a,da)*expr(b,db)*expr(c,dc),dresult);
}
#endif

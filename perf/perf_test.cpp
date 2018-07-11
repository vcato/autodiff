#include "../vec3expr.hpp"
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
#endif

#if 0
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


#if 1
void testFunction(float &dv)
{
  auto a = vec3(expr(2),expr(3,dv),expr(4));
  auto b = vec3(expr(5),expr(6),   expr(7));
  evalAndAddDeriv(dot(a,b),1);
}
#endif

template <typename Function>
inline float finiteDeriv(Function f,float &v)
{
  float h = 1e-3;
  float old_value = v;
  v = old_value - h;
  float value1 = f();
  v = old_value + h;
  float value2 = f();
  v = old_value;
  return (value2-value1)/(2*h);
}

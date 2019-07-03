namespace autodiff {
template <typename Function>
FloatVec3 finiteDeriv(Function f,FloatVec3 &v)
{
  float x = ::finiteDeriv(f,v.x());
  float y = ::finiteDeriv(f,v.y());
  float z = ::finiteDeriv(f,v.z());

  return {x,y,z};
}
}

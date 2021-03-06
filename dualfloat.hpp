#ifndef DUALFLOAT_HPP
#define DUALFLOAT_HPP


namespace autodiff {


struct DualFloat {
  float value;
  float &deriv;
};


inline DualFloat dual(float av,float &da)
{
  return {av,da};
}


}

#endif /* DUALFLOAT_HPP */

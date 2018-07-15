namespace autodiff {


template <typename ExprWrapper,typename Value>
Value evalAndAddDeriv(ExprWrapper &&wrapper,const Value& dresult)
{
  auto e = internal(wrapper);
  Evaluator<decltype(e)> eval(e);
  Value result = eval.value();
  eval.addDeriv(dresult);
  return result;
}


}

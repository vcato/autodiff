namespace autodiff {


template <typename ExprWrapper,typename Value>
auto evalAndAddDeriv(ExprWrapper &&wrapper,const Value& dresult)
{
  auto e = internal(wrapper);
  Evaluator<decltype(e)> eval(e);
  auto result = eval.value();
  eval.addDeriv(dresult);
  return result;
}


}

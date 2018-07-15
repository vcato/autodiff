template <typename ExprWrapper,typename Value>
Value evalAndAddDeriv(const ExprWrapper &mat33_expr,const Value& dresult)
{
  auto e = mat33_expr.expr;
  Evaluator<decltype(e)> eval(e);
  Value result = eval.value();
  eval.addDeriv(dresult);
  return result;
}

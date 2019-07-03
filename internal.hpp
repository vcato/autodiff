namespace autodiff {
template <typename T> auto internal(void *);
  // forward declaration of internal as a template, so that be don't get
  // a compile error.  There is no actual version of this function
  // which takes a void*.  Specific overrides of internal() are created
  // for various types.
}

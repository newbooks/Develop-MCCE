#ifndef DDM__INTERNAL__ANNOTATION_H__INCLUDED
#define DDM__INTERNAL__ANNOTATION_H__INCLUDED

namespace ddm {

static void prevent_opt_elimination() {
  asm volatile("");
}

/**
 * Prevent elimination of variable in compiler optimizations.
 */
template <typename T>
void prevent_opt_elimination(T && var) {
  asm volatile("" : "+r" (var)); // noop-read
}

} // namespace ddm

#endif // DDM__INTERNAL__ANNOTATION_H__INCLUDED

by Bob Cook
professorcook.org
Jan. 5 2017

As far as I know this is the only C99 implementation.
Therefore, certain deficiencies in the standard are addressed.

1. Objects have a capitalized name, values are lower case.
2. Handles that are defined constants are negative, GrB_NULL is a zero handle.
   Handles (>=1) index a GrB_Object array.
3. All overloaded function names have been replaced.
4. A new object GrB_Variable was added to handle overloaded arguments, unaryop, and binaryop.
5. Eliminated all GrB_Plus_INT32 etc. since GrB_Variable encodes types.
6. For methods with a typed constant argument, replaced with int, double, GrB_Variable versions to
   avoid a blizzard of names like opengl.
7. Added a context argument to unaryop and binaryop so that functions could be reentrant and not
   rely on global variables as in the 2nd example.
8. Replaced unaryop and binaryop callbacks with generic GrB_Variable definitions.

Expanding the code would be an interesting project for a good undergrad or graduate student.
Caveat1: some of the functions (e.g. new...) are not thread safe
# GraphBLAS

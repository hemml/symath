# symath

A simple (and rather fast) symbolic math package for Common Lisp. I'm using it in my CFD-code generation project (not published yet), where it deals with huge expressions with thousands of terms, sometimes reducing their sizes to dozens of terms.

This is a really simple package, implements `(symplify expr)` function. Here `expr` is an algebraic expression like `'(/ (+ a b) a)`, a simplified version will be returned.

Expression can contain numbers, symbols, subexpressions (functions) and arrays. 1D arrays will be treated as vectors, square 2D arrays as matrices. Matrices and vectors can be multiplied to numbers and other vectors/matrices. The method `(array-multiply x y)` is exported, overload it to implement tensor multiplication or something else.

In the expressions some special functions can be used:

- `(sqr x)` - multiple x to itself. Works with vectors (returns vector length) and matrices.
- `(vector ...)` - construct a vector from arguments.
- `(aref idx array)` - get array element.

These functions can be transformed during simplification, all other functions will be kept as is (but their arguments will be simplified). Trigonometry functions support is planned, but not yet implemented.

The simplification algorithm is rather complicated, it contains the following main steps:

1. All vectoir/matrix operations (if any) are implemented, producing final vector/marix, and each element is simplified separately
2. All functions with numeric arguments are computed, where possible
3. Term reduction is performed for divisions and subtraction
4. Bracketing common factors performed. Most common factors are extracted first

**WARNING:** The library is in beta stage, it may be buggy, please don't rely on it completely, always test results for validity!

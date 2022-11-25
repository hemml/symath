# symath

A simple symbolic math package for Common Lisp.

This is a really simple package, implements `(symplify expr)` function. Here `expr` is an algebraic expression like `'(/ (+ a b) a)`, a simplified version will be returned.

Expression can contain numbers, symbols, subexpressions (functions) and arrays. 1D arrays will be treated as vectors, square 2D arrays as matrices. Matrices and vectors can be multiplied to numbers and other vectors/matrices. The method `(array-multiply x y)` is exported, overload it to implement tensor multiplication or something else.

In the expressions some special functions can be use:

- `(sqr x)` - multiple x to itself. Works with vectors (returns vector length) and matrices.
- `(vector ...)` - construct a vector from arguments.
- `(aref idx array)` - get array element.

The simplification algorithm is rather complicated, it contains the following main steps:

1. All vectoir/matrix operations (if any) are implemented, producing final vector/marix, and each element is simplified separately.
2. All functions with numeric arguments are computed, where possible.
3. Term reduction is performed for divisions and subtraction.
4. Bracketing common factors performed. Most common factors are extracted first.

Currently, only simple algebraic expressions

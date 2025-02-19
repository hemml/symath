# symath

**WARNING:** The library is in beta stage, it may be buggy, please don't rely on it completely, always test results for validity!


A simple (and rather fast) symbolic math package for Common Lisp. I'm using it in my CFD-code generation project (not published yet), where it deals with huge expressions with thousands of terms, sometimes reducing their sizes to dozens of terms.

This is a really simple package, contains mostly the `(symplify expr &key no-cache)` function. Here `expr` is an algebraic expressions like `'(/ (+ a b) a)`, a simplified version will be returned. If `:no-cache t` is specified, the caching will be disabled (useful for debugging).

Expression can contain numbers, symbols, subexpressions (functions) and arrays. 1D arrays will be treated as vectors, square 2D arrays as matrices. Matrices and vectors can be multiplied to numbers and other vectors/matrices. The method `(array-multiply x y)` is exported, overload it to implement tensor multiplication or something else.

In the expressions some special functions can be used:

- `(sqr x)` - multiple x to itself. Works with vectors (returns vector length) and matrices.
- `(vector ...)` - construct a vector from arguments.
- `(aref idx array)` - get array element.
- `(exp x)` - an exponent
- `(log x &optional b)` - a logarithm
- `(diff expr var)` - a differentiation, see below

These functions can be transformed during simplification, all other functions will be kept as is (but their arguments will be simplified). Trigonometry functions support is planned, but not yet implemented.

## Differentiation

You can differentiate expressions with `(diff expr1 var)` function with `simplify`:
```
  (simplify '(diff (exp (cos x)) x))
  => (* -1 (EXP (COS X)) (SIN X))
```

The simplification algorithm knows how to differentiate the following functions: `+ - * / exp expt log sin cos tan ctan asyn acos atan actan`. Any other functions will be left under the `diff`:
```
  (simplify '(diff (exp (foo x)) x))
  (* (EXP (FOO X)) (DIFF (FOO X) X))
```
But you can add a differentiation rule for your function using the `with-templates` macro. Here you can see, how to use a finite difference (a central one with the step `delta`) for the `foo` function:
```
  (with-templates (((diff (foo $1) $2) =>
                    `(/ (- (foo ,(replace-subexpr $1 $2 `(+ ,$2 delta)))
                           (foo ,(replace-subexpr $1 $2 `(- ,$2 delta))))
                        (* 2 delta))))
    (simplify '(diff (exp (foo x)) x) :no-cache t))
  => (/ (* (- (* 1/2 (FOO (+ X DELTA))) (* 1/2 (FOO (- X DELTA)))) (EXP (FOO X)))
     DELTA)
```
The first argument of `with-templates` is a list of templates. Each template must have the following form:
```
(pattern => &rest code)
```
The symbol `=>` can be ommited. The pattern is an expression, where the special symbols denote:
- `$XX` - any expression, may be $1, $aaa, etc. All the `$`-symbols, except the single `$` can be referred in the code as variables.
- `@XX` - a number or a constant (see below)
- `_XX` - not a number/constant
- `&rest var` - keep the rest terms in the `var`, like: `(+ &rest args)`

The code must return an expression to replace the template or `NIL` to skip the replacement (you can perform any checks in the code to ensure what the expression is really satisfying your pattern). All templates are applied recursively while it is possible to find any pattern in the expression. While debugging your templates you can use `:no-cache t` while calling the `simplify`, to prevent caching wrong results.

You can specify some variables as a constants with `with-constants` macro:
```
  (simplify '(diff (expt x y) x))
  => (* (+ (* (DIFF Y X) (LOG X)) (/ Y X)) (EXPT X Y))
  (with-constants (y) (simplify '(diff (expt x y) x) :no-cache t))
  => (* Y (EXPT X (- Y 1)))
```

Feel free to contribute your diff-templates for specific functions, via issues or pull requests on the github repository.

## Extracting a subexpression(s)

`(extract-subexpr expr subexpr &key expand)` function can be used to isolate specific subexpression. It will convert the `expr` into the form `(subexpr^n)*e1+e2`, and return `(values n e1 e2)`. If it cannot isolate the `subexpr`, it will return `(values 0 0 expr)`. If the `expand` is set to `T`, the `expr` and `subexpr` will be transformed to have a better chance for extraction - all brackets will be opened in `e1` and `e2`.
This function is a very simple utility function and does not perform any transformations to solve the equation, the `subexpr` must be present in the `expr` more or less explicitly:

```
(extract-subexpr '(+ (sqrt x) 1) 'x) ;; will return (values 1/2 1 1)
(extract-subexpr '(+ (sqrt (+ x y)) 1) '(+ x y)) ;; will also return (values 1/2 1 1)
```

`(get-polynome-cfs expr subexpr &key expand)` returns a plist like `((0 . cf0) (1 . cf1) ...)` where `сf0`, `cf1`, etc is a polynomial coefficients of `expr` against `subexpr`. The `expand` argument works in the same way as in `extract-subexpr`. **NB**: The resulting alist will contain only nonzero coefficients!

`(replace-subexpr e e1 e2)` - find and replace all occurrences of expression `e1` in expression `e` and replace them to `e1`, returning the modified expression. The second returned value will be number of replacements. The expression `e1` must be explicitly present in `e` to be replaced, but the function can find and replace, for example, a subexpression `(+ c d)` in `(+ a b c d e)`.

`(split-to-subexprs vcs &key temps (min-weight 0) (gen-tmp #'symath::gen-tmp-var) subst-self)` - extract most common subexpressions from equation system `vcs` and replace them to temporary variables.
The equation system must be in alist form: `((var1 . expr1) (var2 . expr2) ...)`. The function will return an extended equation system like: `((tmp1 . tmp-expr1) (tmp2 . tmp-expr2) ... (var1 . expr1) (var2 . expr2) ...)` where `tmp1...` are symbols (by default) representing temporary variables, holding subexpressions, which occurs multiple times in the original system. The second returned values will be a list of temporary variables. The `min-weight` parameter is a minimal "weight" of the expression to be moved into a temporary variable (see `symath::expr-weight` for details). The `gen-tmp` parameter is a function of single argument, which must return an unique object (a symbol, for example) which will be a new temporary variable. The function will get an expression as a parameter. If `subst-self` is `T` the function will substitute original variables in expressions, which has the same subexpressions inside. **WARNING:** Use with great caution if your equation system is self-depended or can have multiple assignments for the same variable.

The default `gen-tmp` will be the function inside `symath` package, which returns unique interned symbols with names like `symath::tmpXXXX` where `XXXX` is an unique number, increasing each time. To reset the counter before execution of your program, use the following macro:
```
(with-var-cnt-reset
  (split-to-subexprs ...))
```

## How it works

The simplification algorithm is rather complicated, it contains the following main steps:

1. All vectoir/matrix operations (if any) are implemented, producing final vector/marix, and each element is simplified separately
2. Templates are applied
3. The expression is normalized - all subtractions like (- a b c ...) are replaced to (+ a (* -1 (+ b c ...))) and divisions like (/ a b) to (* a (expt b -1))
4. All functions with numeric arguments are computed, where possible
5. Term reduction is performed for divisions and subtraction
6. Bracketing common factors performed. Most common factors are extracted first
7. Repeat 4-6 until the result is stabilized
6. The result is denormailzed - all subtractions and divisions are returned back

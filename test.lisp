(ql:quickload :clunit)

(defpackage :symath-test
  (:use cl symath clunit))

(in-package :symath-test)

(defsuite symath ())

(deftest equal-args (symath)
  (assert-true (symath::equal-args '(1 2.0 3) '(1 2 3)))
  (assert-true (symath::equal-args '((+ x y) 2 3) '((+ y x) 2 3.0)))
  (assert-false (symath::equal-args '(1 2) '(1 2 3)))
  (assert-false (symath::equal-args '(3 1 2) '(1 2 3))))

(deftest equal-arg-sets (symath)
  (assert-true (symath::equal-arg-sets '(3 2 1) '(1 2.0 3)))
  (assert-true (symath::equal-arg-sets '((+ x y) 2 3) '(3 (+ y x) 2)))
  (assert-false (symath::equal-arg-sets '(1 2) '(1 2 3))))

(deftest equal-expr ((symath) (equal-args equal-arg-sets))
  (assert-true (symath::equal-expr '(+ x y z)
                                   '(+ y z x)))
  (assert-true (symath::equal-expr '(* x y z)
                                   '(* y z x)))
  (assert-true (symath::equal-expr '(- x y z)
                                   '(- x z y)))
  (assert-true (symath::equal-expr '(/ x y z)
                                   '(/ x z y)))
  (assert-false (symath::equal-expr '(+ x y z)
                                    '(+ y z x x)))
  (assert-false (symath::equal-expr '(* x y z)
                                    '(* y z)))
  (assert-true (symath::equal-expr '(+ x (* y z))
                                   '(+ x (* z y))))
  (assert-true (symath::equal-expr '(+ x (sqrt (* y z)))
                                   '(+ x (sqrt (* z y)))))
  (assert-true (symath::equal-expr #2A(((+ iv1 iv2 8) 2.0 3) (4 v2 6) (7 8 9))
                                   #2A(((+ iv2 8 iv1) 2 3) (4 v2 6) (7 8 9))))
  (assert-false (symath::equal-expr #2A(((+ iv1 iv2 8) 2 3) (4 v1 6) (7 8.0 9))
                                    #2A(((+ iv2 8 iv1) 2 3) (4 v2 6) (7 8 9))))
  (assert-false (symath::equal-expr #2A(((+ iv1 iv1 8) 2 3) (4 v2 6) (7 8 9))
                                    #2A(((+ iv2 8 iv1) 2 3) (4 v2 6) (7 8 9)))))

(deftest equal-expr-1 ((symath) (equal-expr))
  (assert-true (symath::equal-expr-1 '(+ 1 x)
                                     '(+ -1 (* -1 x))))
  (assert-true (symath::equal-expr-1 '(+ x 2)
                                      '(+ -2 (* -1 x))))
  (assert-true (symath::equal-expr-1 '(+ x 2)
                                     '(* -1 (+ 2 x))))
  (assert-true (symath::equal-expr-1 '(- x 2)
                                     '(- 2 x)))
  (assert-true (symath::equal-expr-1 '(- a b c d)
                                     '(+ b c d (* -1 a))))
  (assert-true (symath::equal-expr-1 '(- x a b c)
                                     '(- (+ a b c) x))))

(deftest collect-exprs ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::collect-exprs '(+ x (+ y (+ z 1) a) b))
                                   '(+ y z x 1 a b)))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::collect-exprs '(* x (* y (* z 1) a) b))
                                   '(* y z x 1 a b))))
(deftest sqr ((symath) (equal-expr))
 (assert-true (symath::equal-expr (symath::sqr 'x) '(expt x 2)))
 (assert-true (symath::equal-expr (symath::sqr #(x y z)) '(+ (expt x 2) (expt y 2) (expt z 2)))))

(deftest extract-nums ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::extract-nums '(+ 0 a b c))
                                   '(+ a b c)))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::extract-nums '(+ 1 a 2 b 3 c))
                                   '(+ a b c 6)))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::extract-nums '(* 5 a 2 b 3 c))
                                   '(* a b c 30)))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::extract-nums '(* 7 a (+ 2 3) b c))
                                   '(* a b c 35)))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::extract-nums '(* 7 a (+ 2 (sqrt (* 4 4)) b c)))
                                   '(* a 7 (+ 6.0 b c))))
  (assert-true (symath::equal (symath::math-rec-funcall #'symath::extract-nums '(* 0 a (+ 2 (sqrt (* 4 4)) b c)))
                              0))
  (assert-true (equal (symath::math-rec-funcall #'symath::extract-nums '(* 1 a b c))
                      '(* a b c)))
  (assert-true (equal (symath::math-rec-funcall #'symath::extract-nums '(* (+ 3 (* 3 -1)) a (+ 2 (sqrt (* 4 4)) b c)))
                      0)))

(deftest sqr ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::sqr 'x) '(expt x 2)))
  (assert-true (symath::equal-expr (symath::sqr #(x y z)) '(+ (expt x 2) (expt y 2) (expt z 2)))))

(deftest norm-expr ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::norm-expr '(+ a b (* c d)))
                                   '(+ a b (* c d))))
  (assert-true (symath::equal-expr (symath::norm-expr '(- a b))
                                   '(+ a (* -1 b))))
  (assert-true (symath::equal-expr (symath::norm-expr '(- a b c d))
                                   '(+ a (* -1 (+ b c d)))))
  (assert-true (symath::equal-expr (symath::norm-expr '(/ a b))
                                   '(* a (expt b -1))))
  (assert-true (symath::equal-expr (symath::norm-expr '(/ a b c d))
                                   '(* a (expt (* b c d) -1))))
  (assert-true (symath::equal-expr (symath::norm-expr '(sqrt (+ a b c d)))
                                   '(expt (+ a b c d) 1/2))))

(deftest denorm-expr-expt ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(expt x -5))
                                   '(/ 1 (expt x 5))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(* (expt x -5) y))
                                   '(/ y (expt x 5))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(* (expt x -5) (expt a -3) y))
                                   '(/ y (* (expt x 5) (expt a 3)))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(* (expt x -5) (expt a -3)))
                                   '(/ 1 (* (expt x 5) (expt a 3)))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(* (expt x -5) (expt a -3) y z))
                                   '(/ (* y z) (* (expt x 5) (expt a 3)))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(* (expt x -5) b y z))
                                   '(/ (* b y z) (expt x 5))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(* (expt x -5) b y (expt q -1) z))
                                   '(/ (* b y z) (* q (expt x 5)))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(* (expt x -1) b y (expt q -1) z))
                                   '(/ (* b y z) (* q x))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::denorm-expr-expt '(* (expt x -1) b y (expt q -1) (expt z 2)))
                                   '(/ (* b y (expt z 2)) (* q x)))))

(deftest denorm-expr-zop ((symath) (equal-expr))
  (assert-true (equal (symath::denorm-expr-zop '(* x))
                      'x))
  (assert-true (equal (symath::denorm-expr-zop '(+ x))
                      'x))
  (assert-true (symath::equal-expr (symath::denorm-expr-zop '(+ x y))
                                   '(+ x y)))
  (assert-true (symath::equal-expr (symath::denorm-expr-zop '(* x y z))
                                   '(* x y z))))

(deftest denorm-expr-minus ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::denorm-expr-minus '(+ x (* -1 y)))
                                   '(- x y)))
  (assert-true (symath::equal-expr (symath::denorm-expr-minus '(+ x (* -5 y)))
                                   '(- x (* 5 y))))
  (assert-true (symath::equal-expr (symath::denorm-expr-minus '(+ x y z (* -1 a b c d)))
                                   '(- (+ x y z) (* a b c d))))
  (assert-true (symath::equal-expr (symath::denorm-expr-minus '(+ x y z (* -1 a) (* -1 b)))
                                   '(- (+ x y z) a b)))
  (assert-true (symath::equal-expr (symath::denorm-expr-minus '(+ x (* -1 a) (* -1 b)))
                                   '(- x a b)))
  (assert-true (symath::equal-expr (symath::denorm-expr-minus '(+ (* -1 a) (* -1 b)))
                                   '(* -1 (+ a b)))))

(deftest denorm-expr-minus1 ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::denorm-expr-minus1 '(+ -5 x y z))
                                   '(- (+ x y z) 5)))
  (assert-true (symath::equal-expr (symath::denorm-expr-minus1 '(+ -5 x))
                                   '(- x 5))))

(deftest denorm-expr-plus-n ((symath) (equal-expr))
  (assert-true (equal (symath::denorm-expr-plus-n '(+ 5 x y z))
                      '(+ x y z 5)))
  (assert-true (equal (symath::denorm-expr-plus-n '(+ 5 x))
                      '(+ x 5))))

(deftest expand-expt1 ((symath) (equal-expr))
  (assert-true (equal (symath::expand-expt1 '(expt (* x y z) 2))
                      '(* (expt x 2) (expt y 2) (expt z 2)))))

(deftest expand-expt2 ((symath) (equal-expr))
  (assert-true (equal (symath::expand-expt2 '(expt x (+ y z)))
                      '(* (expt x y) (expt x z)))))

(deftest collect-expt ((symath) (equal-expr equal-expr-1 extract-nums collect-exprs))
  (assert-true (symath::equal-expr (symath::chain-func-rec (symath::extract-nums symath::collect-exprs symath::collect-expt) '(* x (expt x y) (expt x z) q))
                                   '(* (expt x (+ y z 1)) q)))
  (assert-true (symath::equal-expr (symath::chain-func-rec (symath::extract-nums symath::collect-exprs symath::collect-expt) '(expt (expt x z) y))
                                   '(expt x (* y z))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::collect-expt '(expt (* x y) z))
                                   '(* (expt x z) (expt y z))))
  (assert-true (symath::equal-expr (symath::chain-func-rec (symath::extract-nums symath::collect-exprs symath::collect-expt) '(* (+ (* -1 a) b) (expt (+ a (* -1 b)) -1)))
                                   -1))
  (assert-true (let ((res (symath::chain-func-rec (symath::extract-nums symath::collect-exprs symath::collect-expt) '(* (+ (* -1 a) b) (+ a (* -1 b))))))
                 (or (symath::equal-expr res '(* -1 (expt (+ (* -1 a) b) 2)))
                     (symath::equal-expr res '(* -1 (expt (+ (* -1 b) a) 2))))))
  (assert-true (let ((res (symath::chain-func-rec (symath::extract-nums symath::collect-exprs symath::collect-expt) '(* (+ (* -1 a) b) (+ a (* -1 b)) 10))))
                 (or (symath::equal-expr res '(* -10 (expt (+ (* -1 a) b) 2)))
                     (symath::equal-expr res '(* -10 (expt (+ (* -1 b) a) 2))))))
  (assert-true (let ((res (symath::chain-func-rec (symath::extract-nums symath::collect-exprs symath::collect-expt) '(* (+ (* -1 a) b) (+ a (* -1 b)) 10 (+ a b)))))
                 (or (symath::equal-expr res '(* -10 (expt (+ (* -1 a) b) 2) (+ a b)))
                     (symath::equal-expr res '(* -10 (expt (+ (* -1 b) a) 2) (+ a b)))))))

(deftest expand-mul ((symath) (equal-expr collect-exprs))
  (assert-true (symath::equal-expr (symath::expand-mul '(* (+ x y) (+ a b) q))
                                   '(+ (* x a q) (* x b q) (* y a q) (* y b q))))
  (assert-true (symath::equal-expr (symath::math-rec-funcall #'symath::expand-mul '(* (+ x y) (+ a (* b (+ c d))) q))
                                   '(+ (* x a q) (* x b c q) (* x b d q) (* y a q) (* y b c q) (* y b d q)))))


(deftest collect-common-nums ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::collect-common-nums '(+ (* 2 x) (* 4 y z) (* 6 z)))
                                   '(* 2 (+ x (* 2 y z) (* 3 z)))))
  (assert-true (symath::equal-expr (symath::collect-common-nums '(+ (* 2 x) (* 4.0 y z) (* 6 z)))
                                   '(* 2 (+ x (* 2 y z) (* 3 z))))))


(deftest expand-int-expt ((symath) (equal-expr))
  (let ((symath::*max-degree-expansion* 3))
    (assert-true (symath::equal-expr (symath::expand-int-expt '(expt (+ x y) 3))
                                     '(* (+ x y) (+ x y) (+ x y))))
    (assert-true (symath::equal-expr (symath::expand-int-expt '(expt (+ x y) -3))
                                     '(* (expt (+ x y) -1) (expt (+ x y) -1) (expt (+ x y) -1))))
    (assert-true (symath::equal-expr (symath::expand-int-expt '(expt (+ x y) 5))
                                     '(expt (+ x y) 5)))
    (assert-true (symath::equal-expr (symath::expand-int-expt '(expt x 2))
                                     '(expt x 2)))))

(deftest extract-subexpr ((symath) (equal-expr norm-expr extract-nums expand-mul collect-exprs collect-expt expand-expt1 expand-expt2 expand-int-expt))
  (assert-true (multiple-value-bind (n e1 e2) (extract-subexpr '(+ (sqrt x) 1) 'x)
                 (and (equal n 1/2)
                      (equal e1 1)
                      (equal e2 1))))
  (assert-true (multiple-value-bind (n e1 e2) (extract-subexpr '(+ (sqrt (+ x y)) 1) '(+ x y))
                  (and (equal n 1/2)
                       (equal e1 1)
                       (equal e2 1))))
  (assert-true (multiple-value-bind (n e1 e2) (extract-subexpr '(+ (* (- a b) x) b (* -1 a)) '(- a b))
                  (and (equal n 1)
                       (symath::equal-expr e1 '(+ -1 X))
                       (equal e2 0)))))

(deftest collect-common ((symath) (equal-expr equal-expr-1 extract-nums collect-expt collect-exprs extract-subexpr))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (* x y) (* x z) a))
                                   '(+ (* x (+ y z)) a)))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (* (+ x 10) y) (* (+ x 10) z) a))
                                   '(+ (* (+ x 10) (+ y z)) a)))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (* (+ x 10) (+ y 20)) (+ y 20) a))
                                   '(+ (* (+ y 20) (+ x 11)) a)))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (* x a) (* x b) (* y a) (* y b)))
                                   '(* (+ x y) (+ a b))))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (* (expt x 2) a) (* (expt x 2) b) (* y a) (* y b)))
                                   '(* (+ (expt x 2) y) (+ a b))))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (expt x 2) (expt x 3) (expt x 4)))
                                   '(* (expt x 2) (+ 1 (* x (+ x 1))))))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (expt x 2) (* a (expt x 3)) (* (expt x 4) b)))
                                   '(* (expt x 2) (+ 1 (* x (+ a (* x b)))))))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (* x a q) (* x b q) (* y a q) (* y b q)))
                                   '(* q (+ x y) (+ a b))))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (expt x 2) (* -1 (expt x 2))))
                                   0))
  (assert-true (symath::equal-expr (symath::denorm-expr (symath::collect-common '(+ (* (+ x (* -1 y)) a) (* (+ y (* -1 x)) b) c)))
                                  '(+ (* (- x y) (- a b)) c)))
  (assert-true (symath::equal-expr (symath::denorm-expr (symath::collect-common '(+ (* (+ x (* -1 y)) a) (* (+ (* y 2) (* -2 x)) b) c)))
                                   '(+ (* (- x y) (- a (* 2 b))) c)))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (* x (+ a b c)) a b c d e))
                                   '(+ (* (+ a b c) (+ 1 x)) d e)))
  (assert-true (symath::equal-expr (symath::collect-common '(+ (* x (+ a (* -1 b))) (* -1 a) b))
                                   '(* (+ a (* -1 b)) (+ x -1))))
  (assert-true (or (symath::equal-expr (symath::collect-common '(+ (* (expt (+ x (* -1 y)) 2) aaa) (* (expt (+ y (* -1 x)) 2) bbb)))
                                       '(* (EXPT (+ X (* -1 Y)) 2) (+ BBB AAA)))
                   (symath::equal-expr (symath::collect-common '(+ (* (expt (+ x (* -1 y)) 2) aaa) (* (expt (+ y (* -1 x)) 2) bbb)))
                                       '(* (EXPT (+ Y (* -1 X)) 2) (+ BBB AAA))))))

(deftest array-multiply ((symath) (equal-expr))
  (assert-true (symath::equal-expr (array-multiply 'x 'y)
                                   '(* x y)))
  (assert-true (symath::equal-expr (array-multiply 'x #(a b c))
                                   #((* x a) (* x b) (* x c))))
  (assert-true (symath::equal-expr (array-multiply '(* x y)  #(a b c))
                                   #((* x y a) (* x y b) (* x y c))))
  (assert-true (symath::equal-expr (array-multiply 'x  #((* a q) (* b q) (* c q)))
                                   #((* x a q) (* x b q) (* x c q))))
  (assert-true (symath::equal-expr (array-multiply '(* x y)  #((* a q) (* b q) (* c q)))
                                   #((* x y a q) (* x y b q) (* x y c q))))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                       '((x11 x12 x13)
                                                         (x21 x22 x23)
                                                         (x31 x32 x33)))
                                                   'a)
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((* x11 a) (* x12 a) (* x13 a))
                                          ((* x21 a) (* x22 a) (* x23 a))
                                          ((* x31 a) (* x32 a) (* x33 a))))))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                        '(((* x11 q) (* x12 q) (* x13 q))
                                                          ((* x21 q) (* x22 q) (* x23 q))
                                                          ((* x31 q) (* x32 q) (* x33 q))))
                                                   'a)
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((* x11 a q) (* x12 a q) (* x13 a q))
                                          ((* x21 a q) (* x22 a q) (* x23 a q))
                                          ((* x31 a q) (* x32 a q) (* x33 a q))))))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                        '(((* x11 q) (* x12 q) (* x13 q))
                                                          ((* x21 q) (* x22 q) (* x23 q))
                                                          ((* x31 q) (* x32 q) (* x33 q))))
                                                   '(* a b))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((* x11 b a q) (* x12 b a q) (* x13 b a q))
                                          ((* x21 b a q) (* x22 b a q) (* x23 b a q))
                                          ((* x31 b a q) (* x32 b a q) (* x33 b a q))))))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                        '((x11 x12 x13)
                                                          (x21 x22 x23)
                                                          (x31 x32 x33)))
                                                   '(* a b))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((* x11 a b) (* x12 a b) (* x13 a b))
                                          ((* x21 a b) (* x22 a b) (* x23 a b))
                                          ((* x31 a b) (* x32 a b) (* x33 a b))))))
  (assert-true (symath::equal (array-multiply 0 (make-array '(3 3)
                                                  :initial-contents
                                                   '((x11 x12 x13)
                                                     (x21 x22 x23)
                                                     (x31 x32 x33))))
                              0))
  (assert-true (equal (array-multiply (make-array '(3 3)
                                        :initial-contents
                                         '((x11 x12 x13)
                                           (x21 x22 x23)
                                           (x31 x32 x33)))
                                      0)
                      0))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                        '((x11 x12 x13)
                                                          (x21 x22 x23)
                                                          (x31 x32 x33)))
                                                   1)
                                   (make-array '(3 3)
                                     :initial-contents
                                        '((x11 x12 x13)
                                          (x21 x22 x23)
                                          (x31 x32 x33)))))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                        '((x11 x12 x13)
                                                          (x21 x22 x23)
                                                          (x31 x32 x33)))
                                                   (make-array '(3 3)
                                                     :initial-contents
                                                        '((y11 y12 y13)
                                                          (y21 y22 y23)
                                                          (y31 y32 y33))))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((+ (* x11 y11) (* x12 y21) (* x13 y31))
                                           (+ (* x11 y12) (* x12 y22) (* x13 y32))
                                           (+ (* x11 y13) (* x12 y23) (* x13 y33)))
                                          ((+ (* x21 y11) (* x22 y21) (* x23 y31))
                                           (+ (* x21 y12) (* x22 y22) (* x23 y32))
                                           (+ (* x21 y13) (* x22 y23) (* x23 y33)))
                                          ((+ (* x31 y11) (* x32 y21) (* x33 y31))
                                           (+ (* x31 y12) (* x32 y22) (* x33 y32))
                                           (+ (* x31 y13) (* x32 y23) (* x33 y33)))))))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                        '(((* x11 a) (* x12 a) (* x13 a))
                                                          ((* x21 a) (* x22 a) (* x23 a))
                                                          ((* x31 a) (* x32 a) (* x33 a))))
                                                   (make-array '(3 3)
                                                     :initial-contents
                                                        '((y11 y12 y13)
                                                          (y21 y22 y23)
                                                          (y31 y32 y33))))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((+ (* a x11 y11) (* a x12 y21) (* a x13 y31))
                                           (+ (* a x11 y12) (* a x12 y22) (* a x13 y32))
                                           (+ (* a x11 y13) (* a x12 y23) (* a x13 y33)))
                                          ((+ (* a x21 y11) (* a x22 y21) (* a x23 y31))
                                           (+ (* a x21 y12) (* a x22 y22) (* a x23 y32))
                                           (+ (* a x21 y13) (* a x22 y23) (* a x23 y33)))
                                          ((+ (* a x31 y11) (* a x32 y21) (* a x33 y31))
                                           (+ (* a x31 y12) (* a x32 y22) (* a x33 y32))
                                           (+ (* a x31 y13) (* a x32 y23) (* a x33 y33)))))))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                        '(((* x11 a) (* x12 a) (* x13 a))
                                                          ((* x21 a) (* x22 a) (* x23 a))
                                                          ((* x31 a) (* x32 a) (* x33 a))))
                                                   (make-array '(3 3)
                                                     :initial-contents
                                                        '(((* y11 b) (* y12 b) (* y13 b))
                                                          ((* y21 b) (* y22 b) (* y23 b))
                                                          ((* y31 b) (* y32 b) (* y33 b)))))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((+ (* a b x11 y11) (* a b x12 y21) (* a b x13 y31))
                                           (+ (* a b x11 y12) (* a b x12 y22) (* a b x13 y32))
                                           (+ (* a b x11 y13) (* a b x12 y23) (* a b x13 y33)))
                                          ((+ (* a b x21 y11) (* a b x22 y21) (* a b x23 y31))
                                           (+ (* a b x21 y12) (* a b x22 y22) (* a b x23 y32))
                                           (+ (* a b x21 y13) (* a b x22 y23) (* a b x23 y33)))
                                          ((+ (* a b x31 y11) (* a b x32 y21) (* a b x33 y31))
                                           (+ (* a b x31 y12) (* a b x32 y22) (* a b x33 y32))
                                           (+ (* a b x31 y13) (* a b x32 y23) (* a b x33 y33)))))))
  (assert-true (symath::equal-expr (array-multiply (make-array '(3 3)
                                                     :initial-contents
                                                        '((x11 x12 x13)
                                                          (x21 x22 x23)
                                                          (x31 x32 x33)))
                                                   (make-array '(3 3)
                                                     :initial-contents
                                                        '(((* y11 b) (* y12 b) (* y13 b))
                                                          ((* y21 b) (* y22 b) (* y23 b))
                                                          ((* y31 b) (* y32 b) (* y33 b)))))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((+ (* b x11 y11) (* b x12 y21) (* b x13 y31))
                                           (+ (* b x11 y12) (* b x12 y22) (* b x13 y32))
                                           (+ (* b x11 y13) (* b x12 y23) (* b x13 y33)))
                                          ((+ (* b x21 y11) (* b x22 y21) (* b x23 y31))
                                           (+ (* b x21 y12) (* b x22 y22) (* b x23 y32))
                                           (+ (* b x21 y13) (* b x22 y23) (* b x23 y33)))
                                          ((+ (* b x31 y11) (* b x32 y21) (* b x33 y31))
                                           (+ (* b x31 y12) (* b x32 y22) (* b x33 y32))
                                           (+ (* b x31 y13) (* b x32 y23) (* b x33 y33))))))))

(deftest calc-arrays ((symath) (equal-expr array-multiply))
  (assert-true (symath::equal-expr (symath::calc-arrays `(+ ,(make-array '(3 3)
                                                               :initial-contents
                                                                  '((x11 x12 x13)
                                                                    (x21 x22 x23)
                                                                    (x31 x32 x33)))
                                                           ,(make-array '(3 3)
                                                               :initial-contents
                                                                  '((y11 y12 y13)
                                                                    (y21 y22 y23)
                                                                    (y31 y32 y33)))))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((+ x11 y11) (+ x12 y12) (+ x13 y13))
                                          ((+ x21 y21) (+ x22 y22) (+ x23 y23))
                                          ((+ x31 y31) (+ x32 y32) (+ x33 y33))))))
  (assert-true (symath::equal-expr (symath::calc-arrays `(- ,(make-array '(3 3)
                                                               :initial-contents
                                                                  '((x11 x12 x13)
                                                                    (x21 x22 x23)
                                                                    (x31 x32 x33)))
                                                            ,(make-array '(3 3)
                                                               :initial-contents
                                                                  '((y11 y12 y13)
                                                                    (y21 y22 y23)
                                                                    (y31 y32 y33)))
                                                            ,(make-array '(3 3)
                                                               :initial-contents
                                                                  '((z11 z12 z13)
                                                                    (z21 z22 z23)
                                                                    (z31 z32 z33)))))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((- x11 y11 z11) (- x12 y12 z12) (- x13 y13 z13))
                                          ((- x21 y21 z21) (- x22 y22 z22) (- x23 y23 z23))
                                          ((- x31 y31 z31) (- x32 y32 z32) (- x33 y33 z33))))))
  (assert-true (symath::equal-expr (symath::calc-arrays `(* a ,(make-array '(3 3)
                                                                 :initial-contents
                                                                    '((x11 x12 x13)
                                                                      (x21 x22 x23)
                                                                      (x31 x32 x33)))
                                                              ,(make-array '(3 3)
                                                                 :initial-contents
                                                                    '((y11 y12 y13)
                                                                      (y21 y22 y23)
                                                                      (y31 y32 y33)))))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((+ (* a x11 y11) (* a x12 y21) (* a x13 y31))
                                           (+ (* a x11 y12) (* a x12 y22) (* a x13 y32))
                                           (+ (* a x11 y13) (* a x12 y23) (* a x13 y33)))
                                          ((+ (* a x21 y11) (* a x22 y21) (* a x23 y31))
                                           (+ (* a x21 y12) (* a x22 y22) (* a x23 y32))
                                           (+ (* a x21 y13) (* a x22 y23) (* a x23 y33)))
                                          ((+ (* a x31 y11) (* a x32 y21) (* a x33 y31))
                                           (+ (* a x31 y12) (* a x32 y22) (* a x33 y32))
                                           (+ (* a x31 y13) (* a x32 y23) (* a x33 y33)))))))
  (assert-true (symath::equal-expr (symath::calc-arrays `(/ ,(make-array '(3 3)
                                                               :initial-contents
                                                                  '((x11 x12 x13)
                                                                    (x21 x22 x23)
                                                                    (x31 x32 x33)))
                                                            a))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((* x11 (expt a -1)) (* x12 (expt a -1)) (* x13 (expt a -1)))
                                          ((* x21 (expt a -1)) (* x22 (expt a -1)) (* x23 (expt a -1)))
                                          ((* x31 (expt a -1)) (* x32 (expt a -1)) (* x33 (expt a -1)))))))
  (assert-true (symath::equal-expr (symath::calc-arrays `(* a b ,(make-array '(3 3)
                                                                   :initial-contents
                                                                      '((x11 x12 x13)
                                                                        (x21 x22 x23)
                                                                        (x31 x32 x33)))
                                                                ,(make-array '(3 3)
                                                                   :initial-contents
                                                                      '((y11 y12 y13)
                                                                        (y21 y22 y23)
                                                                        (y31 y32 y33)))))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((+ (* a b x11 y11) (* a b x12 y21) (* a b x13 y31))
                                           (+ (* a b x11 y12) (* a b x12 y22) (* a b x13 y32))
                                           (+ (* a b x11 y13) (* a b x12 y23) (* a b x13 y33)))
                                          ((+ (* a b x21 y11) (* a b x22 y21) (* a b x23 y31))
                                           (+ (* a b x21 y12) (* a b x22 y22) (* a b x23 y32))
                                           (+ (* a b x21 y13) (* a b x22 y23) (* a b x23 y33)))
                                          ((+ (* a b x31 y11) (* a b x32 y21) (* a b x33 y31))
                                           (+ (* a b x31 y12) (* a b x32 y22) (* a b x33 y32))
                                           (+ (* a b x31 y13) (* a b x32 y23) (* a b x33 y33)))))))
  (assert-true (symath::equal-expr (symath::calc-arrays `(expt ,(make-array '(3 3)
                                                                  :initial-contents
                                                                     '((x11 x12 x13)
                                                                       (x21 x22 x23)
                                                                       (x31 x32 x33)))
                                                               2))
                                   (make-array '(3 3)
                                     :initial-contents
                                        '(((+ (* x11 x11) (* x12 x21) (* x13 x31))
                                           (+ (* x11 x12) (* x12 x22) (* x13 x32))
                                           (+ (* x11 x13) (* x12 x23) (* x13 x33)))
                                          ((+ (* x21 x11) (* x22 x21) (* x23 x31))
                                           (+ (* x21 x12) (* x22 x22) (* x23 x32))
                                           (+ (* x21 x13) (* x22 x23) (* x23 x33)))
                                          ((+ (* x31 x11) (* x32 x21) (* x33 x31))
                                           (+ (* x31 x12) (* x32 x22) (* x33 x32))
                                           (+ (* x31 x13) (* x32 x23) (* x33 x33)))))))
  (assert-true (symath::equal-expr (symath::calc-arrays `(expt ,#(x y z) 2))
                                   '(+ (expt x 2) (expt y 2) (expt z 2))))
  (assert-true (symath::equal-expr (symath::calc-arrays `(expt ,#(x y z) 3))
                                   #((* x (+ (expt x 2) (expt y 2) (expt z 2)))
                                     (* y (+ (expt x 2) (expt y 2) (expt z 2)))
                                     (* z (+ (expt x 2) (expt y 2) (expt z 2)))))))

(deftest collect-same-expts ((symath) (equal-expr))
  (assert-true (symath::equal-expr (symath::collect-same-expts '(* a b (expt x 2) (expt y 2)))
                                   '(* a b (expt (* x y) 2))))
  (assert-true (symath::equal-expr (symath::collect-same-expts '(* (expt x 2) (expt y 2)))
                                   '(expt (* x y) 2)))
  (assert-true (symath::equal-expr (symath::collect-same-expts '(* (expt x 2) (expt a 3) (expt y 2) (expt b 3)))
                                   '(* (expt (* x y) 2) (expt (* a b) 3)))))

(deftest simplify ((symath) (equal-expr collect-same-expts calc-arrays collect-common collect-common-nums
                             collect-expt expand-mul expand-expt2 expand-expt1 denorm-expr-plus-n denorm-expr-minus1
                             denorm-expr-minus denorm-expr-zop denorm-expr-expt norm-expr extract-nums collect-exprs)
                            expand-int-expt)
  (symath::clear-simplify-cache)
  (assert-true (symath::equal-expr (simplify '(/ (* (+ x 1) (+ y 1)) (* (+ y 1) (+ z 1) (+ x 1))))
                                   '(/ 1 (+ z 1))))
  (assert-true (symath::equal-expr (simplify '(+ (* x (+ a b c)) a b c d))
                                   '(+ (* (+ x 1) (+ a b c)) d)))
  (assert-true (equal (simplify '(- x (+ x 1)))
                      -1))
  (assert-true (symath::equal-expr (simplify '(/ (- xmi xma) (* (- xma xmi) (- yma ymi) (- zma zmi))))
                                   '(/ -1 (* (- zma zmi) (- yma ymi)))))
  (assert-true (symath::equal-expr (simplify '(+ (expt x b) (expt x b)))
                                   '(* 2 (expt x b)))))

(run-suite 'symath)

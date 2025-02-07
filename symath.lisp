(defpackage :symath
  (:use cl alexandria)
  (:export simplify array-multiply extract-subexpr get-polynome-cfs replace-subexpr split-to-subexprs with-var-cnt-reset))

(in-package :symath)

(defun mequal (a b)
  (or (equal a b)
      (and (numberp a)
           (numberp b)
           (= a b))))

(defun isfunc (f e) (and (listp e) (equal (car e) f)))

(defun filter-func (f l &optional invrs)
  (labels ((flt (e) (isfunc f e)))
    (if invrs
        (remove-if #'flt l)
        (remove-if-not #'flt l))))

(defun all-arraysp (l)
  "Check if all elements in l are arrays"
  (and (arrayp (car l))
       (or (not (cdr l))
           (all-arraysp (cdr l)))))

(defun equal-dimsp (l)
  "Check if all arrays in l has the same dimensions"
  (and (all-arraysp l)
       (apply #'= (mapcar #'array-rank l))
       (notany #'null
               (apply #'mapcar
                 (cons #'=
                       (mapcar #'array-dimensions l))))))

(defun map-array (lam &rest arrs)
  "Like #'mapcar, but for arrays"
  (if (and (all-arraysp arrs)
           (equal-dimsp arrs))
    (let ((arr (make-array (array-dimensions (car arrs)))))
      (loop for i below (apply #'* (array-dimensions (car arrs)))
        do (setf (row-major-aref arr i)
                 (apply lam
                        (mapcar
                          (lambda (a)
                            (row-major-aref a i))
                          arrs))))
      arr)
    (error "Cannot map arrays with different dimensions")))

(defmethod array-multiply (x y)
  (cond ((and (not (arrayp x))
              (not (arrayp y)))
         `(* ,x ,y))
        ((not (arrayp y))
         (array-multiply y x))
        ((mequal x 0) 0)
        ((mequal x 1) y)
        ((not (arrayp x))
         (map-array
           (lambda (e)
             (cond ((or (mequal e 0) (mequal x 0)) 0)
                   ((mequal e 1) x)
                   ((mequal x 1) e)
                   ((and (isfunc '* e) (isfunc '* x)) `(* ,@(cdr e) ,@(cdr x)))
                   ((isfunc '* e) `(* ,@(cdr e) ,x))
                   ((isfunc '* x) `(* ,@(cdr x) ,e))
                   (t `(* ,x ,e))))
           y))
        ((or (> (array-rank x) 2)
             (> (array-rank y) 2))
         (error "Array multiplication for non-2D matrices is not implemented. Overload #'array-multiple to implement."))
        ((and (equal (array-rank x) 2)
              (equal (array-rank y) 2)
              (not (cdr (remove-duplicates (append (array-dimensions x) (array-dimensions y))))))
         (let ((dim (array-dimension x 0)))
           (make-array `(,dim ,dim)
             :initial-contents
               (loop for j below dim collect
                 (loop for i below dim collect
                   `(+ ,@(loop for k below dim
                           collect (let ((v1 (aref x j k))
                                         (v2 (aref y k i)))
                                     (cond ((and (isfunc '* v1) (isfunc '* v2)) `(* ,@(cdr v1) ,@(cdr v2)))
                                           ((isfunc '* v1) `(* ,@(cdr v1) ,v2))
                                           ((isfunc '* v2) `(* ,v1 ,@(cdr v2)))
                                           (t `(* ,v1 ,v2)))))))))))
        (t (error "The '*' operator can be used only for multiplication of vector/matrix to constant or to multiply square matrices. Overload #'array-multiply to extend it."))))

(defun tick ())

(defun boolp (e) (typep e 'boolean))

(defun copy-expr (e)
  (if (listp e)
      (copy-list e)
      e))

(defun mintegerp (x)
  (or (integerp x)
      (and (numberp x)
           (= (floor x) x))))

(defun all-bins (n) ;; 3 => ((0 0 0) (1 0 0) (0 1 0) (1 1 0) (0 0 1) (1 0 1) (0 1 1) (1 1 1))
  (if (> n 1)
      (mapcan
        (lambda (b)
          (list (cons 0 b)
                (cons 1 b)))
        (all-bins (- n 1)))
      `((0) (1))))

(defun all-list-decs (l &key (min 1) (max (length l)))
  (let ((ll (length l)))
    (sort
      (mapcar
        (lambda (b)
          (mapcar #'cdr
            (remove-if (lambda (x) (= (car x) 0))
                       (mapcar #'cons b l))))
        (mapcar #'cdr
          (remove-if-not
            (lambda (x) (and (>= (car x) min)
                             (<= (car x) max)))
            (mapcar (lambda (x) (cons  (count 1 x) x))
                    (all-bins ll)))))
      #'>
      :key #'length)))

(defun sqr (x)
  "Multiple argument by itself. Works for vectors too!"
  (if (arrayp x)
      (if (cdr (array-dimensions x))
          (array-multiply x x)
          `(+ ,@(map 'list (lambda (x) `(expt ,x 2)) x)))
      `(expt ,x 2)))

(defun math-rec-funcall (fn arg &optional deep)
  "Call the function fn recursively"
  (labels ((rfn (x)
             (if (equal deep 0)
                 x
                 (math-rec-funcall fn x (if deep (- deep 1))))))
    (multiple-value-call fn (cond
                              ((listp arg)
                               (cons (car arg) (mapcar #'rfn (cdr arg))))
                              ((arrayp arg)
                               (map-array #'rfn arg))
                              (t arg)))))

(defparameter *return-best-variant* nil) ;; If T, defun-stable-expr functions will return less weighted result. USE WITH CAUTION!

(defmacro defun-stable-expr (name args &rest cod)
  (let ((h (gensym))
        (r (gensym))
        (arg (gensym)))
    `(defun ,name (,arg &optional ,h)
       (let* ((,r (funcall (lambda ,args ,@cod) ,arg)))
         (if (position ,r ,h :test #'equal-expr)
             (if *return-best-variant*
                 (caar (stable-sort (mapcar #'cons ,h (mapcar #'expr-weight ,h)) #'< :key #'cdr))
                 ,r)
             (,name ,r (cons ,r ,h)))))))

(defmacro chain-func-rec (fl arg)
  "Recursively call function chain"
  (let ((rfl (reverse fl)))
    (reduce (lambda (r x)
              (if (listp x)
                  `(multiple-value-call #'math-rec-funcall #',(cadr x) ,r ,(car x))
                  `(multiple-value-call #'math-rec-funcall #',x ,r)))
            (cdr rfl)
            :initial-value (if (listp (car rfl))
                               `(multiple-value-call #'math-rec-funcall #',(cadr (car rfl)) ,arg ,(car (car rfl)))
                               `(multiple-value-call #'math-rec-funcall #',(car rfl) ,arg)))))

(defmacro def-expr-cond (&rest margs)
  "(expr-cond func-name expr-name :numberp cod :boolp cod :op (op cod) ... :defop cod)"
  (let* ((func-name (car margs))
         (e1 (cadr margs))
         (expr-name e1)
         (pars (loop for x on (cddr margs) by #'cddr collect (cons (car x) (cadr x))))
         (ops (mapcar #'cdr (remove-if-not (lambda (x) (equal (car x) :op)) pars)))
         (defop (cdar (remove-if-not (lambda (x) (equal (car x) :defop)) pars)))
         (non-ops (remove-if-not (lambda (x)
                                   (not (or (equal (car x) :op)
                                            (equal (car x) :defop)))) pars)))
     `(defun ,func-name ,(if (listp e1) e1 (list e1))
        (tick)
        (cond ,@(mapcar (lambda (x)
                          `((,(intern (symbol-name (car x))) ,expr-name) ,(cdr x)))
                        non-ops)
              ,@(if (not (or (position :arrayp non-ops)
                             (position :vectorp non-ops)))
                    `(((arrayp ,expr-name)
                       (map-array #',func-name ,expr-name))))
              ((listp ,expr-name)
               (cond
                ,@(mapcar (lambda (x)
                            (let ((op (car x))
                                  (cod (cadr x)))
                                `(,(if (listp op)
                                      `(or ,@(mapcar (lambda (x)
                                                       `(equal (car ,expr-name) ',x)) op))
                                      `(equal (car ,expr-name) ',op))
                                   ,cod)))
                          ops)
                ,@(if defop `((t ,defop)) `((t ,expr-name)))))
              (t ,expr-name)))))

; (defun equal-expr (e1 e2)  ;; Just to suppress warnings about undefined function
;   (declare (ignore e1 e2))
;   (error "STUB"))

(defun equal-args (l1 l2)
  "Check if l1 and l2 are lists of the same expressions"
  (and (equal (length l1) (length l2))
       (labels ((e-a (l1 l2)
                  (or (and (not l1) (not l2))
                      (and (equal-expr (car l1) (car l2))
                           (e-a (cdr l1) (cdr l2))))))
         (e-a l1 l2))))

(defparameter *equarg-uniq-sym* (gensym))

(defun equal-arg-sets (l1 l2)
  "Check if l1 and l2 are lists of the same sets of expressions"
  (and (equal (length l1) (length l2))
    (let* ((tl2 (copy-list l2)))
      (map nil
        (lambda (e)
          (let ((pos (position e tl2 :test #'equal-expr)))
            (if pos
                (setf (elt tl2 pos) *equarg-uniq-sym*)
                (return-from equal-arg-sets nil))))
        l1)
      (equal *equarg-uniq-sym* (car (remove-duplicates tl2))))))

;; equal-expr and equal-expr-1 are not checks for true equality in math sense,
;; expressions must have similar forms to be treated as equal.
;; For example, exprs (+ a b c) and (+ a (+ b c)) will not be found as equal,
;; but (+ a b c) and (+ b c a) will. The equal-expr-1 a bit more tricky, but
;; also don't trust it too much.
;; A better (but much slower) way to check exprs equality is:
;;   (or (equal 1 (symplify '(/ e1 e2))) (equal 0 (symplify '(- e1 e2))))

(defun equal-expr (e1 e2)
  "Check if e1 and e2 a the same exprs"
  (or (tree-equal e1 e2) ;; Quick result if we are lucky
      (mequal e1 e2)
      (and (equal (type-of e1) (type-of e2))
           (cond ((symbolp e1) (and (equal (symbol-name e1) (symbol-name e2)) ;; It is possible to have different symbols like with same names!
                                    (equal (symbol-package e1) (symbol-package e2))))
                 ((listp e1) (and (equal (car e1) (car e2))
                                  (and (= (length e1) (length e2))
                                       (if (or (equal (car e1) '*) (equal (car e1) '+))
                                           (equal-arg-sets (cdr e1) (cdr e2))
                                           (if (or (equal (car e1) '/) (equal (car e1) '-))
                                               (and (equal-expr (cadr e1) (cadr e2))
                                                    (equal-arg-sets (cddr e1) (cddr e2)))
                                               (equal-args (cdr e1) (cdr e2)))))))
                 ((arrayp e1)
                  (and (equal (array-dimensions e1) (array-dimensions e2))
                       (loop for i below (apply #'* (array-dimensions e1))
                             when (not (equal-expr (row-major-aref e1 i) (row-major-aref e2 i))) return nil
                             finally (return t))))
                 (t nil)))))

(defun equal-expr-1 (e1 e2)
  "Check if e1 is -e2"
  (or (and (numberp e1)
           (numberp e2)
           (= e1 (* -1 e2)))
      (or (and (listp e1) (listp e2)
               (equal (car e1) (car e2))
               (equal (length e1) (length e2))
               (equal (car e1) '+)
               (labels ((eq-sum (l1 l2)
                          (or (not l1)
                              (let ((pos (position (car l1) l2 :test #'equal-expr-1)))
                                (and pos
                                     (eq-sum (cdr l1)
                                             (append (subseq l2 0 pos)
                                                     (subseq l2 (+ pos 1)))))))))
                 (eq-sum (cdr e1) (cdr e2))))
          (and (or (isfunc '- e1) (isfunc '- e2))
               (destructuring-bind (me oe) (if (isfunc '- e1)
                                               (list e1 e2)
                                               (list e2 e1))
                 (or (and (isfunc '- oe)  ;; a-b = b-a ; a-b-c-d = (b+c+d)-a
                          (= 3 (length oe))
                          (equal-expr (caddr oe) (cadr me))
                          (or (and (= 3 (length me))
                                   (equal-expr (cadr oe) (caddr me)))
                              (equal-expr (cons '+ (cddr me))
                                          (cadr oe))))
                     (and (isfunc '+ oe) ;; a-b-c-d=b+c+d-a
                          (let ((pos (position-if (lambda (e)
                                                    (equal-expr-1 e (cadr me)))
                                                  oe)))
                            (and pos
                                 (equal-arg-sets (cddr me)
                                                 `(,@(subseq oe 1 pos)
                                                   ,@(subseq oe (+ pos 1))))))))))
          (and (or (isfunc '* e1) (isfunc '* e2) (isfunc '/ e1) (isfunc '/ e2))
               (multiple-value-bind (me oe) (if (or (isfunc '* e1) (isfunc '/ e1))
                                                (values e1 e2)
                                                (values e2 e1))
                 (labels ((get-n-nn (e)
                            (values (apply #'* (remove-if-not #'numberp (cdr e)))
                                    (remove-if #'numberp (cdr e))))
                          (count-pm (l1 l2 &key (nm 0))
                            (if (and l1 l2)
                                (if-let ((pos (position (car l1) l2 :test #'equal-expr)))
                                  (count-pm (cdr l1)
                                            (append (subseq l2 0 pos) (subseq l2 (+ pos 1)))
                                            :nm nm)
                                  (if-let ((pos (position (car l1) l2 :test #'equal-expr-1)))
                                    (count-pm (cdr l1)
                                              (append (subseq l2 0 pos) (subseq l2 (+ pos 1)))
                                              :nm (+ nm 1))
                                    0))
                                (if (or (and l1 (mequal (car l1) -1))
                                        (and l2 (mequal (car l2) -1)))
                                    (+ nm 1)
                                    (if (and (not l1) (not l2))
                                        nm
                                        0)))))
                    (multiple-value-bind (mn mnn) (get-n-nn me)
                      (or (and (or (isfunc '* oe) (isfunc '/ oe))
                               (multiple-value-bind (on onn) (get-n-nn oe)
                                 (oddp (count-pm (cons mn mnn) (cons on onn)))))
                          (and (not (cdr mnn))
                               (or (and (= mn -1)
                                        (equal-expr (car mnn) oe))
                                   (and (= mn 1)
                                        (equal-expr-1 (car mnn) oe))))))))))))

(defun expr-weight (e &key (weight 0) max-weight)
  "Calculate a difficulty to calculate the expression e. If the weight exceeds max-weight, don't calculate more."
  (if (and max-weight (> weight max-weight))
      (return-from expr-weight weight))
  (if (listp e)
      (+ (apply #'+ (mapcar #'expr-weight (cdr e)))
         (* (length (cdr e))
            (case (car e)
                  ((+ -) 1)
                  (abs 2)
                  (cast 3) ;; cast type, like integer -> float, etc
                  (* 4)
                  (exp 50)
                  (/ 30)
                  ((sqrt expt) 100)
                  ((sin cos tan) 60)
                  ((asin acos atan) 80)
                  (otherwise 200))))
      (if (arrayp e)
          (loop for i below (apply #'* (array-dimensions e))
                sum (row-major-aref e i) into s
                when (> (+ s weight) max-weight) return s
                finally (return s))
          0.1)))

; (defmacro defun-stable-expr (name args &rest cod)
;   "Define the function of one argument, which will call itself with its result as a argument, until already a seen result will be returned"
;   (let ((h (gensym))
;         (r (gensym))
;         (arg (gensym)))
;     `(defun ,name (,arg &optional ,h)
;        (let* ((,r (funcall (lambda ,args ,@cod) ,arg)))
;          (if (position ,r ,h :test #'equal-expr)
;              ,r ; (caar (stable-sort (mapcar #'cons ,h (mapcar #'expr-eqight ,h)) #'< :key #'cdr))
;              (,name ,r (cons ,r ,h)))))))

(defun collect-exprs (e)
  "Collecting * and + expressions: (+ a (+ b c)) => (+ a b c)"
  (if (listp e)
      (if (position (car e) '(+ *))
          (cons (car e)
                (mapcan
                  (lambda (e1)
                    (if (listp e1)
                        (if (equal (car e) (car e1))
                            (cdr e1)
                            (list e1))
                        (list e1)))
                  (copy-list (cdr e))))
          e)
      e))

(def-expr-cond extract-nums e ;; Compute everything where it is possible
  :op ((+ *)
       (let* ((op (car e))
              (args (cdr e))
              (num (apply (symbol-function op) (remove-if-not #'numberp args)))
              (not-nums (remove-if #'numberp args)))
         (if not-nums
             (if (equal op '*)
                 (cond ((mequal num 0)
                        0)
                       ((mequal num 1)
                        (if (cdr not-nums)
                            (cons op not-nums)
                            (car not-nums)))
                       (t `(,op ,num ,@not-nums)))
                 (if (mequal num 0)
                     (if (cdr not-nums)
                         (cons op not-nums)
                         (car not-nums))
                     `(,op ,num ,@not-nums)))
             num)))
  :op (expt (if (or (mequal 0 (cadr e))
                    (mequal 1 (cadr e)))
                (cadr e)
                (if (numberp (caddr e))
                    (cond ((= (caddr e) 0) 1)
                          ((= (caddr e) 1) (cadr e))
                          ((numberp (cadr e)) (expt (cadr e) (caddr e)))
                          ((equal '&math-e (cadr e)) (exp (caddr e)))
                          (t e))
                    e)))
  :op (cast (if (numberp (caddr e))
                (cond
                  ((equal 'real*8 (cadr e)) (float (caddr e)))
                  ((equal 'real*4 (cadr e)) (float (caddr e)))
                  ((equal 'integer (cadr e))
                   (multiple-value-bind (f r) (floor (caddr e))
                     (if (= 0.0 r) f e)))
                  (t e))
                e))
  :defop (if (and (fboundp (car e))
                  (every #'numberp (cdr e)))
             (apply (car e) (cdr e))
             e))

(def-expr-cond norm-expr e ;; Convert every / and - to * and +, sqrt to (expt x 1/2)
  :op (- (if (= (length e) 2)
            `(* -1 ,(cadr e))
            `(+ ,(cadr e)
                 (* -1 ,(if (= (length e) 3)
                            (caddr e)
                           `(+ ,@(cddr e)))))))
  :op (/ `(* ,(cadr e) (expt ,(if (= (length e) 3)
                                  (caddr e)
                                 `(* ,@(cddr e)))
                             -1)))
  :op (sqrt `(expt ,(cadr e) 1/2))
  :op (exp `(expt &math-e ,(cadr e)))
  :op (expt (if (and (arrayp (cadr e))
                     (numberp (caddr e))
                     (evenp (caddr e))
                     (equal 1 (array-rank (cadr e))))
                (cond ((equal (caddr e) 0)
                       1)
                      (t `(expt ,(sqr (cadr e)) (/ (caddr e) 2))))
                e)))

(def-expr-cond denorm-expr-expt e ;; The reverse of norm-expr for / and expt
  :op (* (labels ((dvp (x)
                    (and (listp x)
                         (equal (length x) 3)
                         (equal (car x) '/)
                         (mequal (cadr x) 1))))
           (let* ((divs (remove-if-not #'dvp (cdr e)))
                  (ndivs (remove-if #'dvp (cdr e)))
                  (dvop (if divs
                            (if (cdr divs)
                                `(/ 1 (* ,@(mapcar #'caddr divs)))
                                (car divs))))
                  (ndvop (if ndivs
                             (if (cdr ndivs)
                                 `(* ,@ndivs)
                                 (car ndivs)))))
             (if divs
                 (if ndivs
                     `(/ ,ndvop ,(caddr dvop))
                     dvop)
                 e))))
  :op (expt (if (numberp (caddr e))
                (cond ((= (caddr e) 1/2)
                       `(sqrt ,(cadr e)))
                      ((= (caddr e) 1)
                       (cadr e))
                      ((= (caddr e) 0)
                       1)
                      ((< (caddr e) 0)
                       `(/ 1 ,(if (mequal (caddr e) -1)
                                  (cadr e)
                                  `(expt ,(cadr e) ,(* -1 (caddr e))))))

                      (t e))
                (if (equal '&math-e (cadr e))
                    `(exp ,(caddr e))
                    e))))

(def-expr-cond denorm-expr-zop e ;; Replace every (* x) and (+ x) to x
  :op ((+ *) (if (cddr e)
                 e
                 (cadr e))))

(def-expr-cond denorm-expr-minus e ;; Convert (+ ... (* -1 ....) ...)  to (- ...)
  :op (+ (labels ((isneg (x)
                    (and (isfunc '* x)
                         (cddr x)
                         (find-if #'minusp (remove-if-not #'numberp (cdr x)))))
                  (deneg (l)
                    (if l
                        (if (minusp (car l))
                            (if (mequal (car l) -1)
                                (cdr l)
                                (cons (* -1 (car l)) (cdr l)))
                            (cons (car l) (deneg (cdr l)))))))
           (let* ((negs (remove-if-not #'isneg (cdr e)))
                  (poss (remove-if #'isneg (cdr e)))
                  (pnegs (mapcar
                           (lambda (x)
                             (let ((dn (deneg (cdr x))))
                               (if (cdr dn)
                                 `(* ,@dn)
                                 (car dn))))
                           negs)))
             (if negs
                 (if poss
                     `(- ,(if (cdr poss)
                              `(+ ,@poss)
                              (car poss))
                         ,@pnegs)
                     `(* -1 (+ ,@pnegs)))
                 e)))))

(def-expr-cond denorm-expr-minus1 e ;; (+ -5 x) => (- x 5)
  :op (+ (if (and (numberp (cadr e))
                  (< (cadr e) 0))
            `(- ,(if (cdddr e)
                     `(+ ,@(cddr e))
                     (caddr e))
                ,(* -1 (cadr e)))
             e)))

(def-expr-cond denorm-expr-plus-n e ;; Move numbers to the end of (+ ...)
  :op (+ (if (and (numberp (cadr e))
                  (cddr e))
             `(+ ,@(cddr e) ,(cadr e))
              e)))

(defun denorm-expr (e)
  (chain-func-rec
    (denorm-expr-expt
     denorm-expr-plus-n
     denorm-expr-zop
     denorm-expr-minus1
     denorm-expr-minus)
    e))

(def-expr-cond expand-expt1 e ;; (expt (* a ...) n) => (* (expt a n) ...)
  :op (expt (if (isfunc '* (cadr e))
                (cons '* (mapcar
                           (lambda (x)
                             `(expt ,x ,(caddr e)))
                           (cdr (cadr e))))
                e)))

(def-expr-cond expand-expt2 e ;; (expt x (+ n1 ...)) => (* (expt x n) ...)
  :op (expt (if (isfunc '+ (caddr e))
                (cons '* (mapcar
                           (lambda (x)
                             `(expt ,(cadr e) ,x))
                           (cdr (caddr e))))
                e)))

(def-expr-cond collect-expt e ;; Reverse of expand-expt*
  :op (* (if (cddr e)
           (let ((args (list))
                 (cf 1))
             (labels ((process-arg (e)
                        (let ((e1 (if (isfunc 'expt e) (cadr e) e))
                              (d (if (isfunc 'expt e) (caddr e) 1)))
                          (if-let ((r (assoc e1 args :test #'equal-expr)))
                            (setf (cdr r) `(+ ,(cdr r) ,d))
                            (if-let ((r (assoc e1 args :test #'equal-expr-1)))
                              (progn
                                (setf (cdr r) `(+ ,(cdr r) ,d))
                                (setf cf (* cf -1)))
                              (push (cons e1 d) args))))))
               (map nil #'process-arg (cdr e))
               (let ((arg2 `(,@(if (= cf 1) nil (list cf))
                             ,@(mapcar
                                 (lambda (x)
                                   (if (mequal (cdr x) 1) (car x) `(expt ,(car x) ,(cdr x))))
                                 args))))
                 (if (cdr arg2)
                     `(* ,@arg2)
                     (car arg2)))))
           e))
  :op (expt (cond ((isfunc 'expt (cadr e))
                   `(expt ,(cadadr e) (* ,(caddr e) ,(car (cddadr e)))))
                  ((isfunc '* (cadr e))
                   `(* ,@(mapcar (lambda (x)
                                   (if (isfunc 'expt x)
                                      `(expt ,(cadr x) (* ,(caddr e) ,(caddr x)))
                                      `(expt ,x ,(caddr e))))
                                 (cdadr e))))
                  (t e))))

(defparameter *max-degree-expansion* 5)

(def-expr-cond expand-int-expt e ;; (expt x 3) => (* x x x)
  :op (expt (if (and (mintegerp (caddr e))
                     (listp (cadr e))
                     (> (abs (caddr e)) 1)
                     (<= (abs (caddr e)) *max-degree-expansion*))
                (cons '* (loop for i below (abs (caddr e)) collect
                            (if (< (caddr e) 0)
                                `(expt ,(cadr e) -1)
                                (cadr e))))
                e)))

(def-expr-cond expand-mul e ;; (* (+ ...) (+ ...)) => (+ (* ...) (* ...))
  :op (* (let* ((args (cdr e))
                (plus (mapcar #'cdr (mapcar #'collect-exprs (filter-func '+ args))))
                (notplus (filter-func '+ args t)))
           (if plus
               (cons '+
                 (labels ((comb (l1 l2)
                            (loop for x in l1 append
                               (loop for y in l2 collect (cons y x))))
                          (comb-all (l)
                            (if (cdr l)
                                (comb-all (copy-tree (cons (comb (car l) (cadr l)) (cddr l))))
                                (car l))))
                   (loop for x in (comb-all (cons (list nil) plus)) collect
                     (collect-exprs (cons '* (concatenate 'list notplus x))))))
               e))))

(def-expr-cond collect-common-nums e ;; (+ (* 2 x) (* 4 y)) => (* 2 (+ x (* 2 y)))
  :op (+ (let* ((nums (mapcar
                        (lambda (x)
                          (if (numberp x)
                              x
                              (if (and (isfunc '* x)
                                       (numberp (cadr x)))
                                  (cadr x)
                                  1)))
                        (cdr e))))
           (if (not (every #'mintegerp nums))
             e
             (let* ((mcnt (length (remove-if-not #'minusp nums)))
                    (pcnt (- (length nums) mcnt))
                    (g (apply #'gcd (mapcar #'truncate nums)))
                    (g1 (* g (if (> mcnt pcnt) -1 1))))
               (if (> g 1)
                   `(* ,g1 (+ ,@(mapcar
                                  (lambda (x)
                                    (if (numberp x)
                                        (/ x g1)
                                        (if (and (isfunc '* x)
                                                 (numberp (cadr x)))
                                           (if (mequal (cadr x) g1)
                                             (if (cdddr x)
                                               `(* ,@(cddr x))
                                                (caddr x))
                                             `(* ,(/ (cadr x) g1) ,@(cddr x)))
                                           (error "Internal error!!!"))))
                                  (cdr e))))
                   e))))))

(defun is-int (e)
  (or (rationalp e)
      (and (floatp e) (mequal (floor e) e))))

(defun is-int-expt (e)
  (and (isfunc 'expt e)
       (is-int (caddr e))))

(defun eqf (e1 e2 &key int-expt) ;; helper function for extract-subexpr and collect-one-common
  (let ((expt-fn (if int-expt
                     #'is-int-expt
                     (lambda (x) (isfunc 'expt x)))))
    (or (equal-expr e1 e2)
        (equal-expr-1 e1 e2)
        (and (funcall expt-fn e2)
             (or (equal-expr e1 (cadr e2))
                 (equal-expr-1 e1 (cadr e2))))
        (and (funcall expt-fn e1)
             (or (equal-expr e2 (cadr e1))
                 (equal-expr-1 e2 (cadr e1))))
        (and (funcall expt-fn e1)
             (funcall expt-fn e2)
             (or (equal-expr (cadr e1) (cadr e2))
                 (equal-expr-1 (cadr e1) (cadr e2)))))))


(defun extract-subexpr-norm (e subex) ;; convert e (must be normalized) to the form: (subex^n)*e1+e2, return (values n e1 e2)
  (labels
    ((equal-expr-2 (e1 e2)
       (or (equal-expr e1 e2)
           (equal-expr-1 e1 e2)))
     (get-deg (x &key in-mul)
       (cond ((equal-expr subex x)
              (list 1 1))
             ((equal-expr-1 subex x)
              (list 1 -1))
             ((and (isfunc 'expt x)
                   (equal-expr subex (cadr x)))
              (list (caddr x) 1))
             ((and (isfunc 'expt x)
                   (equal-expr-1 (cadr x) subex))
              (list (caddr x) `(expt -1 ,(caddr x))))
             ((and (not in-mul)
                   (isfunc '* x))
              (map nil
                (lambda (x1)
                  (if-let ((res (get-deg x1 :in-mul t)))
                    (return-from get-deg res)))
                (cdr x))))))
    (let* ((args (if (isfunc '+ e)
                     (if (isfunc '+ subex)
                         (labels ((rmse (le ls &optional res m)
                                    (if (and le ls)
                                        (let ((pos (position (car le) ls :test (if m #'equal-expr-1 #'equal-expr))))
                                          (multiple-value-call #'rmse
                                            (cdr le)
                                            (if pos
                                                `(,@(subseq ls 0 pos) ,@(subseq ls (1+ pos)))
                                                ls)
                                            (if pos
                                                res
                                                (cons (car le) res))
                                            m))
                                        (if ls
                                            `(,@(reverse res) ,@le)
                                            (values `(,@(reverse res) ,(if m `(* -1 ,subex) subex) ,@le) t)))))
                            (multiple-value-bind (args fnd) (rmse (cdr e) (cdr subex))
                              (if fnd
                                  args
                                  (rmse (cdr e) (cdr subex) nil t))))
                         (cdr e))
                     (list e)))
           (dgsl (mapcar #'cons
                         (mapcar #'get-deg args)
                         args))
           (ee-args (remove-if #'null dgsl :key #'car))
           (nee-args (mapcar #'cdr (remove-if-not #'null dgsl :key #'car)))
           (degs (remove-duplicates (mapcar #'caar ee-args) :test #'equal-expr))
           (ee-deg (if (and degs (every #'is-int degs))
                       (apply #'min degs)
                       (if (cdr degs)
                           0
                           (if degs
                               (car degs)
                               0))))
           (eargs (mapcar
                    (lambda (de)
                      (let ((dg (caar de))
                            (ml (cadr (car de)))
                            (ex (cdr de)))
                        (cond
                          ((and (mequal dg ee-deg)
                                (equal-expr-2 ex subex))
                           ml)
                          ((and (isfunc 'expt ex)
                                (equal-expr-2 (cadr ex) subex))
                           `(* ,ml (expt ,(cadr ex) (+ ,(caddr ex) (* -1 ,ee-deg)))))
                          ((isfunc '* ex)
                           (labels ((rpe (exl)
                                      (if exl
                                          (let ((ex2 (car exl)))
                                            (if (eqf ex2 subex)
                                              (cond ((mequal dg ee-deg)
                                                     (cons ml (cdr exl)))
                                                    ((isfunc 'expt ex2)
                                                     (cons
                                                       `(* ,ml (expt ,(cadr ex2) (+ ,(caddr ex2) (* -1 ,ee-deg))))
                                                       (cdr exl)))
                                                    (t (cons
                                                         `(* ,ml (expt ,ex2 (+ ,dg (* -1 ,ee-deg))))
                                                         (cdr exl))))
                                              (cons ex2 (rpe (cdr exl)))))
                                          (error (format nil "Internal error: common term not found!")))))
                             (cons '* (rpe (cdr ex)))))
                          (t (error "Internal error: cannot find expr to collect")))))
                    ee-args)))
      (values ee-deg
              (if (cdr eargs)
                  (cons '+ eargs)
                  (if eargs
                      (car eargs)
                      0))
              (if (cdr nee-args)
                  (cons '+ nee-args)
                  (if nee-args
                      (car nee-args)
                      0))))))

(defun extract-subexpr (e subex &key expand) ;; normailze -> extract-subexpr-norm -> denormalize, if expant is T, expand all muls
  (let* ((en (math-rec-funcall #'norm-expr e))
         (en (if expand
                 (math-rec-funcall #'extract-nums
                   (math-rec-funcall #'collect-exprs
                     (math-rec-funcall #'expand-mul
                       (math-rec-funcall #'expand-int-expt
                         (math-rec-funcall #'extract-nums
                           (math-rec-funcall #'collect-exprs
                             (math-rec-funcall #'extract-nums
                               (math-rec-funcall #'collect-expt
                                 (math-rec-funcall #'expand-expt2
                                  (math-rec-funcall #'expand-expt1
                                    (math-rec-funcall #'collect-exprs en)))))))))))
                 en))
         (en (if expand
                 (math-rec-funcall #'extract-nums (cons (car en) (mapcar #'collect-expt (cdr en))))
                 en))
         (sn (math-rec-funcall #'norm-expr subex))
         (sn (if expand
                 (math-rec-funcall #'extract-nums
                   (math-rec-funcall #'collect-exprs
                     (math-rec-funcall #'expand-mul
                       (math-rec-funcall #'expand-int-expt
                         (math-rec-funcall #'extract-nums
                           (math-rec-funcall #'collect-exprs
                             (math-rec-funcall #'extract-nums
                               (math-rec-funcall #'collect-expt
                                 (math-rec-funcall #'expand-expt2
                                  (math-rec-funcall #'expand-expt1
                                    (math-rec-funcall #'collect-exprs sn)))))))))))
                 sn))
         (sn (if (and expand (listp sn))
                 (math-rec-funcall #'extract-nums (cons (car en) (mapcar #'collect-expt (cdr sn))))
                 sn)))
    (multiple-value-bind (n e1 e2) (extract-subexpr-norm en sn)
      (if (mequal n 0)
          (values 0 0 e)
          (labels ((den (e)
                     (funcall (if expand #'simplify #'identity)
                              (chain-func-rec (extract-nums denorm-expr) e))))
            (values (den n)
                    (den e1)
                    (den e2)))))))

(defun get-polynome-cfs (e v &key expand) ;; Convert e to polynome against v and return alist of coefficients, like ((0 . zero-cf) (1 . 1-cf) (2 . 2-cf) ...)
  (multiple-value-bind (n e1 e2) (extract-subexpr e v :expand expand)
    (cond ((> n 0) (cons (cons 0 e2)
                         (mapcar (lambda (cf) (cons (+ (car cf) n) (cdr cf))) (get-polynome-cfs e1 v :expand expand))))
          ((= n 0) (list (cons 0 e2)))
          (t (error (format nil "Strange polynomial coefficient found: ~A" n))))))

(defun collect-one-common (e) ;; extract one most common expression from e. Don't ever try to improve, if not understand it completely.
  (if (isfunc '+ e)
    (let ((e (math-rec-funcall #'collect-exprs e))
          (ecnt-l nil)
          (eqf (lambda (e1 e2) (eqf e1 e2 :int-expt t))))
      (labels
        ((count-expr (e1 nt)
           (if (not (numberp e1))
             (labels ((inc-hash-test (e)
                        (if (not (numberp e))
                            (if-let ((vk (rassoc e ecnt-l :test eqf)))
                                    (pushnew nt (car vk))
                                    (push (cons (list nt) e) ecnt-l)))))
               (map nil
                 (lambda (e2)
                   (if (is-int-expt e2)
                       (inc-hash-test (cadr e2))
                       (inc-hash-test e2)))
                 (if (isfunc '* e1)
                     (cons e1 (cdr e1))
                     (list e1)))))))
        (map nil
             #'count-expr
             (cdr e)
             (loop for i below (length (cdr e)) collect i))
        (let* ((mxc (apply #'max (mapcar (lambda (x) (length (car x))) ecnt-l)))
               (subexs (if (= mxc 1)
                           (mapcar
                             (lambda (x)
                               (cons '+ x))
                             (let ((sum-terms nil))
                               (map nil
                                 (lambda (e)
                                   (when (isfunc '* e)
                                     (map nil (lambda (e)
                                                (when (isfunc '+ e)
                                                  (map nil (lambda (e)
                                                             (pushnew e sum-terms :test (disjoin #'equal-expr #'equal-expr-1)))
                                                           (cdr e))))
                                              (cdr e))))
                                 (cdr e))
                               (let ((terms (remove-if-not (rcurry #'find sum-terms :test (disjoin #'equal-expr #'equal-expr-1)) (cdr e))))
                                 (when (cdr terms)
                                   (all-list-decs terms :min 2 :max (max 2 (length terms))))))))))
          (if (= mxc 1)
              (map nil
                   #'count-expr
                   subexs
                   (loop for i below (length subexs) collect (+ i (length (cdr e))))))
          (if (> (apply #'max (mapcar (lambda (x) (length (car x))) ecnt-l)) 1)
            (let* ((ee (cdr (car (stable-sort (reverse ecnt-l) #'> :key (lambda (x) (length (car x)))))))
                   (ee (if (is-int-expt ee)
                           (cadr ee)
                           ee)))
              (multiple-value-bind (dg e1 e2) (extract-subexpr-norm e ee)
                (if (not (mequal dg 0))
                    (values (let* ((eedg (if (mequal dg 1)
                                             ee
                                             `(expt ,ee ,dg)))
                                   (ee1 (if (mequal e1 1)
                                            eedg
                                            `(* ,eedg ,e1))))
                              (if (mequal e2 0)
                                  ee1
                                  `(+ ,ee1 ,e2)))
                            t)
                    (values e nil))))
            (values e nil)))))
    (values e nil)))

(defun-stable-expr collect-common (e) ;; collect all comvon subexprs
  (labels ((lp (e &optional (deep 0))
             (if (> deep 1000) (error "Too deep recursion in collect-common"))
             (let* ((e (math-rec-funcall #'collect-exprs
                         (math-rec-funcall #'collect-common-nums
                           (math-rec-funcall #'extract-nums e)))))
               (if (listp e)
                   (let ((e (cons (car e) (mapcar (lambda (x) (lp x (1+ deep))) (cdr e)))))
                     (if (isfunc '+ e)
                         (multiple-value-bind (e1 r) (collect-one-common e)
                           (if r
                               (lp e1 (1+ deep))
                               e))
                         e))
                   e))))
    (lp e)))

(def-expr-cond calc-arrays e
   :op ((+ -) (if (some #'arrayp (cdr e))
                  (let* ((zar (make-array (array-dimensions (find-if #'arrayp (cdr e)))
                                          :initial-element 0))
                         (me (mapcar
                               (lambda (x)
                                 (if (mequal x 0)
                                     zar
                                     x))
                               (cdr e))))
                    (if (equal-dimsp me)
                        (apply #'map-array
                          (cons (lambda (&rest els)
                                  `(,(car e) ,@els))
                                me))
                        (error "Cannot mix elements with different ranks in +/-")))
                  e))
   :op (* (let ((arrs (remove-if-not #'arrayp (cdr e)))
                (not-arrs (remove-if #'arrayp (cdr e))))
            (if arrs
                (reduce
                  #'array-multiply
                  (cdr arrs)
                  :initial-value (if not-arrs
                                     (array-multiply (if (cdr not-arrs)
                                                         `(* ,@not-arrs)
                                                          (car not-arrs))
                                                     (car arrs))
                                     1))
                e)))
   :op (/ (calc-arrays `(* ,(cadr e) (expt ,(if (cdddr e)
                                                `(* ,@(cddr e))
                                                (caddr e))
                                           -1))))
   :op (expt (if (arrayp (caddr e))
                 (error "Array cannot be an exponent, sorry")
                 (if (arrayp (cadr e))
                     (let ((ar (cadr e))
                           (ex (caddr e)))
                       (if (and (numberp ex)
                                (integerp ex)
                                (>= ex 0))
                           (cond ((mequal ex 0) 1)
                                 ((mequal ex 1) ar)
                                 (t (if (cdr (array-dimensions ar))
                                        (reduce
                                          #'array-multiply
                                          (loop repeat (- ex 1) collect ar)
                                          :initial-value ar)
                                        (if (evenp ex)
                                            `(+ ,@(loop for i across ar collect `(expt ,i ,ex)))
                                            (array-multiply `(+ ,@(loop for i across ar collect `(expt ,i ,(- ex 1)))) ar)))))
                           (error "(expt array e): e must be integer >= 0")))
                     e)))
   :op (vector (make-array (list (length (cdr e)))
                 :initial-contents (cdr e)))
   :op (aref (apply #'aref (cdr e)))
   :op (sqr (sqr (cadr e))))

(def-expr-cond collect-same-expts e ;; (* (expt x 2) (expt y 2)) => (expt (* x y) 2)
   :op (* (labels ((is-expt (x) (isfunc 'expt x))
                   (is-the-expt (ex) (lambda (e) (equal-expr (caddr e) ex))))
            (let* ((expts (remove-if-not #'is-expt (cdr e)))
                   (non-expts (remove-if #'is-expt (cdr e)))
                   (expt-exs (remove-duplicates (mapcar #'caddr expts) :test #'equal-expr)))
              (if (and (cdr expts)
                       (< (length expt-exs) (length expts)))
                  (let ((expt-terms (mapcar
                                      (lambda (ex)
                                        (let ((expts-ex (remove-if-not (is-the-expt ex) expts)))
                                          `(expt ,(if (cdr expts-ex)
                                                      `(* ,@(mapcar #'cadr expts-ex))
                                                      (cadr (car expts-ex)))
                                                 ,ex)))
                                      expt-exs)))
                    (if non-expts
                        `(* ,@non-expts ,@expt-terms)
                        (if (cdr expt-terms)
                            `(* ,@expt-terms)
                            (car expt-terms))))
                  e)))))

(defmacro chain-debug (expr)
  (if (listp expr)
      (let* ((ll (cddr expr))
             (fn (if ll (cadr expr) (car expr)))
             (op (if ll (caddr expr) (cadr expr)))
             (res (gensym)))
        `(let ((,res (chain-debug ,op)))
           (format t "~&(~A ~A)~%" ',fn ,res)
           (,@(if ll '(math-rec-funcall)) ,fn ,res)))
      expr))

(defun-stable-expr simplify-expr3 (e) ;; Simplify normalized expression
  ;;(chain-debug
  (math-rec-funcall #'extract-nums
    (math-rec-funcall #'collect-same-expts
      (math-rec-funcall #'collect-expt
        (collect-common
          (math-rec-funcall #'extract-nums
            (math-rec-funcall #'expand-mul
              (math-rec-funcall #'expand-int-expt
                (math-rec-funcall #'extract-nums
                  (math-rec-funcall #'collect-exprs
                    (math-rec-funcall #'extract-nums
                      (math-rec-funcall #'collect-expt
                        (math-rec-funcall #'expand-expt2
                         (math-rec-funcall #'expand-expt1
                           (math-rec-funcall #'collect-exprs e)))))))))))))))

(defparameter *simplify-cache* (list))

(defun clear-simplify-cache ()
  (setf *simplify-cache* (list)))

(defun simplify-expr2 (e) ;; normalize -> simplify -> denormalize
  (if-let ((r (assoc e *simplify-cache* :test #'equal-expr)))
    (cdr r)
    (let ((res (math-rec-funcall #'denorm-expr
                 (simplify-expr3 (math-rec-funcall #'norm-expr e)))))
      (push (cons e res) *simplify-cache*)
      res)))

(defun simplify (e)
  (let ((e1 (math-rec-funcall #'calc-arrays e)))
    (if (arrayp e1)
        (map-array #'simplify-expr2 e1)
        (simplify-expr2 e1))))

(defun count-subexprs (e &key (hash (make-hash-table)))
  (labels ((asc (e)
             (let ((v (gethash e hash)))
               (if v
                   (values e v)
                   (maphash (lambda (k v)
                              (when (equal-expr e k)
                                (return-from asc (values k v))))
                            hash))
               nil))
           (find-se (e &optional no-variants)
             (when (listp e)
               (multiple-value-bind (k v) (asc e)
                 (if k
                     (setf (gethash k hash) (1+ v))
                     (progn
                       (setf (gethash e hash) 1)
                       (when (not no-variants)
                         (if (or (equal (car e) '*)
                                 (equal (car e) '+))
                             (loop for se in (symath::all-list-decs (cdr e) :max (1- (length (cdr e))))
                               if (cdr se) do (find-se (cons (car e) se) t)
                                  else do (find-se (car se))
                               (cdr e))
                             (map nil #'find-se (cdr e))))))))))
    (find-se e)
    hash))

(defun replace-subexpr (e e1 e2)
  (let ((cnt 0))
    (labels ((rep-e-args (args &optional es (rem-e (cdr e1)))
               (if (not rem-e)
                   (progn
                     (incf cnt)
                     (cons e2 (rep-e-args args)))
                   (if args
                       (let ((pos (position (car args) rem-e :test #'equal-expr)))
                         (if pos
                             (rep-e-args (cdr args) (cons (car args) es) (append (subseq rem-e 0 pos) (subseq rem-e (1+ pos))))
                             (cons (car args) (rep-e-args (cdr args) es rem-e))))
                       es)))
             (rep-e (e)
               (if (symath::equal-expr e e1)
                   (progn
                     (incf cnt)
                     e2)
                   (if (listp e)
                       (cons (car e)
                             (if (and (listp e1)
                                      (equal (car e) (car e1))
                                      (or (equal '+ (car e))
                                          (equal '* (car e))))
                                 (rep-e-args (mapcar #'rep-e (cdr e)))
                                 (mapcar #'rep-e (cdr e))))
                       e))))
      (values (rep-e e) cnt))))

(defvar *tmp-cnt* 0)
(defun gen-tmp-var (&rest args)
  (declare (ignore args))
  (intern (format nil "tmp~A" (incf *tmp-cnt*))))

(defmacro with-var-cnt-reset (&rest code)
  `(let ((*tmp-cnt* 0))
     ,@code))

(defun split-to-subexprs (vcs &key temps (min-weight 0) (gen-tmp #'gen-tmp-var) subst-self)
  (let ((orig-vrs (mapcar #'car vcs)))
    (labels ((apply-reps (e reps)
                (if reps
                    (apply-reps (replace-subexpr e (cdar reps) (caar reps))
                                (cdr reps))
                    e))
             (apply-reps-ve (reps)
               (lambda (ve)
                 (cons (car ve) (apply-reps (cdr ve) reps))))
             (simp (ve)
               (cons (car ve) (simplify (cdr ve))))
             (rep-exs (el rpl tmps)
               (if rpl
                   (let ((te (apply-reps (car rpl) tmps)))
                     (if (and (listp te) (or (= 0 min-weight) (> (expr-weight te) min-weight)))
                         (let* ((tve (cons (funcall gen-tmp te) te))
                                (ltve (list tve)))
                           (multiple-value-call #'rep-exs
                             (mapcar (apply-reps-ve ltve) el)
                             (cdr rpl)
                             (cons tve tmps)))
                         (multiple-value-call #'rep-exs el (cdr rpl) tmps)))
                   (values (append (mapcar #'simp (reverse tmps))
                                   (mapcar #'simp el))
                           (reverse (mapcar #'car tmps)))))
             (orig-dep (e)
               (if (position e orig-vrs :test #'equal-expr)
                   (return-from orig-dep t)
                   (if (listp e)
                       (map nil (lambda (e)
                                  (when (orig-dep e)
                                    (return-from orig-dep t)))
                                (cdr e))))))
      (let* ((hs (make-hash-table))
             (vcs (mapcar (apply-reps-ve temps) vcs))
             (vcs (if subst-self
                      (labels ((rep (vl &optional rl)
                                 (if vl
                                     (cons (cons (caar vl) (apply-reps (cdar vl) rl))
                                           (rep (cdr vl) (cons (car vl) rl))))))
                        (rep vcs))
                      vcs)))
        (mapcar (lambda (e) (count-subexprs (cdr e) :hash hs)) vcs)
        (multiple-value-call #'rep-exs
          vcs
          (mapcar #'car (sort (loop for k being the hash-key of hs and v being the hash-value of hs
                                     when (and (> v 1)
                                               (or (= 0 min-weight) (> (expr-weight k) min-weight))
                                               (not (orig-dep k)))
                                     collect (cons k v))
                              #'>
                              :key #'cdr))
          (reverse temps))))))

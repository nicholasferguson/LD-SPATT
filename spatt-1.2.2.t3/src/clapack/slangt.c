#include "blaswrap.h"
#include "f2c.h"

doublereal slangt_(char *norm, integer *n, real *dl, real *d__, real *du)
{
/*  -- LAPACK auxiliary routine (version 3.1) --   
       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
       November 2006   


    Purpose   
    =======   

    SLANGT  returns the value of the one norm,  or the Frobenius norm, or   
    the  infinity norm,  or the  element of  largest absolute value  of a   
    real tridiagonal matrix A.   

    Description   
    ===========   

    SLANGT returns the value   

       SLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum),   
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
    normF  denotes the  Frobenius norm of a matrix (square root of sum of   
    squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in SLANGT as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, SLANGT is   
            set to zero.   

    DL      (input) REAL array, dimension (N-1)   
            The (n-1) sub-diagonal elements of A.   

    D       (input) REAL array, dimension (N)   
            The diagonal elements of A.   

    DU      (input) REAL array, dimension (N-1)   
            The (n-1) super-diagonal elements of A.   

    =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2, r__3, r__4, r__5;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer i__;
    static real sum, scale;
    extern logical lsame_(char *, char *);
    static real anorm;
    extern /* Subroutine */ int slassq_(integer *, real *, integer *, real *, 
	    real *);


    --du;
    --d__;
    --dl;

    /* Function Body */
    if (*n <= 0) {
	anorm = 0.f;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	anorm = (r__1 = d__[*n], dabs(r__1));
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    r__2 = anorm, r__3 = (r__1 = dl[i__], dabs(r__1));
	    anorm = dmax(r__2,r__3);
/* Computing MAX */
	    r__2 = anorm, r__3 = (r__1 = d__[i__], dabs(r__1));
	    anorm = dmax(r__2,r__3);
/* Computing MAX */
	    r__2 = anorm, r__3 = (r__1 = du[i__], dabs(r__1));
	    anorm = dmax(r__2,r__3);
/* L10: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = dabs(d__[1]);
	} else {
/* Computing MAX */
	    r__3 = dabs(d__[1]) + dabs(dl[1]), r__4 = (r__1 = d__[*n], dabs(
		    r__1)) + (r__2 = du[*n - 1], dabs(r__2));
	    anorm = dmax(r__3,r__4);
	    i__1 = *n - 1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		r__4 = anorm, r__5 = (r__1 = d__[i__], dabs(r__1)) + (r__2 = 
			dl[i__], dabs(r__2)) + (r__3 = du[i__ - 1], dabs(r__3) );
		anorm = dmax(r__4,r__5);
/* L20: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	if (*n == 1) {
	    anorm = dabs(d__[1]);
	} else {
/* Computing MAX */
	    r__3 = dabs(d__[1]) + dabs(du[1]), r__4 = (r__1 = d__[*n], dabs(
		    r__1)) + (r__2 = dl[*n - 1], dabs(r__2));
	    anorm = dmax(r__3,r__4);
	    i__1 = *n - 1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		r__4 = anorm, r__5 = (r__1 = d__[i__], dabs(r__1)) + (r__2 = 
			du[i__], dabs(r__2)) + (r__3 = dl[i__ - 1], dabs(r__3) );
		anorm = dmax(r__4,r__5);
/* L30: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.f;
	sum = 1.f;
	slassq_(n, &d__[1], &c__1, &scale, &sum);
	if (*n > 1) {
	    i__1 = *n - 1;
	    slassq_(&i__1, &dl[1], &c__1, &scale, &sum);
	    i__1 = *n - 1;
	    slassq_(&i__1, &du[1], &c__1, &scale, &sum);
	}
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of SLANGT */

} /* slangt_ */

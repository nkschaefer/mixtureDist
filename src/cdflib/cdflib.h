/**
 * cdflib made available by the authors:
 *
 * Barry Brown, James Lovato, Kathy Russell,
 * Department of Biomathematics,
 * University of Texas,
 * Houston, Texas. 
 * 
 * At this website: https://people.sc.fsu.edu/~jburkardt/c_src/cdflib/cdflib.html
 * 
 * Under the GNU LGPL (below):
 * 
 * ===============
 * 
 * GNU LESSER GENERAL PUBLIC LICENSE
 * 
 * Version 3, 29 June 2007
 * 
 * Copyright © 2007 Free Software Foundation, Inc. <https://fsf.org/>
 * 
 * Everyone is permitted to copy and distribute verbatim copies of this license document, 
 * but changing it is not allowed.
 * 
 * This version of the GNU Lesser General Public License incorporates the terms and 
 * conditions of version 3 of the GNU General Public License, supplemented by the 
 * additional permissions listed below.
 * 0. Additional Definitions.
 * 
 * As used herein, "this License" refers to version 3 of the GNU Lesser General Public 
 * License, and the "GNU GPL" refers to version 3 of the GNU General Public License.
 * 
 * "The Library" refers to a covered work governed by this License, other than an 
 * Application or a Combined Work as defined below.
 * 
 * An "Application" is any work that makes use of an interface provided by the Library, 
 * but which is not otherwise based on the Library. Defining a subclass of a class 
 * defined by the Library is deemed a mode of using an interface provided by the Library.
 * 
 * A "Combined Work" is a work produced by combining or linking an Application with 
 * the Library. The particular version of the Library with which the Combined Work was 
 * made is also called the "Linked Version".
 * 
 * The "Minimal Corresponding Source" for a Combined Work means the Corresponding 
 * Source for the Combined Work, excluding any source code for portions of the Combined 
 * Work that, considered in isolation, are based on the Application, and not on the 
 * Linked Version.
 * 
 * The "Corresponding Application Code" for a Combined Work means the object code 
 * and/or source code for the Application, including any data and utility programs 
 * needed for reproducing the Combined Work from the Application, but excluding the 
 * System Libraries of the Combined Work.
 * 1. Exception to Section 3 of the GNU GPL.
 * 
 * You may convey a covered work under sections 3 and 4 of this License without being 
 * bound by section 3 of the GNU GPL.
 * 2. Conveying Modified Versions.
 * 
 * If you modify a copy of the Library, and, in your modifications, a facility refers 
 * to a function or data to be supplied by an Application that uses the facility 
 * (other than as an argument passed when the facility is invoked), then you may convey 
 * a copy of the modified version:
 * 
 *   a) under this License, provided that you make a good faith effort to ensure that, 
 *   in the event an Application does not supply the function or data, the facility 
 *   still operates, and performs whatever part of its purpose remains meaningful, or
 *    b) under the GNU GPL, with none of the additional permissions of this License 
 *    applicable to that copy.
 *
 * 3. Object Code Incorporating Material from Library Header Files.
 * 
 * The object code form of an Application may incorporate material from a header file 
 * that is part of the Library. You may convey such object code under terms of your 
 * choice, provided that, if the incorporated material is not limited to numerical 
 * parameters, data structure layouts and accessors, or small macros, inline functions 
 * and templates (ten or fewer lines in length), you do both of the following:
 * 
 *   a) Give prominent notice with each copy of the object code that the Library is 
 *   used in it and that the Library and its use are covered by this License.
 *   b) Accompany the object code with a copy of the GNU GPL and this license document.
 *
 * 4. Combined Works.
 * 
 * You may convey a Combined Work under terms of your choice that, taken together, 
 * effectively do not restrict modification of the portions of the Library contained 
 * in the Combined Work and reverse engineering for debugging such modifications, if 
 * you also do each of the following:
 * 
 *   a) Give prominent notice with each copy of the Combined Work that the Library is 
 *   used in it and that the Library and its use are covered by this License.
 *   b) Accompany the Combined Work with a copy of the GNU GPL and this license 
 *   document.
 *   c) For a Combined Work that displays copyright notices during execution, include 
 *   the copyright notice for the Library among these notices, as well as a reference 
 *   directing the user to the copies of the GNU GPL and this license document.
 *   d) Do one of the following:
 *       0) Convey the Minimal Corresponding Source under the terms of this License, 
 *       and the Corresponding Application Code in a form suitable for, and under 
 *       terms that permit, the user to recombine or relink the Application with a 
 *       modified version of the Linked Version to produce a modified Combined Work, 
 *       in the manner specified by section 6 of the GNU GPL for conveying 
 *       Corresponding Source.
 *       1) Use a suitable shared library mechanism for linking with the Library. A 
 *       suitable mechanism is one that (a) uses at run time a copy of the Library 
 *       already present on the user's computer system, and (b) will operate properly 
 *       with a modified version of the Library that is interface-compatible with the 
 *       Linked Version.
 *   e) Provide Installation Information, but only if you would otherwise be required 
 *   to provide such information under section 6 of the GNU GPL, and only to the 
 *   extent that such information is necessary to install and execute a modified 
 *   version of the Combined Work produced by recombining or relinking the Application 
 *   with a modified version of the Linked Version. (If you use option 4d0, the 
 *   Installation Information must accompany the Minimal Corresponding Source and 
 *   Corresponding Application Code. If you use option 4d1, you must provide the 
 *   Installation Information in the manner specified by section 6 of the GNU GPL 
 *   for conveying Corresponding Source.)
 *
 * 5. Combined Libraries.
 * 
 * You may place library facilities that are a work based on the Library side by side 
 * in a single library together with other library facilities that are not Applications 
 * and are not covered by this License, and convey such a combined library under terms 
 * of your choice, if you do both of the following:
 * 
 *   a) Accompany the combined library with a copy of the same work based on the Library, 
 *   uncombined with any other library facilities, conveyed under the terms of this License.
 *   b) Give prominent notice with the combined library that part of it is a work based 
 *   on the Library, and explaining where to find the accompanying uncombined form of the same work.
 *
 * 6. Revised Versions of the GNU Lesser General Public License.
 * 
 * The Free Software Foundation may publish revised and/or new versions of the GNU Lesser 
 * General Public License from time to time. Such new versions will be similar in spirit 
 * to the present version, but may differ in detail to address new problems or concerns.
 * 
 * Each version is given a distinguishing version number. If the Library as you received 
 * it specifies that a certain numbered version of the GNU Lesser General Public License 
 * "or any later version" applies to it, you have the option of following the terms and 
 * conditions either of that published version or of any later version published by the 
 * Free Software Foundation. If the Library as you received it does not specify a version 
 * number of the GNU Lesser General Public License, you may choose any version of the GNU 
 * Lesser General Public License ever published by the Free Software Foundation.
 * 
 * If the Library as you received it specifies that a proxy can decide whether future 
 * versions of the GNU Lesser General Public License shall apply, that proxy's public 
 * statement of acceptance of any version is permanent authorization for you to choose 
 * that version for the Library.
 */
#ifndef CDFLIB_H
#define CDFLIB_H
#ifdef __cplusplus
extern "C"{
#endif

double algdiv ( double *a, double *b );
double alnrel ( double *a );
double apser ( double *a, double *b, double *x, double *eps );
double bcorr ( double *a0, double *b0 );
double beta ( double a, double b );
double beta_asym ( double *a, double *b, double *lambda, double *eps );
double beta_frac ( double *a, double *b, double *x, double *y, double *lambda,
  double *eps );
void beta_grat ( double *a, double *b, double *x, double *y, double *w,
  double *eps,int *ierr );
void beta_inc ( double *a, double *b, double *x, double *y, double *w,
  double *w1, int *ierr );
void beta_inc_values ( int *n_data, double *a, double *b, double *x, double *fx );
double beta_log ( double *a0, double *b0 );
double beta_pser ( double *a, double *b, double *x, double *eps );
double beta_rcomp ( double *a, double *b, double *x, double *y );
double beta_rcomp1 ( int *mu, double *a, double *b, double *x, double *y );
double beta_up ( double *a, double *b, double *x, double *y, int *n, double *eps );
void binomial_cdf_values ( int *n_data, int *a, double *b, int *x, double *fx );
void cdfbet ( int *which, double *p, double *q, double *x, double *y,
  double *a, double *b, int *status, double *bound );
void cdfbin ( int *which, double *p, double *q, double *s, double *xn,
  double *pr, double *ompr, int *status, double *bound );
void cdfchi ( int *which, double *p, double *q, double *x, double *df,
  int *status, double *bound );
void cdfchn ( int *which, double *p, double *q, double *x, double *df,
  double *pnonc, int *status, double *bound );
void cdff ( int *which, double *p, double *q, double *f, double *dfn,
  double *dfd, int *status, double *bound );
void cdffnc ( int *which, double *p, double *q, double *f, double *dfn,
  double *dfd, double *phonc, int *status, double *bound );
void cdfgam ( int *which, double *p, double *q, double *x, double *shape,
  double *scale, int *status, double *bound );
void cdfnbn ( int *which, double *p, double *q, double *s, double *xn,
  double *pr, double *ompr, int *status, double *bound );
void cdfnor ( int *which, double *p, double *q, double *x, double *mean,
  double *sd, int *status, double *bound );
void cdfpoi ( int *which, double *p, double *q, double *s, double *xlam,
  int *status, double *bound );
void cdft ( int *which, double *p, double *q, double *t, double *df,
  int *status, double *bound );
void chi_noncentral_cdf_values ( int *n_data, double *x, double *lambda, 
  int *df, double *cdf );
void chi_square_cdf_values ( int *n_data, int *a, double *x, double *fx );
void cumbet ( double *x, double *y, double *a, double *b, double *cum,
  double *ccum );
void cumbin ( double *s, double *xn, double *pr, double *ompr,
  double *cum, double *ccum );
void cumchi ( double *x, double *df, double *cum, double *ccum );
void cumchn ( double *x, double *df, double *pnonc, double *cum,
  double *ccum );
void cumf ( double *f, double *dfn, double *dfd, double *cum, double *ccum );
void cumfnc ( double *f, double *dfn, double *dfd, double *pnonc,
  double *cum, double *ccum );
void cumgam ( double *x, double *a, double *cum, double *ccum );
void cumnbn ( double *s, double *xn, double *pr, double *ompr,
  double *cum, double *ccum );
void cumnor ( double *arg, double *result, double *ccum );
void cumpoi ( double *s, double *xlam, double *cum, double *ccum );
void cumt ( double *t, double *df, double *cum, double *ccum );
double dbetrm ( double *a, double *b );
double dexpm1 ( double *x );
double dinvnr ( double *p, double *q );
void dinvr ( int *status, double *x, double *fx,
  unsigned long *qleft, unsigned long *qhi );
double dlanor ( double *x );
double dpmpar ( int *i );
void dstinv ( double *zsmall, double *zbig, double *zabsst,
  double *zrelst, double *zstpmu, double *zabsto, double *zrelto );
double dstrem ( double *z );
void dstzr ( double *zxlo, double *zxhi, double *zabstl, double *zreltl );
double dt1 ( double *p, double *q, double *df );
void dzror ( int *status, double *x, double *fx, double *xlo,
  double *xhi, unsigned long *qleft, unsigned long *qhi );
void erf_values ( int *n_data, double *x, double *fx );
double error_f ( double *x );
double error_fc ( int *ind, double *x );
double esum ( int *mu, double *x );
double eval_pol ( double a[], int *n, double *x );
double exparg ( int *l );
void f_cdf_values ( int *n_data, int *a, int *b, double *x, double *fx );
void f_noncentral_cdf_values ( int *n_data, int *a, int *b, double *lambda, 
  double *x, double *fx );
double fifdint ( double a );
double fifdmax1 ( double a, double b );
double fifdmin1 ( double a, double b );
double fifdsign ( double mag, double sign );
long fifidint ( double a );
long fifmod ( long a, long b );
double fpser ( double *a, double *b, double *x, double *eps );
void ftnstop ( char *msg );
double gam1 ( double *a );
void gamma_inc ( double *a, double *x, double *ans, double *qans, int *ind );
void gamma_inc_inv ( double *a, double *x, double *x0, double *p, double *q,
  int *ierr );
void gamma_inc_values ( int *n_data, double *a, double *x, double *fx );
double gamma_ln1 ( double *a );
double gamma_log ( double *a );
void gamma_rat1 ( double *a, double *x, double *r, double *p, double *q,
  double *eps );
void gamma_values ( int *n_data, double *x, double *fx );
double gamma_x ( double *a );
double gsumln ( double *a, double *b );
int ipmpar ( int *i );
void negative_binomial_cdf_values ( int *n_data, int *f, int *s, double *p, 
  double *cdf );
void normal_cdf_values ( int *n_data, double *x, double *fx );
void poisson_cdf_values ( int *n_data, double *a, int *x, double *fx );
double psi ( double *xx );
void psi_values ( int *n_data, double *x, double *fx );
double rcomp ( double *a, double *x );
double rexp ( double *x );
double rlog ( double *x );
double rlog1 ( double *x );
void student_cdf_values ( int *n_data, int *a, double *x, double *fx );
double stvaln ( double *p );
void timestamp ( void );

#ifdef __cplusplus
}
#endif
#endif

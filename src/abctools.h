
#ifndef ___ABCTOOLS_H___
#define ___ABCTOOLS_H___


void mycpyd(double *a, int *len, double *b);
void mycpyi(int *a, int *len, int *b);
void myrevd(double *dx, int *n, double *dy);
void myrevi(int *a, int *la, int *b);
void mysortd(double *a, int *la, double *sorted, int *order, int *inc);
void nnk(double *x, int *lx, int *cols, int *k, double *D);
void nnone(double *x, int *lx, int *k, double *D);
void distanceij(double *x, int *i, int *j, int *p, double *d);
double norm2(double *x, int *n);


#endif

#include <R.h> 
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <time.h>

#include "abctools.h"


void mycpyd(double *a, int *len, double *b){

int i;

for(i=0;i<*len;i++){
*(b+i)=*(a+i);
}

}

/* ********** */

void mycpyi(int *a, int *len, int *b){

int i;

for(i=0;i<*len;i++){
*(b+i)=*(a+i);
}

}

/* ********** */

void myrevd(double *dx, int *n, double *dy){

int incx=-1, incy=1;

F77_CALL(dcopy)(n,dx,&incx,dy,&incy);

}

/* ********** */

void myrevi(int *a, int *la, int *b){

int *tmp=a+(*la-1);  /* *tmp points to the last element of "a" */
int i;

for(i=0;i<*la;i++){
    *b++=*tmp--;
}

}

/* ********** */

void mysortd(double *a, int *la, double *sorted, int *order, int *inc){

int i,*o=calloc(*la,sizeof(int));
double *s=calloc(*la,sizeof(double));



for(i=0;i<*la;i++){
    *(o+i)=(i+1);
}

mycpyd(a,la,s);
rsort_with_index(s,o,*la);

if(*inc==0){
    myrevi(o,la,order);
    myrevd(s,la,sorted);
}
else{
    mycpyi(o,la,order);
    mycpyd(s,la,sorted);
}
free(o);
free(s);
}

/****************************************/

void nnk(double *x, int *lx, int *cols, int *k, double *D){

int i,j, l,m, last=*k-1;

/* MAN: unused variables one, zero 24/6/14
one=1, zero=0, */
int en=*lx, p=*cols; 

double *nearest=calloc(*k,sizeof(double));
double tot=0;
double newval=0;
double xmax=p*5000;		/* something ridiculous for initialization */
double halfp=(double) p; 
	halfp*=0.5;
double *d=calloc(en,sizeof(double));


for(i=0;i<en;i++){
	/* initialize */
	
	for(j=0;j<*k;j++){
		*(nearest+j)=xmax;
	}

	/* work out the distance from x_i to x_j (in R^cols) */
	/* missing out the distance for x_i */
	for(j=0;j<en;j++){
		if(i==j){
			continue;	/* avoid distance computation of i with itself */
		}	
		distanceij(x,&i,&j,&p,&newval);
	
		if (newval < *(nearest+last)){		/* insert val into nearest vector, in one of k places */
			for (l = 0; l < *k; l++){
		    		if (newval < *(nearest+l)) {
		    		/* shift bigger values up one space */
					for (m = last; m > l; m--) {
					    *(nearest+m) = *(nearest+m - 1);
					}
				/* and insert distance in correct place */	
				*(nearest+l) = newval;
				/* we can stop doing the checking now, so 'break' out of for loop  */
				break;			
				}
			}
		}	
	}
	*(d+i)=*(nearest+last);		/* kth nearest neighbour distance */
        tot+= (*(d+i)<=0) ? 0 : log(*(d+i));
}

free(d);
free(nearest);
tot=p*tot/en;

*D=tot +log(R_pow(M_PI,halfp)/gammafn(halfp + 1)) - digamma(*k) +log(en); 

}

/********************************/

void nnone(double *x, int *lx, int *k, double *D){

/* this is a fast nnk for 1 dimension */

int i,j, one=1; 
/* MAN 24/6/14 unused variable lxlow
lxlow=*lx-1;	*/


double *dvec,*sn;
double *sx=calloc(*lx,sizeof(double));
int *o=calloc(*lx,sizeof(int));
int *on;
int en=*lx, left =0, right =0;
double tot=0;
int nn;
double *d=calloc(en,sizeof(double));



/* for each element in x, we only need to check the k neighbours each side, like
in getnbrs.  Also, since we only want the distances and not the associated positions,
we can just sort x first, and then work out distances (one dimension only). */

mysortd(x,lx,sx,o,&one);
free(o);

for(i=0;i<*lx;i++){
	left=0;
	right=0;
	for(j=0;j<*k;j++){
		left+= ((i-j-1)>=0) ? 1 : 0;
	}
	for(j=0;j<*k;j++){
		right+= ((i+j+1)<= (en-1)) ? 1 : 0;
	}
	nn=left+right;
	dvec=calloc(nn,sizeof(double));
        for(j=0;j<left;j++){
                *(dvec+j)=fabs(*(sx+i)-*(sx+i-j-1));
        }
        for(j=0;j<right;j++){
                *(dvec+left+j)=fabs(*(sx+i)-*(sx+i+j+1));
        }
       	sn=calloc(nn,sizeof(double));
	on=calloc(nn,sizeof(int));
        mysortd(dvec,&nn,sn,on,&one);
        *(d+i)=*(sn+*k-1);
        tot+= (*(d+i)<=0) ? 0 : log(*(d+i));
        free(dvec);
        free(sn);
        free(on);
}

free(d);
free(sx);
tot=tot/en;

*D=tot +log(2) - digamma(*k) +log(en); 

}

/***********************************/
/* temporary distance function */

void distanceij(double *x, int *i, int *j, int *p, double *d){

double *v=calloc(*p,sizeof(double));


int k;

/* p-vector of differences */

for(k=0;k<*p;k++){
	*(v+k)=*(x+(*i**p)+k)-*(x+(*j**p)+k);
}

*d=norm2(v,p);
free(v); 

}

/***********************************/

double norm2(double *x, int *n){

int one=1;

/* DP: add return here v0.3-2 */
return F77_CALL(dnrm2)(n,x,&one);

}

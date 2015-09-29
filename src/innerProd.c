#include <Rmath.h>
#include <R.h>
#include <math.h>

/*C code for calculating:
a CUSUM statistic (Inner Products - IP) given y;
finding the argmax of IP;
combining CUSUMs across scales;
repeating for random intervals
*/
void finner_prod_maxp(double *out, int *n, int *max_b, double *p, double *max_inner){
	int iter,s0,e0;
	double aux,maxipi;
	max_inner[0]=0;
	max_b[0] = 0;
	s0 = floor((1-(*p))* (*n)) ;
	e0 = ceil((*p) * (*n));
	maxipi = 0;
	for(iter=s0-1; iter < e0-1; iter++){
		aux = fabs(out[iter]);
		if(aux>maxipi){
			max_b[0] = iter+1;
			maxipi = aux;
		}
	}
	max_inner[0] = out[max_b[0]-1];
}

void fast_inner_prod(double *x, int *n, double *out, double *I_plus, double *I_minus) {
     double sum_of_x, i_a,I_plusinv,fct,npw2,n_inv;
     int iter;
     double n_a = (double)(*n);
     i_a = (double)0;
     I_plusinv=0;
     fct=0;
     n_inv = 1 / n_a;
     npw2 = n_a * n_a;
     sum_of_x = 0;
        for (iter=1; iter < *n; iter++) sum_of_x += x[iter];

        I_minus[0] = 1/sqrt(npw2 - n_a) * sum_of_x;
        I_plus[0] = sqrt(1 - n_inv) * x[0];
        out[0] = I_plus[0] - I_minus[0];

        for (iter=1; iter < (*n) - 1; iter++) {
            i_a = (double)iter;
            I_plusinv =  1/(i_a+1);
            fct = sqrt((n_a-i_a-1) * i_a * I_plusinv / ( n_a-i_a) );
            I_plus[iter] = x[iter] * sqrt(I_plusinv - n_inv) + I_plus[iter-1] * fct;
            I_minus[iter] = I_minus[iter-1] / fct - x[iter] / sqrt(npw2 * I_plusinv - n_a);
            out[iter] = (I_plus[iter] - I_minus[iter]);
        }

        for (iter=0; iter < (*n)-1 ; iter++) out[iter] = out[iter] / pow(sum_of_x/ n_a,1.0);//divide by the mean
}

void across_fip(double *x_a, int *n,  int *d, double *res_a,double *tau,int *epp, double *aux_res, double *aux_res2,
                double *I_plus, double *I_minus, double *Ts,int *p1,int *max_b, double *p, double *max_inner) {
 int i, j, ilow,n_2, scales_temp;
 double temp, crit;
 int TF;
 n_2 = *n;
 ilow=0;
 scales_temp = *d;

 for (j=0; j < scales_temp; j++){
        ilow = j * n_2;
        for (i=0; i < n_2; i++) aux_res[i] = x_a[i+ilow];
    fast_inner_prod(aux_res,n,aux_res2,I_plus,I_minus);

      for (i=0; i < n_2-1; i++){
            TF = ((*p1)<epp[j]);
            if (TF==1) {
             res_a[i]=res_a[i]+0;
             continue;
            }
            temp=fabs(aux_res2[i]);
            crit=tau[j]*((*Ts));
            if (temp >  crit) res_a[i]=res_a[i]+fabs(aux_res2[i]);//{
         }
      }
finner_prod_maxp(res_a, n, max_b, p, max_inner);
 }


void multi_across_fip(double *x_a, int *n,  int *d, double *tau,int *epp, double *Ts,int *min_draw, int *M, double *p,int *max_b,
                      double *max_inner,int *out1,double *out2,int *out3,int *out4) {
 int p1=0,p2=0,m,i,k=0,ilow=0,temp=0,len,j,len2=0, scales_temp;
 int Maux=*M, *len3=&len2, *p11=&p1;
 len = *n;
 scales_temp = *d;

GetRNGstate();
 for (m=0; m < Maux; m++){



        p1=  ceil(runif(1,(len-2) + 1));
        p2=  ceil(runif(p1,len));


   // p1= rand() % (len-2) + 1;
   //p2= rand() % (len-p1-1) + ((p1)+1);
      if (m==0) {
            p1=1+epp[(*d)-1];
            p2=len-epp[(*d)-1];
        }
    if ((p2-(p1))< *min_draw) continue;
    //Rprintf("%d \n",p2);
    out3[m]=p1;
    out4[m]=p2;
    temp=(p2-(p1)+1)*scales_temp;
    len2=p2-(p1)+1;
    double *dyn_arr = (double*)malloc(temp*sizeof(double));
    double *I_plus = (double*)calloc((p2-(p1)),sizeof(double));
    double *res_a = (double*)calloc((p2-(p1)),sizeof(double));
    double *aux_res = (double*)calloc((p2-(p1)+1),sizeof(double));
    double *aux_res2 = (double*)calloc((p2-(p1)),sizeof(double));
    double *I_minus = (double*)calloc((p2-(p1)),sizeof(double));

           k=0;
           for (j=0; j < scales_temp; j++){
           ilow = j * len;
           for (i=p1-1; i < p2; i++) {
           dyn_arr[k]=x_a[i+ilow];
           k=k+1;
           }
            }
            if (m==0) across_fip(dyn_arr, len3 , d, res_a, tau, epp, aux_res, aux_res2, I_plus, I_minus, Ts, p11,max_b, &p[0], max_inner);
            else across_fip(dyn_arr, len3 , d, res_a, tau, epp, aux_res, aux_res2, I_plus, I_minus, Ts, p11,max_b, &p[1], max_inner);
    out1[m]=*max_b+p1-1;
    out2[m]=*max_inner;
//Free memory
free(dyn_arr);
free(res_a);
free(I_plus);
free(aux_res);
free(aux_res2);
free(I_minus);

}
  PutRNGstate();
}


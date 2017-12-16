#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#define PI 3.141592653589793
#define SQH 0.707106781186547  /* square root of 2 */
#define SWAP(a,b)  tempr=(a); (a) = (b); (b) = tempr
#define IMAGE_SIZE 512
#define QF 50
#define MSB(a) a < 0 ? 1:0


static void fft1();

void dct1(x,n)
float *x;
int n;
{
  int i,ii,nn,mm;
  float tc,ts,sqn,temp2;
  double temp1;
  float *v;
  void fft1();

  nn = n >> 1;
  mm = n << 1;
  sqn = (float)sqrt((double)n);

  v = (float *) calloc (mm,sizeof(float));
  if (v == NULL) {
    printf("allocation failure\n");
    exit(1);
  }

  for (i=0;i<nn;i++) {
    ii = i << 1;
    v[ii] = x[ii];
    v[ii+1] = 0.0;
  }
  for (i=nn;i<n;i++) {
    ii = i << 1; 
    v[ii] = x[mm-ii-1];
    v[ii+1] = 0.0;
  }

  fft1(v-1,n,1);

  temp2 = SQH/sqn;
  x[0] = v[0]/sqn;
  for (i=1;i<=nn;i++) {
    ii = i << 1;
    temp1 = (double)(PI*i/mm);
    tc = (float) cos(temp1);
    ts = (float) sin(temp1);
    x[i] = 2.0*(tc*v[ii] + ts*v[ii+1])*temp2;
    x[n-i] = 2.0*(ts*v[ii] - tc*v[ii+1])*temp2;
  }

  free(v);
}
/* -------------------------------------------------- */

void idct1(x,n)
float *x;
int n;
{
  int i,ii,mm,nn;
  float *v;
  float temp2,tc,ts,sqn;
  double temp1;
  void fft1();

  nn = n >> 1;
  mm = n << 1;
  sqn = (float)sqrt((double)n);

  v = (float *) calloc (mm,sizeof(float));
  if (v == NULL) {
    printf("allocation failure\n");
    exit(1);
  }

  temp2 = sqn/SQH;  
  v[0] = x[0]*sqn;
  v[1] = 0.0;
  for (i=1;i<n;i++) {
    ii = i << 1;
    temp1 = (double)(PI*i/mm);
    tc = (float)cos(temp1);
    ts = (float)sin(temp1);
    v[ii] = 0.5*(tc*x[i] + ts*x[n-i])*temp2;
    v[ii+1] = 0.5*(ts*x[i] - tc*x[n-i])*temp2;
  }

  fft1(v-1,n,-1);

  for (i=0;i<nn;i++) {
    ii = i << 1;
    x[ii] = v[ii];
  }
  for (i=nn;i<n;i++) {
    ii = i << 1;
    x[mm-ii-1] = v[ii];
  }
  free(v);
}
/* -------------------------------------------------- */

/* 1-D fft program  */
/* Replace data by its discrete Fourier ransform if isign is input as 1,
   or replace data by its inverse discrete Fourier transform if
   isign is input as -1.  "data" is a complex array of length nn, input
   as a real array data [1..2*nn], nn must be an integer power of 2.     */

/* If your data array is zero-offset, that is the range of data is 
   [0..2*nn-1], you have to decrease the pointer to data by one when
   fft1 is invoked, for example fft1(data-1,256,1).
   The real part of the first output will now be return in data[0],
   the imaginary part in data[1] and so on.                              */

static void fft1(data,nn,isign)
float *data;
int nn,isign;
{
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;
    n = nn << 1;
    j = 1;
    for (i=1;i<n;i+=2) {
       if (j>i) {
          SWAP(data[j],data[i]);
          SWAP(data[j+1],data[i+1]);
       }
       m = n >> 1;
       while (m>=2 && j>m) {
         j -= m;
         m >>= 1;
       }
       j += m;
    }
    mmax = 2;
    while (n>mmax) {
         istep = 2*mmax;
         theta = 6.28318530717959/(-isign*mmax);
         wtemp = sin(0.5*theta);
         wpr = -2.0*wtemp*wtemp;
         wpi = sin(theta);
         wr = 1.0;
         wi = 0.0;
         for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
              j = i+mmax;
              tempr = wr*data[j]-wi*data[j+1];
              tempi = wr*data[j+1]+wi*data[j];
              data[j] = data[i]-tempr;
              data[j+1] = data[i+1]-tempi;
              data[i] += tempr;
              data[i+1] += tempi;
            }
            wr = (wtemp=wr)*wpr-wi*wpi+wr;
            wi = wi*wpr+wtemp*wpi+wi;
         }
         mmax = istep;
     }
    
     if (isign == -1) {
        for (i=1;i<=n;++i)
          data[i] = data[i]/nn;
     }
}

void dct2(x,n)
float **x;
int n;
{
  int i,j;
  float *y;
  void dct1();

  y = (float *) calloc (n,sizeof(float));
  if (y == NULL) {
   printf("allocation failure\n");
   exit(1);
  }
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      y[j] = x[j][i];
    dct1(y,n);
    for (j=0;j<n;j++) 
      x[j][i] = y[j];
  }   /* end of loop i */

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      y[j] = x[i][j];
    dct1(y,n);
    for (j=0;j<n;j++) 
      x[i][j] = y[j];
  }   /* end of loop i */

  free(y);
}

/* ----------------------------------------------- */

void idct2(x,n)
float **x;
int n;
{
  int i,j;
  float *y;
  void idct1();

  y = (float *) calloc (n,sizeof(float));
  if (y == NULL) {
   printf("allocation failure\n");
   exit(1);
  }
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      y[j] = x[j][i];
    idct1(y,n);
    for (j=0;j<n;j++) 
      x[j][i] = y[j];
  }   /* end of loop i */

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      y[j] = x[i][j];
    idct1(y,n);
    for (j=0;j<n;j++) 
      x[i][j] = y[j];
  }   /* end of loop i */

  free(y);
}

/* ----------------------------------------------- */

FILE *fout;
FILE *fin;
char temp_byte;
int read_counter = 0;
const int dcHuffmanValues[12] = { 0,2,3,4,5,6,14,30,62,126,254,510 };
const int dcHuffmanLength[12] = { 2,3,3,3,3,3,4,5,6,7,8,9 };

char get_bit();
int dc_get_ssss();
int get_diff_val(int);

int main(int argc,char **argv){

    char* input_fname = argv[1];

    fin = fopen(input_fname,"rb");
    fout = fopen(argv[2],"wb");
    int i,j;

    /*
    for(i=0 ; i< 15 ;i++){
        printf("bit:%d\n",get_bit());
    }
    */

}

char get_bit(){
    if(read_counter == 0){
        fread(&temp_byte, sizeof(char), 1, fin);
    }
    char result = MSB( temp_byte );
    temp_byte <<= 1;
    read_counter++;
    if(read_counter == 8){
        read_counter = 0;
    }
    return result;
}

int dc_get_ssss(){
    int tmp = 0;
    int code_len = 0;
    int arr_ptr = 0;
    int FLAG_SUCESS = 0;
    while(1)
    {
        tmp += get_bit();
        code_len++;
        while( code_len == dcHuffmanLength[arr_ptr]){
            if(tmp == dcHuffmanValues[arr_ptr]){
                FLAG_SUCESS = 1;
                break;
            }
            arr_ptr++;
        }
        tmp <<=1;
        if(FLAG_SUCESS)
            break;
    }
    return arr_ptr;
}

int get_diff_val(int ssss){
    int diff_tab_len = 1 << ssss;
    int diff_index = 0;
    int diff_val;
    int i;
    for(i = 0 ; i < ssss ; i++){
        diff_index <<= 1;
        diff_index += get_bit();
    }
    if( diff_index >= (diff_tab_len / 2) ){
        diff_val = diff_index;
    }
    else{
        diff_val = diff_index - diff_tab_len + 1;
    }

    if(ssss == 0)
        return 0;
    else
        return diff_val;
}
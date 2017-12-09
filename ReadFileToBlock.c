#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#define PI 3.141592653589793
#define SQH 0.707106781186547  /* square root of 2 */
#define SWAP(a,b)  tempr=(a); (a) = (b); (b) = tempr
#define IMAGE_SIZE 512
#define QF 50


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


unsigned char image_byte[IMAGE_SIZE][IMAGE_SIZE];

unsigned char s_matrix[8][8] = {{16,11,10,16,24,40,51,61},
                                {12,12,14,19,26,58,60,55},
                                {14,13,16,24,40,57,69,56},
                                {14,17,22,29,51,87,80,62},
                                {18,22,37,56,68,109,103,77},
                                {24,36,55,64,81,104,113,92},
                                {49,64,78,87,103,121,120,101},
                                {72,92,95,98,112,100,103,99}};
FILE *fout;

const int acHuffmanLength[256] = {
     4, 2, 2, 3, 4, 5, 7, 8,
    10,16,16, 0, 0, 0, 0, 0,
     0, 4, 5, 7, 9,11,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 5, 8,10,12,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 6, 9,12,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 6,10,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 7,11,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 7,12,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 8,12,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 9,15,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 9,16,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0, 9,16,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0,10,16,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0,10,16,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0,11,16,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
     0,16,16,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0,
    11,16,16,16,16,16,16,16,
    16,16,16, 0, 0, 0, 0, 0
};

const int acHuffmanTable[256] = {
    0x000a,0x0000,0x0001,0x0004,0x000b,0x001a,0x0078,0x00f8,
    0x03f6,0xff82,0xff83,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x000c,0x001b,0x0079,0x01f6,0x07f6,0xff84,0xff85,
    0xff86,0xff87,0xff88,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x001c,0x00f9,0x03f7,0x0ff4,0xff89,0xff8a,0xff8b,
    0xff8c,0xff8d,0xff8e,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x003a,0x01f7,0x0ff5,0xff8f,0xff90,0xff91,0xff92,
    0xff93,0xff94,0xff95,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x003b,0x03f8,0xff96,0xff97,0xff98,0xff99,0xff9a,
    0xff9b,0xff9c,0xff9d,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x007a,0x07f7,0xff9e,0xff9f,0xffa0,0xffa1,0xffa2,
    0xffa3,0xffa4,0xffa5,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x007b,0x0ff6,0xffa6,0xffa7,0xffa8,0xffa9,0xffaa,
    0xffab,0xffac,0xffad,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x00fa,0x0ff7,0xffae,0xffaf,0xffb0,0xffb1,0xffb2,
    0xffb3,0xffb4,0xffb5,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x01f8,0x7fc0,0xffb6,0xffb7,0xffb8,0xffb9,0xffba,
    0xffbb,0xffbc,0xffbd,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x01f9,0xffbe,0xffbf,0xffc0,0xffc1,0xffc2,0xffc3,
    0xffc4,0xffc5,0xffc6,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x01fa,0xffc7,0xffc8,0xffc9,0xffca,0xffcb,0xffcc,
    0xffcd,0xffce,0xffcf,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x03f9,0xffd0,0xffd1,0xffd2,0xffd3,0xffd4,0xffd5,
    0xffd6,0xffd7,0xffd8,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x03fa,0xffd9,0xffda,0xffdb,0xffdc,0xffdd,0xffde,
    0xffdf,0xffe0,0xffe1,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0x07f8,0xffe2,0xffe3,0xffe4,0xffe5,0xffe6,0xffe7,
    0xffe8,0xffe9,0xffea,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x0000,0xffeb,0xffec,0xffed,0xffee,0xffef,0xfff0,0xfff1,
    0xfff2,0xfff3,0xfff4,0x0000,0x0000,0x0000,0x0000,0x0000,
    0x07f9,0xfff5,0xfff6,0xfff7,0xfff8,0xfff9,0xfffa,0xfffb,
    0xfffc,0xfffd,0xfffe,0x0000,0x0000,0x0000,0x0000,0x0000
};

//最左上為0,0
//最右下為64,64
assign_block(float** target,int x,int y){
    int i,j;
    // minus 128
    for(i = 0; i < 8 ; i++){
        for(j = 0 ; j < 8 ; j++){
            target[i][j] = (float)image_byte[8*y+i][8*x+j] - 128 ;
        }
    }
}

float** get_block(int x,int y){
    //allocation for block
    float **block;
    block = (float **) malloc(sizeof(float *)*8);
    int i,j;
    for(i = 0 ; i < 8 ; i++){
        block[i] = (float *) malloc(sizeof(float) * 8);
    }
    assign_block(block,0,0);

    return block;
} 

quantization(float** block){
  float factor;
  if(QF < 50)
    factor = 5000 / QF;
  else if (QF >= 50)
    factor = 200 - 2 * QF;

  float q_matrix[8][8];

  int i,j;
  //printf("factor:%f\n",factor);
  for(i = 0; i < 8 ; i++){
        for(j = 0 ; j < 8 ; j++){
            q_matrix[i][j] = s_matrix[i][j] * factor / 100;
            block[i][j] = block[i][j] / q_matrix [i][j];
        }
  }

}

unsigned char byte_for_output = 0;
int byte_count = 0;

void output(unsigned char byte){
  fwrite(&byte,1,sizeof(byte),fout);
}

void bit_output(char bit){
  byte_for_output <<= 1;
  byte_for_output += bit - '0';
  byte_count++;
  // for output
  if(byte_count == 8)
  {
    output(byte_for_output);
    byte_count = 0;
    byte_for_output = 0;
  }
}

void output_preproc(unsigned int val,int len){
  char bit_str[32];
  int i;
  for(i = 31 ; i >= 0 ; i--){
    bit_str[i] = (val & 1) + '0';
    val >>= 1;
  }

  for(i = 32-len ; i < 32 ; i++ ){
    bit_output(bit_str[i]);
  }

}

int main(int argc,char **argv){

    // read and save a byte arr

    char* input_fname = argv[1];
    FILE *fin;

    fin = fopen(input_fname,"rb");
    fout = fopen(argv[2],"wb");
    int i,j;
    /*
    unsigned char temp;
    for( i = 0 ; i < IMAGE_SIZE ; i++){
        for( j = 0 ; j < IMAGE_SIZE ; j++){
            fread(&temp, sizeof(char), 1, fin);
            image_byte[i][j] = temp;
        }
    }
       
    float **block = get_block(0,0);
    
    dct2(block,8);

    quantization(block);

    for(i = 0; i < 8 ; i++){
        for(j = 0 ; j < 8 ; j++){
            printf("%f\t",* (*(block+i)+j));
        }
        printf("\n");
    }
    */

    output_preproc(0xff,8);
    output_preproc(0x08,4);
    output_preproc(0x00,4);


}



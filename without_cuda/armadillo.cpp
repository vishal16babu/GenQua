#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;
typedef struct Matrix
{
    int width;
    int height;
    double* elements;
} Matrix;
extern "C" void MatMul(Matrix A, Matrix B, Matrix C,int tid,int NUM_THREADS)
  {
  	//printf("created thread %d\n",tid );
  	  	if(NUM_THREADS==0)
{
	printf("SHIT\n");
}
  	int factor = (B.width)/NUM_THREADS;
  mat A_mat(A.height,A.width);
  mat B_mat(B.height,factor);
  mat C_mat(C.height,factor);
for(int i=0;i<A.height;i++){
for(int j=0;j<A.width;j++){
A_mat(i,j) = A.elements[i*A.width+j];
}
}
for(int i=0;i<B.height;i++){
for(int j=0;j<factor;j++){
B_mat(i,j) = B.elements[i*B.width+j+factor*tid];
}
}
C_mat = A_mat*B_mat;
for(int i=0;i<C.height;i++){
for(int j=0;j<factor;j++){
C.elements[i*C.width+j+(factor*tid)]=C_mat(i,j) ;
}
}

  return ;
  }
extern "C" void fillMul(Matrix A, Matrix B, Matrix C,int NUM_THREADS)
  {
    //printf("created thread %d\n",tid );
        if(NUM_THREADS==0)
{
  printf("SHIT\n");
}
    int factor = (B.width)/NUM_THREADS;
    int offset = B.width - (factor*NUM_THREADS);
  mat A_mat(A.height,A.width);
  mat B_mat(B.height,offset);
  mat C_mat(C.height,offset);
for(int i=0;i<A.height;i++){
for(int j=0;j<A.width;j++){
A_mat(i,j) = A.elements[i*A.width+j];
}
}
for(int i=0;i<B.height;i++){
for(int j=0;j<offset;j++){
B_mat(i,j) = B.elements[i*B.width+j+factor*NUM_THREADS];
}
}
C_mat = A_mat*B_mat;
for(int i=0;i<C.height;i++){
for(int j=0;j<offset;j++){
C.elements[i*C.width+j+(factor*NUM_THREADS)]=C_mat(i,j) ;
}
}

  return ;
  }

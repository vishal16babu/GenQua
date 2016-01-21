#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <pthread.h>
#include <cuda.h>
#define NUM_THREADs 5
#define BLOCK_SIZE 16
#define PI 3.141592654
#define MEGEXTRA 1000000

typedef struct Matrix
{
    int width;
    int height;
    double* elements;
} Matrix;


__global__ void find1elt(double *m, double *rs, int n)
{
    int rownum = blockIdx.x; // this thread will handle row # rownum
    double sum = 0;
    for (int k = 0; k < n; k++)
        sum += m[rownum*n+k];
    rs[rownum] = sum;
}
__global__ void Vector_Sub( double *dev_a , double *dev_b , double *dev_c,int n)
{
    //Get the id of thread within a block
//      unsigned short tid = threadIdx.x ;
    int tid=blockIdx.x * blockDim.x + threadIdx.x;
    if ( tid < n )  // check the boundry condition for the threads{
    {
        dev_c [tid] = (dev_a[tid]*1.0) - (dev_b[tid]*0.5) ;
// printf("%f %f\n",dev_a[tid],dev_b[tid]);
    }

}


extern "C" void get_host(Matrix A, Matrix D,double* Host_c , int n, int m,int dev_no){
    int nDevices;

  cudaGetDeviceCount(&nDevices);
  if(nDevices!=6){
    printf("shit\n");
    nDevices = 6;
  }
    dev_no = ((dev_no)%(nDevices));
        dev_no = nDevices - dev_no -1 ;

    cudaSetDevice(dev_no);
           double *dm, // device matrix
               *hrs, // host rowsums
               *drs; // device rowsums
        //int n;
        double msize = n * m * sizeof(double); // size of matrix in bytes
        // allocate space for host matrix
        // as a test, fill matrix with consecutive integers
        // allocate space for device matrix
        cudaMalloc((void **)&dm,msize);
        // copy host matrix to device matrix
        cudaMemcpy(dm,D.elements,msize,cudaMemcpyHostToDevice);
        // allocate host, device rowsum arrays
        double rssize = n * sizeof(double);
        hrs = (double *) malloc(rssize);
        cudaMalloc((void **)&drs,rssize);
        // set up parameters for threads structure
        dim3 dimGrid(n,1); // n blocks
        dim3 dimBlock(1,1,1); // 1 thread per block
        // invoke the kernel
        find1elt<<<dimGrid,dimBlock>>>(dm,drs,m);
        // wait for kernel to finish
        cudaThreadSynchronize();
        // copy row vector from device to host
        cudaMemcpy(hrs,drs,rssize,cudaMemcpyDeviceToHost);
        // check results
        /*
        printf("Sum(A*Q.xA),2)=\n");
        if (n < 100) for(int i=0; i<n; i++) printf("%f\n",hrs[i]);*/
        // clean up
        double *dm1, // device matrix
               *hrs1, // host rowsums
               *drs1; // device rowsums
        //int n;
        double msize1 = n * m * sizeof(double); // size of matrix in bytes
        // allocate space for host matrix
        // as a test, fill matrix with consecutive integers
        //  t = 0,i,j;

        // allocate space for device matrix
        cudaMalloc((void **)&dm1,msize1);
        // copy host matrix to device matrix
        cudaMemcpy(dm1,A.elements,msize1,cudaMemcpyHostToDevice);
        // allocate host, device rowsum arrays
        double rssize1 = n * sizeof(double);
        hrs1 = (double *) malloc(rssize1);
        cudaMalloc((void **)&drs1,rssize1);
        // set up parameters for threads structure
        dim3 dimGrid1(n,1); // n blocks
        dim3 dimBlock1(1,1,1); // 1 thread per block
        // invoke the kernel
        find1elt<<<dimGrid1,dimBlock1>>>(dm1,drs1,m);
        // wait for kernel to finish
        cudaThreadSynchronize();
        // copy row vector from device to host
        cudaMemcpy(hrs1,drs1,rssize1,cudaMemcpyDeviceToHost);
        // check results
        /* printf("Sum(A,2)=\n");
             if (n < 100) for(int i=0; i<n; i++) printf("%f\n",hrs1[i]);*/

        //Device array
        double* dev_c1 ;

        //Allocate the memory on the GPU
        cudaMalloc((void **)&dev_c1 , n*sizeof(double) ) ;


        dim3 dimGrid2(n,1); // n blocks
        dim3 dimBlock2(1,1,1); // 1 thread per block
        //Make a call to GPU kernel
        //Vector_Sub<<< dimGrid2 ,dimBlock2 >>> (dev_a1 , dev_b1 , dev_c1,n ) ;
        Vector_Sub<<< dimGrid2 ,dimBlock2 >>> (drs1 , drs , dev_c1,n ) ;
        //Copy back to Host array from Device array
        cudaMemcpy(Host_c , dev_c1 , n*sizeof(double) , cudaMemcpyDeviceToHost);

        //Display the result
        //printf("The F:=\n");
      /*  for ( int i = 0; i<n; i++ )
            printf ("%f\n", Host_c[i] ) ; */

        //Free the Device array memory
        cudaFree (dev_c1) ;

        cudaFree(dm);
        free(hrs);
        cudaFree(drs);
        cudaFree(dm1);
        free(hrs1);
        cudaFree(drs1);
        return;
}
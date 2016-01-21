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

__global__ void MatElMulKernel(Matrix, Matrix, Matrix);
extern "C" void MatElMul(Matrix AA, Matrix B, Matrix C, int dev_no)
{int nDevices;

  cudaGetDeviceCount(&nDevices);
  if(nDevices!=6){
    printf("shit\n");
    nDevices = 6;
  }
    dev_no = ((dev_no)%(nDevices));
    dev_no = nDevices - dev_no -1 ;
        cudaSetDevice(dev_no);

    // loading matrix A
    Matrix d_A;
    d_A.width=AA.width;
    d_A.height=AA.height;
    size_t size=AA.width * AA.height * sizeof(double);
    cudaMalloc(&d_A.elements, size);
    cudaMemcpy(d_A.elements, AA.elements, size, cudaMemcpyHostToDevice);
    //loading matrix B
    Matrix d_B;
    d_B.width=B.width;
    d_B.height=B.height;
    size= B.width * B.height * sizeof(double);
    cudaMalloc(&d_B.elements, size);
    cudaMemcpy(d_B.elements, B.elements, size, cudaMemcpyHostToDevice);

    // Allocate C in device memory
    Matrix d_C;
    d_C.width = C.width;
    d_C.height = C.height;
    size = C.width * C.height * sizeof(double);
    cudaMalloc(&d_C.elements, size);

    // Invoke kernel
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid((B.width + dimBlock.x - 1)/dimBlock.x, (AA.height + dimBlock.y - 1)/dimBlock.y);
    //printf("%d=\n",(B.width + dimBlock.x - 1)/dimBlock.x);
    //printf("%d=\n",(A.height + dimBlock.y - 1)/dimBlock.y);
    MatElMulKernel<<<dimGrid, dimBlock>>>(d_A, d_B, d_C);
    cudaThreadSynchronize();

    // Read C from device memory
    cudaMemcpy(C.elements, d_C.elements, size, cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_A.elements);
    cudaFree(d_B.elements);
    cudaFree(d_C.elements);
}

__global__ void MatElMulKernel(Matrix AA, Matrix B, Matrix C)
{
    double Cvalue = 0;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    if(row>AA.height||col>B.width)
        return;
//  for (int e = 0; e < A.width; ++e)
    Cvalue = (AA.elements[row * C.width+ col ]) * (B.elements[row*C.width + col]);
    C.elements[row * C.width + col] = Cvalue;

}
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <iostream>
#include <armadillo>
#define NUM_THREADS 48
using namespace std;
using namespace arma;
typedef struct Matrix
{
    int width;
    int height;
    double* elements;
} Matrix;
Matrix A2;
Matrix B2;
Matrix C2;
FILE* fptr2;

void MatelMul(Matrix A2, Matrix B2, Matrix C2,int tid)
  {
    //printf("created thread %d\n",tid );
        if(NUM_THREADS==0)
{
    printf("SHIT\n");
}
    int factor = (B2.width)/NUM_THREADS;
  mat A_mat(A2.height,A2.width);
  mat B_mat(B2.height,factor);
  mat C_mat(C2.height,factor);
for(int i=0;i<A2.height;i++){
for(int j=0;j<A2.width;j++){
A_mat(i,j) = A2.elements[i*A2.width+j];
}
}
for(int i=0;i<B2.height;i++){
for(int j=0;j<factor;j++){
B_mat(i,j) = B2.elements[i*B2.width+j+factor*tid];
}
}
C_mat = A_mat%B_mat;
for(int i=0;i<C2.height;i++){
for(int j=0;j<factor;j++){
C2.elements[i*C2.width+j+(factor*tid)]=C_mat(i,j) ;
}
}

  return ;
  }
 void fillelMul(Matrix A2, Matrix B2, Matrix C2)
  {
    //printf("created thread %d\n",tid );
        if(NUM_THREADS==0)
{
  printf("SHIT\n");
}
    int factor = (B2.width)/NUM_THREADS;
    int offset = B2.width - (factor*NUM_THREADS);
  mat A_mat(A2.height,A2.width);
  mat B_mat(B2.height,offset);
  mat C_mat(C2.height,offset);
for(int i=0;i<A2.height;i++){
for(int j=0;j<A2.width;j++){
A_mat(i,j) = A2.elements[i*A2.width+j];
}
}
for(int i=0;i<B2.height;i++){
for(int j=0;j<offset;j++){
B_mat(i,j) = B2.elements[i*B2.width+j+factor*NUM_THREADS];
}
}
C_mat = A_mat%B_mat;
for(int i=0;i<C2.height;i++){
for(int j=0;j<offset;j++){
C2.elements[i*C2.width+j+(factor*NUM_THREADS)]=C_mat(i,j) ;
}
}

  return ;
  }



void *form_el_partial(void *threadid)
{
    long tid;
    tid = (long)threadid;
    MatelMul(A2,B2,C2,tid);


    pthread_exit((void*) threadid);
}


extern "C" void matelmulti(Matrix AA, Matrix BB, Matrix CC)
{

    //cout<<y.n_rows<<" "<<y.n_cols<<endl;
         int m=5000;
int factor = AA.height/m;
int offset = AA.height - factor*m;
    int n=AA.width;
    int l=BB.width;
    int i;
    int j;
    fptr2 =fopen("Q_matrix_fin2.txt","w");
    fseek(fptr2, 0, SEEK_SET);
B2 = BB;
int num_elts = 0;
    for(int iter =0;iter<factor;iter++){
    A2.height = m;
    A2.width = n;
    A2.elements = (double*)malloc(A2.height*A2.width*sizeof(double));
   
    time_t t2;
    srand(time(NULL));
    for(int i=0;i<A2.height;i++){
for(int j=0;j<A2.width;j++){
 A2.elements[i*A2.width+j]=AA.elements[(i+m*iter)*AA.width+j];
}
}

    //print_matrix(A);
    //print_matrix(B);
    C2.height = A2.height ;
    C2.width = B2.width ;
    C2.elements = (double*)malloc(C2.height*C2.width*sizeof(double)); 

//printf("here\n");

    
 struct timespec start, finish;
 double elapsed;
 clock_gettime(CLOCK_MONOTONIC, &start);
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    size_t stacksize;
    int rc;
    long t;
    void *status;

    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    long temp_t;
    for(t=0; t<NUM_THREADS; t++)
    {
        temp_t = t;
        //printf("In main: creating thread %ld\n", temp_t);
        rc = pthread_create(&threads[temp_t], &attr,form_el_partial, (void *)temp_t);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            //exit(-2);
        }

    }
    pthread_attr_destroy(&attr);
    for(t=0; t<NUM_THREADS; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            //exit(-2);
        }
        //printf("Main: completed join with thread %ld having a status of %ld\n",t,(long)status);
    }
    fillelMul(A2,B2,C2);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 2000000000.0;
         
              printf("%f\n", elapsed);
//printf("hereagainn\n");

//print_matrix(C);
    free(A2.elements);
     for(int i=0;i<C2.height;i++){
for( j=0;j<C2.width;j++){
     CC.elements[i*C2.width+j + num_elts] = C2.elements[i*C2.width+j];
}

}
num_elts = num_elts + C2.height*C2.width;
        free(C2.elements);

}

A2.height = offset;
    A2.width = n;
    A2.elements = (double*)malloc(A2.height*A2.width*sizeof(double));
   
    time_t t2;
    srand(time(NULL));
     for(int i=0;i<A2.height;i++){
for(int j=0;j<A2.width;j++){
 A2.elements[i*A2.width+j]=AA.elements[(i+m*factor)*AA.width+j];
}
}

    //print_matrix(A);
    //print_matrix(B);
    C2.height = A2.height ;
    C2.width = B2.width ;
    C2.elements = (double*)malloc(C2.height*C2.width*sizeof(double)); 

printf("here\n");

    
 struct timespec start, finish;
 double elapsed;
 clock_gettime(CLOCK_MONOTONIC, &start);
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    size_t stacksize;
    int rc;
    long t;
    void *status;

    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    long temp_t;
    for(t=0; t<NUM_THREADS; t++)
    {
        temp_t = t;
        //printf("In main: creating thread %ld\n", temp_t);
        rc = pthread_create(&threads[temp_t], &attr,form_el_partial, (void *)temp_t);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            //exit(-2);
        }

    }
    pthread_attr_destroy(&attr);
    for(t=0; t<NUM_THREADS; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            //exit(-2);
        }
        //printf("Main: completed join with thread %ld having a status of %ld\n",t,(long)status);
    }
    fillelMul(A2,B2,C2);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 2000000000.0;
         
              printf("%f\n", elapsed);
//printf("hereagainn\n");

//print_matrix(C);
    free(A2.elements);
     for(int i=0;i<C2.height;i++){
for( j=0;j<C2.width;j++){
    //printf("%f\n", C.elements[i*C.width+j]);
     CC.elements[i*C2.width+j + num_elts] = C2.elements[i*C2.width+j];
         //printf("%f\n", CC.elements[i*C.width+j + num_elts] );

}

}
num_elts = num_elts + C2.height*C2.width;
        free(C2.elements);

    free(B2.elements);
fclose(fptr2);


    return ;

}

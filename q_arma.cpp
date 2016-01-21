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
Matrix A1;
Matrix B1;
Matrix C1;
FILE* fptr;

void MatMulti(Matrix A1, Matrix B1, Matrix C1,int tid)
  {
    //printf("created thread %d\n",tid );
        if(NUM_THREADS==0)
{
    printf("SHIT\n");
}
    int factor = (B1.width)/NUM_THREADS;
  mat A_mat(A1.height,A1.width);
  mat B_mat(B1.height,factor);
  mat C_mat(C1.height,factor);
for(int i=0;i<A1.height;i++){
for(int j=0;j<A1.width;j++){
A_mat(i,j) = A1.elements[i*A1.width+j];
}
}
for(int i=0;i<B1.height;i++){
for(int j=0;j<factor;j++){
B_mat(i,j) = B1.elements[i*B1.width+j+factor*tid];
}
}
C_mat = A_mat*B_mat;
for(int i=0;i<C1.height;i++){
for(int j=0;j<factor;j++){
C1.elements[i*C1.width+j+(factor*tid)]=C_mat(i,j) ;
}
}

  return ;
  }
 void fillMul(Matrix A1, Matrix B1, Matrix C1)
  {
    //printf("created thread %d\n",tid );
        if(NUM_THREADS==0)
{
  printf("SHIT\n");
}
    int factor = (B1.width)/NUM_THREADS;
    int offset = B1.width - (factor*NUM_THREADS);
  mat A_mat(A1.height,A1.width);
  mat B_mat(B1.height,offset);
  mat C_mat(C1.height,offset);
for(int i=0;i<A1.height;i++){
for(int j=0;j<A1.width;j++){
A_mat(i,j) = A1.elements[i*A1.width+j];
}
}
for(int i=0;i<B1.height;i++){
for(int j=0;j<offset;j++){
B_mat(i,j) = B1.elements[i*B1.width+j+factor*NUM_THREADS];
}
}
C_mat = A_mat*B_mat;
for(int i=0;i<C1.height;i++){
for(int j=0;j<offset;j++){
C1.elements[i*C1.width+j+(factor*NUM_THREADS)]=C_mat(i,j) ;
}
}

  return ;
  }



void *form_partial(void *threadid)
{
    long tid;
    tid = (long)threadid;
    MatMulti(A1,B1,C1,tid);


    pthread_exit((void*) threadid);
}
void print_matrix(Matrix x)
{
    int i;
    int j;
    for(i=0; i<(x.height); i++)
    {
        for(j=0; j<(x.width); j++)
        {
            printf("%lf\t",(x.elements[i*(x.width)+(j)]));
        }
        printf("\n");
    }
    printf("\n");
}

extern "C" void MatMul(Matrix AA, Matrix BB, Matrix CC)
{

    //cout<<y.n_rows<<" "<<y.n_cols<<endl;
         int m=5000;
int factor = AA.height/m;
int offset = AA.height - factor*m;
    int n=AA.width;
    int l=BB.width;
    int i;
    int j;
    fptr =fopen("Q_matrix_fin2.txt","w");
    fseek(fptr, 0, SEEK_SET);
B1 = BB;
int num_elts = 0;
    for(int iter =0;iter<factor;iter++){
    A1.height = m;
    A1.width = n;
    A1.elements = (double*)malloc(A1.height*A1.width*sizeof(double));
   
    time_t t1;
    srand(time(NULL));
    for(int i=0;i<A1.height;i++){
for(int j=0;j<A1.width;j++){
 A1.elements[i*A1.width+j]=AA.elements[(i+m*iter)*AA.width+j];
}
}

    //print_matrix(A);
    //print_matrix(B);
    C1.height = A1.height ;
    C1.width = B1.width ;
    C1.elements = (double*)malloc(C1.height*C1.width*sizeof(double)); 

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
        rc = pthread_create(&threads[temp_t], &attr,form_partial, (void *)temp_t);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            //exit(-1);
        }

    }
    pthread_attr_destroy(&attr);
    for(t=0; t<NUM_THREADS; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            //exit(-1);
        }
        //printf("Main: completed join with thread %ld having a status of %ld\n",t,(long)status);
    }
    fillMul(A1,B1,C1);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
         
              printf("%f\n", elapsed);
//printf("hereagainn\n");

//print_matrix(C);
    free(A1.elements);
     for(int i=0;i<C1.height;i++){
for( j=0;j<C1.width;j++){
     CC.elements[i*C1.width+j + num_elts] = C1.elements[i*C1.width+j];
}

}
num_elts = num_elts + C1.height*C1.width;
        free(C1.elements);

}

A1.height = offset;
    A1.width = n;
    A1.elements = (double*)malloc(A1.height*A1.width*sizeof(double));
   
    time_t t1;
    srand(time(NULL));
     for(int i=0;i<A1.height;i++){
for(int j=0;j<A1.width;j++){
 A1.elements[i*A1.width+j]=AA.elements[(i+m*factor)*AA.width+j];
}
}

    //print_matrix(A);
    //print_matrix(B);
    C1.height = A1.height ;
    C1.width = B1.width ;
    C1.elements = (double*)malloc(C1.height*C1.width*sizeof(double)); 

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
        rc = pthread_create(&threads[temp_t], &attr,form_partial, (void *)temp_t);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            //exit(-1);
        }

    }
    pthread_attr_destroy(&attr);
    for(t=0; t<NUM_THREADS; t++)
    {
        rc = pthread_join(threads[t], &status);
        if (rc)
        {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            //exit(-1);
        }
        //printf("Main: completed join with thread %ld having a status of %ld\n",t,(long)status);
    }
    fillMul(A1,B1,C1);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
         
              printf("%f\n", elapsed);
//printf("hereagainn\n");

//print_matrix(C);
    free(A1.elements);
     for(int i=0;i<C1.height;i++){
for( j=0;j<C1.width;j++){
    //printf("%f\n", C.elements[i*C.width+j]);
     CC.elements[i*C1.width+j + num_elts] = C1.elements[i*C1.width+j];
         //printf("%f\n", CC.elements[i*C.width+j + num_elts] );

}

}
num_elts = num_elts + C1.height*C1.width;
        free(C1.elements);

    free(B1.elements);
fclose(fptr);


    return ;

}

extern "C" void MatElMul(Matrix AA, Matrix BB, Matrix CC)
{
    for(int i=0;i<CC.height;i++){
        for(int j=0;j<CC.width;j++){
            CC.elements[i*CC.width+j] = AA.elements[i*CC.width+j]*BB.elements[i*CC.width+j];
        }
    }

}
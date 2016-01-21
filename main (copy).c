#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <pthread.h>
#include <cuda.h>
#include <string.h>
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

typedef struct pack
{
    int index;
    
    double value;
} pack;

double** Pool;
int pool_size = 2000;
Matrix A;
pthread_attr_t attr;
int np;
int nn;
int NUM_THREADS;
int* range;
FILE* fit_acc;


void get_host(Matrix A, Matrix D,double* Host_c , int n, int m, int dev_no);
void MatMul(Matrix AA, Matrix B, Matrix C, int dev_no);
void MatElMul(Matrix AA, Matrix B, Matrix C, int dev_no);

double linear_kernel(double* x_i , double* x_j, int d)
{
    double sum=0;
    int i;
    for(i=0; i<d; i++)
    {
        sum=sum+(x_i[i]*x_j[i]);
    }
    return sum;
}

double matrix_linear_kernel(double* x_elem , double* x_j,int row, int d)
{
    double sum=0;
    int i;
    for(i=0; i<d; i++)
    {
        sum=sum+(x_elem[row*(d)+i]*x_j[i]);
    }
    return sum;
}

double gaussian_kernel(double* x_i , double* x_j, int d)
{
    double sum=0;
    int i;
    for(i=0; i<d; i++)
    {
        sum=sum+((x_i[i]-x_j[i])*(x_i[i]-x_j[i]));
    }
    sum = exp((-1.0)*((1.0)/(d*1.0))*(sum));
    return sum;
}

void transpose_matrix( Matrix* AA , Matrix* A_T)
{
    int i;
    int j;
    for(i=0; i<(AA->height); i++)
    {
        for(j=0; j<(AA->width); j++)
        {
            (A_T->elements)[j*(AA->height)+i] = (AA->elements)[i*(AA->width) + j];
        }
    }
    return;
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
}




void uniform_distinct_rand(int * AA , int M , int N)
{

    unsigned char is_used[N];
    /* flags */
    int ii;
    for(ii=0; ii<N; ii++)
    {
        is_used[ii]=0;
    }
    int in, im;
    im = 0;

    for (in = N - M; in < N && im < M; ++in)
    {
        int r = rand() % (in + 1); /* generate a random number 'r' */

        if (is_used[r])
            /* we already have 'r' */
            r = in; /* use 'in' instead of the generated number */

        //assert(!is_used[r]);
        AA[im++] = r ; /* +1 since your range begins from 1 */
        is_used[r] = 1;
    }
    return;
}

void get_row(double* AA,int n ,double* arr, int row)
{
    int i;
    for(i=0; i<n; i++)
    {
        arr[i]=AA[row*n+i];
    }
    return;
}

double sum_of_elements(double* arr,int start,int end)
{
   // printf("sum_of_elements\n");
    double sum;
    sum=0;
    int i;

    for( i=start; i<end; i++)
    {
        sum+=arr[i];
    }
    return sum;
}

double gaussrand()
{
    static double U, V;
    static int phase = 0;
    double Z;

    if(phase == 0)
    {
        U = (rand() + 1.) / (RAND_MAX + 2.);
        V = rand() / (RAND_MAX + 1.);
        Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
    }
    else
    {
        Z = sqrt(-2 * log(U)) * cos(2 * PI * V);
    }

    phase = 1 - phase;

    return Z;
}

void fill_up_A(double* AA,int row,int d)
{
    int C= rand()%d + 1;
   // printf("C is %d\n", C);
   // printf("n is %d\n", (np+nn));
   // printf("d is %d\n", d);

    int sp= rand()%np + 1;
    int rp[sp];
    double alpha_p[sp];
    int i;
    /* for(i=0; i<sp; i++)
     {
         rp[i]=(int)(rand()%np);
     }*/
    uniform_distinct_rand(rp,sp,np);
    double sum=0;
    for(i=0; i<sp; i++)
    {
        alpha_p[i]=fabs(gaussrand());
        sum+=alpha_p[i];
    };
    double fp=((np+nn)*C*1.0)/(4*sum);
    for(i=0; i<sp; i++)
    {
        alpha_p[i]=alpha_p[i]*fp;
        AA[row*(nn+np) + rp[i]]=alpha_p[i];
    };


    int sn= rand()%nn + 1;
    int rn[sn];
    double alpha_n[sn];
    /* for(i=0; i<sn; i++)
     {
         rn[i]=(int)(rand()%nn);
     }*/
    uniform_distinct_rand(rn,sn,nn);

    sum=0;
    for(i=0; i<sn; i++)
    {
        alpha_n[i]=fabs(gaussrand());
        sum+=alpha_n[i];
    };
    double fn=((np+nn)*C*1.0)/(4*sum);
    for(i=0; i<sn; i++)
    {
        alpha_n[i]=alpha_n[i]*fn;
        AA[row*(nn+np) + np + rn[i]]=alpha_n[i];
    };


    return;
}

double average_of(double* F, int size)
{
    double avg = (sum_of_elements(F,0,size)/size);
    return avg;
}


void balance(double* alpha)
{
    double sum1 = sum_of_elements(alpha,0,np);
    double sum2 = sum_of_elements(alpha,np,np+nn);
    if(sum2>sum1)
    {
        double epsilon = (sum2-sum1)*0.5 ;


        int l_p = (int)(rand()%np);
        alpha[l_p] = alpha[l_p]+epsilon;
        while(epsilon>0.0000001)
        {
            int l_n = (int)(rand()%nn);
            if(alpha[np+l_n]>=epsilon)
            {
                alpha[np+l_n]=alpha[np+l_n]-epsilon;
                epsilon=0;
                break;
            }
            else
            {
                epsilon=epsilon - alpha[np+l_n];
                alpha[np+l_n]=0;

            }
        }
    }
    else
    {

        double epsilon = (sum1-sum2)*0.5 ;

        int l_n = (int)(rand()%nn);
        alpha[np+l_n] = alpha[np+l_n]+epsilon;
        while(epsilon>0.0000001)
        {
            int l_p = (int)(rand()%np);
            if(alpha[l_p]>=epsilon)
            {
                alpha[l_p]=alpha[l_p]-epsilon;
                epsilon=0;
                break;
            }
            else
            {
                epsilon=epsilon - alpha[l_p];
                alpha[l_p]=0;

            }
        }


    }
   
    return;
}







void *parallel_crossover(void *threadid)
{
   long tid;
   tid = (long)threadid;
   //printf("Hello World! It's me, thread #%ld!\n", tid);
  /* int* temp_row_arr;
   temp_row_arr = (int*)malloc(2*sizeof(int));
if(temp_row_arr ==NULL)
{
printf("temp_row_arr/n");
}

//printf("temp_row_arr1/n");
   uniform_distinct_rand(temp_row_arr,2,pool_size);*/
int a;int i;int index;int j;
a = ((int)rand())%pool_size;
for(i=0;i<A.height;i++){
    if(a<=range[i]){
        index = i;
        break;
    }
}
int right = range[index];
int left;
if(index==0){
     left = 0;
}
else{
     left = range[index-1];
}
int offset = right - left;
int temp_range = pool_size - offset;
int b;
do{
 b = ((int)rand())%temp_range;
if(b<=left){
    b=b;
}
else{
    b = b+offset;
}

j=0;
for(i=0;i<A.width;i++){
    if(Pool[a][i]==Pool[b][i]){
        j++;
    }
}
}while(j==A.width);
double* old1;
old1= (double*)malloc(A.width*sizeof(double));
double* old2;
old2= (double*)malloc(A.width*sizeof(double));
double* new1;
new1= (double*)malloc(A.width*sizeof(double));
double* new2;
new2= (double*)malloc(A.width*sizeof(double));
double* new3;
new3= (double*)malloc(A.width*sizeof(double));
double* new4;
new4= (double*)malloc(A.width*sizeof(double));

for(i=0;i<A.width;i++){
    old1[i]=Pool[a][i];
    old2[i]=Pool[b][i];
}
j=0;

int r_p = (int)(rand()%np);
    int *kp;
    kp = (int*)malloc(r_p * sizeof(int));
    uniform_distinct_rand(kp , r_p , np);
    double temp_swap;
    for(i=0; i<r_p; i++)
    {
        temp_swap = old1[kp[i]];
        old1[kp[i]] = old2[kp[i]];
        old2[kp[i]] = temp_swap;
    }
    free(kp);

    int r_n = (int)(rand()%nn);
    int *kn;
    kn = (int*)malloc(r_n * sizeof(int));
    uniform_distinct_rand(kn , r_n , nn);
    for(i=0; i<r_n; i++)
    {
        temp_swap = old1[np+kn[i]];
        old1[np+kn[i]] = old2[np+kn[i]];
        old2[np+kn[i]] = temp_swap;
    }
    free(kn);
for(i=0; i<np; i++)
    {
        new1[i]=old1[i];
    }
    for(i=0; i<nn; i++)
    {
        new1[i+np]=old1[i+np];
    }
    balance(new1);

    for(i=0; i<np; i++)
    {
        new2[i]=old1[i];
    }
    for(i=0; i<nn; i++)
    {
        new2[i+np]=old2[i+np];
    }
    balance(new2);

    for(i=0; i<np; i++)
    {
        new3[i]=old2[i];
    }
    for(i=0; i<nn; i++)
    {
        new3[i+np]=old1[i+np];
    }
    balance(new3);

    for(i=0; i<np; i++)
    {
        new4[i]=old2[i];
    }
    for(i=0; i<nn; i++)
    {
        new4[i+np]=old2[i+np];
    }
    balance(new4);
for(i=0;i<A.width;i++){
    A.elements[(((4*tid)+2)*A.width)+i] = new1[i];
}
for(i=0;i<A.width;i++){
    A.elements[(((4*tid)+3)*A.width)+i] = new2[i];
}
for(i=0;i<A.width;i++){
    A.elements[(((4*tid)+4)*A.width)+i] = new3[i];
}
for(i=0;i<A.width;i++){
    A.elements[(((4*tid)+5)*A.width)+i] = new4[i];
}
free(old1);
free(old2);
free(new1);
free(new2);
free(new3);
free(new4);
//free(temp_row_arr);
   pthread_exit((void*) threadid);
}



void print_accuracy(double* arr,double* y_elem,double* x_elem,double**v,int v_size,int m,int d){
    double ans=0;
    int val =0;
    int i;int j;
    for(j=0;j<v_size;j++){
        ans = 0;
for(i=0;i<m;i++){
    ans = ans + ((arr[i])*(y_elem[i])*(matrix_linear_kernel(x_elem,v[j],i,d)));
}
if(ans*v[j][d]>=0){
    val++;
}
}
ans = (val*1.0)/(v_size*1.0);
printf("best accuracy is\t%f\n",ans);
return;
}


void mat_accuracy(Matrix* x , Matrix* v , double* arr ,double* y, Matrix* acc ,int rank){
    int i;
Matrix temp;
temp.width =np+nn;
temp.height = 1 ;
temp.elements = (double*)malloc(temp.width*sizeof(double));
for(i=0;i<temp.width;i++){
    temp.elements[i] = arr[i]*y[i];
}
Matrix C;
 C.height = x->height;
    C.width = v->width;
    C.elements = (double*)malloc(C.width * C.height * sizeof(double));
MatMul(*x,*v,C,rank);
Matrix P;
 P.height = temp.height;
    P.width = C.width;
    P.elements = (double*)malloc(P.width * P.height * sizeof(double));
    MatMul(temp,C,P,rank);
    double val =0;
    for(i=0;i<P.width;i++){
if((acc->elements[i]*P.elements[i]>=0)){
val = val+1.0;
}
    }
    val = (val/(P.width*1.0));
   // fprintf(fit_acc,"%f\n",val);
    free(temp.elements);
    free(C.elements);
    free(P.elements);
    return ;
}





void rescale_to_standard(double* arr,int size){
int i;double min;double max;double sum;
sum=0;
min = arr[0];
max= arr[0];
for(i=0;i<size;i++){
    if(arr[i]>max){
        max = arr[i];
    }
    sum+=arr[i];
}
for(i=0;i<size;i++){
    if(arr[i]<min){
        min = arr[i];
    }
}
double diff = max - min;
if(diff!=0){
for(i=0;i<size;i++){
    arr[i] = ((arr[i]-min)/diff);
}
}
else{
    for(i=0;i<size;i++){
        arr[i] = 1.0;
    }
}
   //     fprintf(fit_acc,"%f\n",((sum)/(size)));
    //    fprintf(fit_acc,"%f\n",max);


return;
}


int comparator(const void *p, const void *q) 
{
    int l = ((struct pack *)p)->value;
    int r = ((struct pack *)q)->value; 
    return (r - l);
}



int main(int argc, char** argv)
{





int rank;int numtasks;

MPI_Status stat;
MPI_Datatype rowtype;

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

printf("rank is %d\n", rank);
    FILE * file_p = fopen("../data_xp2.csv","r");
    FILE * file_n = fopen("../data_xn2.csv","r");
    FILE * file_v = fopen("../data_vt2.csv","r");
/*FILE* fptr =fopen("output2/time_table1.txt","w");
  fit_acc= fopen("output2/fitness_table1.txt","w");
  FILE* alpha_array = fopen("output2/bigshake.txt","w");
  fseek(fptr, 0, SEEK_SET);
  fseek(fit_acc,0,SEEK_SET);
     fseek(alpha_array,0,SEEK_SET);
        if(fptr ==  NULL){
            printf("file time opening failed\n");
            return 0;
        }*/
    if(file_p == NULL)
    {
        printf("\n file_p opening failed");
        return 0;
    }
    if(file_n == NULL)
    {
        printf("\n file_n opening failed");
        return 0;
    }
    if(file_v == NULL)
    {
        printf("\n file_v opening failed");
        return 0;
    }
    int d;
    np=0;
    nn=0;
    d=1;
    int v_size;
    v_size =0;
    char line[1000000];
    fgets(line, sizeof(line), file_p) ;
    int i;
    int j;
    for(i=0; i<1000000; i++)
    {
        if(line[i]==',')
        {
            d=d+1;
        }
        if(line[i]=='\n'||line[i]=='\r')
        {
            break;
        }
    }
    fseek(file_p, 0, SEEK_SET);
    while(fgets(line, sizeof(line), file_p))
    {
        np=np+1;
    }
    fseek(file_n, 0, SEEK_SET);
    while(fgets(line, sizeof(line), file_n))
    {
        nn=nn+1;
    }
    fseek(file_v, 0, SEEK_SET);
    while(fgets(line, sizeof(line), file_v))
    {
        v_size=v_size+1;
    }
    Matrix* x;
    x = (Matrix *)malloc(sizeof(Matrix));
    x->width = d;
    x->height = np+nn;

    x->elements = (double *)malloc(d *(np+nn) * sizeof(double));

  /*  v = (double**)malloc((v_size)*sizeof(double*));
    for(i=0;i<v_size;i++){
        v[i] = (double*)malloc((d+1)*sizeof(double));
    }*/
    Matrix* v;
    v = (Matrix *)malloc(sizeof(Matrix));
    v->width = d;
    v->height =v_size ;
    v->elements = (double *)malloc(v_size*d* sizeof(double));

Matrix* acc;
  acc = (Matrix *)malloc(sizeof(Matrix));
    acc->width = v_size;
    acc->height =1 ;
    acc->elements = (double *)malloc(v_size* sizeof(double));
    char* record;
    i=0;
    fseek(file_p, 0, SEEK_SET);

    while(fgets(line, sizeof(line),file_p)&&i<np)
    {
        record = strtok(line,",");
        j=0;
        while(record != NULL)
        {
            sscanf(record, "%lf", &(x->elements[i*(x->width)+(j++)]));
            record = strtok(NULL,",");
        }
        free(record);
        ++i;
    }
    fseek(file_n, 0, SEEK_SET);
    i=np;
    while(fgets(line, sizeof(line),file_n)&&i<(np+nn))
    {
        record = strtok(line,",");
        j=0;
        while(record != NULL)
        {
            sscanf(record, "%lf", &(x->elements[i*(x->width)+(j++)]));
            record = strtok(NULL,",");
        }
        free(record);
        ++i;
    }
    fseek(file_v, 0, SEEK_SET);

    i=0;
    while(fgets(line, sizeof(line),file_v)&&i<(v_size))
    {
        record = strtok(line,",");
        j=0;
        while(record != NULL)
        {
            if(j==d){
                sscanf(record, "%lf", &(acc->elements[i]));
                record = strtok(NULL,",");
            }
            else{
            sscanf(record, "%lf", &(v->elements[i*(v->width)+(j++)]));
            record = strtok(NULL,",");
        }
        }
        free(record);
        ++i;
    }
   // printf("j is %d\n", j);
    /*
    for(i=0;i<(np+nn);i++){
    for(j=0;j<d;j++){
    printf("%lf\t",x[i][j]);
    }
    printf("\n");
    }*/
    // printf("x is \n" );
    // print_matrix(*x);
    Matrix* x_t;
    x_t = (Matrix *)malloc(sizeof(Matrix));
    x_t->width = np+nn;
    x_t->height = d;

    x_t->elements = (double *)malloc(d *(np+nn) * sizeof(double));
    transpose_matrix(x , x_t);
     Matrix* v_t;
    v_t = (Matrix *)malloc(sizeof(Matrix));
    v_t->width = v_size;
    v_t->height = d;

    v_t->elements = (double *)malloc(d *v_size * sizeof(double));
    transpose_matrix(v , v_t);
// printf("x_t is \n" );
    // print_matrix(*x_t);

//printf("%lf\n",x[np][0]);
    fclose(file_p);
    fclose(file_n);
    fclose(file_v);
///created the array x[np+nn][d]
    /*
        */
///created the array y[np+nn]
    /**/
    int n;
    n=np+nn;
    Matrix Q;
    Q.width = n;
    Q.height = n;
    Q.elements = (double *)malloc(n * n * sizeof(double));


    Matrix* K_X;
    K_X =  (Matrix *)malloc(sizeof(Matrix));
    K_X->width = n;
    K_X->height = n;
    K_X->elements = (double *)malloc(n * n * sizeof(double));

    MatMul(*x , *x_t , *K_X ,rank);
//printf("K_X is \n" );
   // free(x);
    free(x_t->elements);
    free(x_t);
    // print_matrix(K_X);
    double y_value_for_p = 1.0;
    double y_value_for_n = -1.0;
    Matrix* y;
    y = (Matrix *)malloc( sizeof(Matrix));
    y->width= (np+nn) ;
    y->height = 1;
    y->elements = (double *)malloc((np+nn) * sizeof(double));
    for(i=0; i<np; i++)
    {
        (y->elements)[i]=y_value_for_p;
    }
    for(i=0; i<nn; i++)
    {
        (y->elements)[i+np]=y_value_for_n;
    }


    Matrix* y_t;
    y_t = (Matrix *)malloc( sizeof(Matrix));
    y_t->width= 1;
    y_t->height = (np+nn) ;
    y_t->elements = (double *)malloc((np+nn) * sizeof(double));
    transpose_matrix(y,y_t);
    Matrix* K_Y;
    K_Y =  (Matrix *)malloc(sizeof(Matrix));
    K_Y->width = n;
    K_Y->height = n;
    K_Y->elements = (double *)malloc(n * n * sizeof(double));


    MatMul(*y , *y_t , *K_Y ,rank);
//printf("K_Y is \n" );
    // print_matrix(K_Y);
   // free(y);
    free(y_t->elements);
    free(y_t);
    MatElMul(*K_Y , *K_X ,   Q,rank);
    free(K_X);
    free(K_Y);
   
struct 
    {
        double values[np+nn];
    } alpha;

    MPI_Datatype alpha_type;

     MPI_Type_contiguous ( (np+nn), MPI_DOUBLE, &alpha_type );
     MPI_Type_commit ( &alpha_type );
if(rank!=0){
    int jj= 0;
    int m = 100;
    int string_arg_m = m;
    NUM_THREADS = ((m-4)/4);
    A.height=m;
    A.width=np+nn;
    A.elements = (double*)malloc(A.width * A.height * sizeof(double));
    
      Pool = (double**)malloc(pool_size*sizeof(double*));
    for(i=0; i<pool_size; i++)
    {
        Pool[i] = (double*)malloc(A.width*sizeof(double));
    }
// print_matrix(Q);
        Matrix C,D;

 C.height = A.height;
    C.width = A.width;
    C.elements = (double*)malloc(C.width * C.height * sizeof(double));
    D.height = A.height;
    D.width = A.width;
    D.elements = (double*)malloc(D.width * D.height * sizeof(double));
            double *Host_c=(double*)malloc(A.height*sizeof(double));
            range = (int*)malloc(A.height*sizeof(int));
   
    double temp;
 
    int g = 10;
 
   
int iter;
clock_t start;
clock_t end;
double cpu_time_used ;
    for(jj=0;jj<5;jj++){

    
for( i=0; i<(A.height * A.width); i++)
    {
        A.elements[i]=0;
    }

    time_t t1;
    srand ( time(NULL) );

    int row;
    for( row=0; row<m; ++row)
    {
        fill_up_A(A.elements,row,d);
    };
for(iter = 0;iter<1;iter++){

start = clock();
    int gi;
    for(gi = 0; gi<g; gi++)
    {//printf("%d\n", gi);
   // cudaSetDevice(1);

        int i;
        MatMul(A, Q, C,rank);
        /*
            printf("Matrix A=\n");
            for(int i=0;i<A.height;i++){
                for(int j=0;j<A.width;j++)
                    printf("%f ", A.elements[i*A.width + j]);
                printf("\n");
            }
            printf("\n");
            printf("Matrix Q=\n");
            for(int i=0;i<B.height; i++){
                for(int j=0;j<B.width;j++)
                    printf("%f ", B.elements[i*B.width + j]);
                printf("\n");
            }
            printf("\n");
            printf("\nThe P=A*Q=\n");
            for(int i=0;i<C.height; i++){
                for(int j=0;j<C.width; j++)
                    printf("%f ", C.elements[i*C.width + j]);
                printf("\n");
            }
            printf("\n");*/
        MatElMul(C,A,D,rank);/*
    printf("\nThe P.*A=\n");
    for(int i=0;i<D.height; i++){
        for(int j=0;j<D.width; j++)
            printf("%f ", D.elements[i*D.width + j]);
        printf("\n");
    }*/
       // printf("\n");
        n = A.height; // number of matrix rows/cols
        m = A.width;
 get_host(A,D,Host_c,n,m,rank);

///filter the values more than average

        m = string_arg_m;
        rescale_to_standard(Host_c,m);
        bool* present;
        present = (bool*)malloc(m*sizeof(bool));
        double avg_F;
        double sum_present_F=0;
        avg_F = average_of(Host_c , m) ;
        int num_present_F = 0;
        for(i=0; i<m; i++)
        {
            if(Host_c[i]>=avg_F)
            {
                present[i]=true;
                num_present_F = num_present_F + 1;
            }
            else
            {
                present[i]=false;
            }
        }

        pack* packs;
        packs = (struct pack*)malloc(num_present_F*sizeof(pack));
        j=0;
        for(i=0; i<m; i++)
        {if(present[i])
            {
                packs[j].index = i;
                packs[j].value = Host_c[i];
                                sum_present_F = sum_present_F + i+1;

                j++;
            }
        }
           qsort(packs, num_present_F, sizeof(pack), comparator);

        int* count;
        count = (int*)malloc(m*sizeof(int));




        for(i=0; i<m; i++)
        {
            count[i]=0;
        }
                int sum_temp = 0;

        for(i=0; i<num_present_F; i++)
        {
          
                temp=(((packs[i].index+1)*pool_size)/sum_present_F);
                temp = temp+0.5;
                count[packs[i].index] = (int)temp;
          
        }
        /*for(i=0; i<m; i++)
        {
            range[i] = count[i] + sum_temp;
                sum_temp = range[i];
        }*/
        i=0;
        int pool_count = 0;
        while(i<m)
        {
           
                for(j=0; (j<count[i]&&pool_count<pool_size); j++)
                {int k;
                    for(k=0; k<A.width; k++)
                    {
                        Pool[pool_count][k]=A.elements[i*A.width+k];
                    }
                    pool_count++;
                }
          
            range[i] = sum_temp + j;
                sum_temp = range[i];
    i++;
        }
        int max_index =0;
        int second_max_index=0;
        bool initial = true;
        for(i=0; i<m; i++)
        {
            if(present[i])
            {
                if(initial)
                {
                    initial = false;
                    max_index = i;
                    second_max_index = i;
                }
                else
                {

                    if(Host_c[i]>Host_c[max_index])
                    {
                        second_max_index = max_index;
                        max_index = i;
                    }
                    else
                    {
                        if(Host_c[i]>Host_c[second_max_index])
                        {
                            second_max_index = i;
                        }
                    }

                }
            }
        }
        double* max_alpha_A ;
        max_alpha_A = (double*)malloc(A.width*sizeof(double));
        double* second_max_alpha_A ;
        second_max_alpha_A = (double*)malloc(A.width*sizeof(double));
       // printf("max_fitness:\t%f\n",Host_c[max_index]);
        get_row(A.elements,A.width,max_alpha_A,max_index);
        get_row(A.elements,A.width,second_max_alpha_A,second_max_index);
for(i=0;i<A.width;i++){
A.elements[i]=max_alpha_A[i];
}
for(i=0;i<A.width;i++){
A.elements[A.width+i]=second_max_alpha_A[i];
}
///code for doing crossover using pool
//print_accuracy(max_alpha_A,y->elements,x->elements,v,v_size,(np+nn),d);
//mat_accuracy(x,v_t,max_alpha_A,y->elements, acc,rank);
free(max_alpha_A);

free(second_max_alpha_A);

        pthread_t threads[NUM_THREADS];
           pthread_attr_t attr;
        size_t stacksize;
        int rc;
        long t;
        void *status;

        /* Initialize and set thread detached attribute */
       pthread_attr_init(&attr);
        pthread_attr_getstacksize (&attr, &stacksize);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

long temp_t;
        for(t=0; t<NUM_THREADS; t++)
        {temp_t = t;
            //printf("In main: creating thread %ld\n", temp_t);
            rc = pthread_create(&threads[temp_t], &attr, parallel_crossover, (void *)temp_t);
            if (rc)
            {
                printf("ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }

            }
               pthread_attr_destroy(&attr);
          for(t=0; t<NUM_THREADS; t++) {
      rc = pthread_join(threads[t], &status);
      if (rc) {
         printf("ERROR; return code from pthread_join() is %d\n", rc);
         exit(-1);
         }
      //printf("Main: completed join with thread %ld having a status of %ld\n",t,(long)status);
      }
    //printf("A.height%d\n",A.height );
   // printf("d is %d\n", d);
    for ( i = 0; i < A.width; ++i)
    {
        A.elements[(A.height-2)*A.width+i]=0;
    }
     for ( i = 0; i < A.width; ++i)
    {
        A.elements[(A.height-1)*A.width+i]=0;
    }
            fill_up_A(A.elements,(A.height-2),d);
            fill_up_A(A.elements,(A.height-1),d);
free(present);

free(count);
        

        };
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        
             printf( "%f\n", cpu_time_used);
                        MPI_Recv ( &alpha,1, alpha_type, 0, rank, MPI_COMM_WORLD,  &stat );
        for(i=0;i<A.width;i++){
alpha.values[i] = A.elements[i];
        }
        MPI_Send ( &alpha,1, alpha_type, 0, rank, MPI_COMM_WORLD );
    //printf("sent from %d\n", rank);

        //fprintf(alpha_array,"\n");
    }

//pthread_exit(NULL);


    }
    free(range);

        
        free(C.elements);
        free(D.elements);
        free(Host_c);
        // printf("final alpha\n");
        
        free(A.elements);
   free(x->elements);
        free(x);
        free(y->elements);
        free(y);
        free(v->elements);
        free(v);
         free(v_t->elements);
        free(v_t);
        free(acc->elements);
        free(acc);
        free(Q.elements);
        
        for(i=0;i<pool_size;i++){
            free(Pool[i]);
        }

        free(Pool);

   // fclose(fptr);
    //fclose(fit_acc);
   // fclose(alpha_array);
               
                //printf("here\n");   
}

else{
    int jj =0;
    A.width = (np+nn);
    A.height = 100 ;
      A.elements = (double*)malloc(A.width * A.height * sizeof(double));
      int row_num =0;
for(jj=0;jj<5;jj++){
for(i=1;i<numtasks;i++){
MPI_Send ( &alpha,1, alpha_type, i, i, MPI_COMM_WORLD ); 
  // printf("sent from master to%d\n", i);

            MPI_Recv ( &alpha,1, alpha_type, i, i, MPI_COMM_WORLD,  &stat );

     // printf("recved from %d to master\n", i);
            int ii =0;
for(ii = 0 ;ii < A.width ; ii ++){
    A.elements[row_num*A.width + ii] = alpha [ii] ;
}
row_num++;
}
}


}


  MPI_Finalize();
          return 0;
    }

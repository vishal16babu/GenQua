#include <stdio.h>
#include <stdlib.h>
typedef struct Matrix
{
    int width;
    int height;
    double* elements;
} Matrix;

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
void MatMul(Matrix AA, Matrix BB, Matrix CC);
 void MatElMul(Matrix AA, Matrix BB, Matrix CC);

int main(int argc, char* argv[]){
	Matrix A;
	Matrix B;
	Matrix C;
	 int m=2;
    int n=2;
    int l=2;
    int i;
    int j;
         
    A.height = m;
    A.width = n;
    A.elements = (double*)malloc(A.height*A.width*sizeof(double));
    B.height = n;
    B.width = l;
    B.elements = (double*)malloc(B.height*B.width*sizeof(double));
    time_t t1;
    srand(time(NULL));
    for( i=0; i<(A.height*A.width); i++)
    {
        A.elements[i] = (double)(rand()%10) ;
    }
    for( i=0; i<(A.height*A.width); i++)
    {
        B.elements[i] = (double)(rand()%10) ;
    }
    print_matrix(A);
    print_matrix(B);
    C.height = A.height ;
    C.width = B.width ;
    C.elements = (double*)malloc(C.height*C.width*sizeof(double)); 
MatElMul(A,B,C);
print_matrix(C);
return 0;
}
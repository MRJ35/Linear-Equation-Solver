#include<stdio.h>
#include<conio.h>

void printMatrix(int r,int c,double** a,double* b)	//Prints the augmented matrix
{
	int i,j;
	printf("\n\nMatrix is \n");
	for (i = 0; i < r; ++i)
	{
		for (j = 0; j < c; ++j)
		{
			printf("%f  ",a[i][j] );
		}
		printf("%f\n",b[i]);
	}
	printf("\n\n");
}
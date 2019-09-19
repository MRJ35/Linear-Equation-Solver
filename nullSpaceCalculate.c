#include<stdio.h>
#include<stdlib.h>


void rref(int row,int col,double** a,int* columnSwaps);

int nullSpaceFinder(int row,int col,double** a)
{
	int* columnSwaps; 
	columnSwaps = (int*)calloc(col,sizeof(int));	//Array to keep track of columns swapped during rref
	int i;
	for (i = 0; i < col; ++i)
	{
		columnSwaps[i] = -1;
	}

	
	rref(row,col,a,columnSwaps);
	
	free(a);
	free(columnSwaps);
	return 0;
}
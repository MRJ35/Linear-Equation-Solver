#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include<math.h>
int n,pt=-1;

int compare(double x,double y);
void nullSpaceMatrix(int row,int col,double** a,double* F);
void rref(int row,int col,double** a,int* columnSwaps);
int checkColumnExchange(int pos,int row,int col,double** a,int* columnSwaps);
int checkRow(int pos,int row,int col,double** a,int* columnSwaps);
void swaprow(int r1,int r2,int row,int col,double** a,int* columnSwaps);
int shiftRow(int pos,int row,int col,double** a,int* columnSwaps);
void columnExchange(int col1,int col2,int row,int col,double** a,int* columnSwaps);
void extractF(int pos,int row,int col,double** a,int* columnSwaps);
void printMatrixx(int r,int c,double** a);
int isZero(int row,int col,double** a);


int isZero(int row,int col,double** a)	//Checks if matrix is zero
{
	int cc=0,i,j;
	for (i = 0; i < row; ++i)
	{
		for (j = 0; j < col; ++j)
		{
			if(compare(a[i][j],0)) 
			{
				cc++;
			}
		}
	}
	if(cc==(row*col))
		return 1;
	else 
		return 0;
}


void printMatrixx(int r,int c,double** a)	//prints rref form of matrix
{
	printf("\nMatrix is \n");
	int i,j;
	for (i = 0; i < r; ++i)
	{
		for (j = 0; j < c; ++j)
		{
			printf("%f  ",a[i][j] );
		}
		printf("\n");
	}
	printf("\n");
}

int checkColumnExchange(int pos,int row,int col,double** a,int* columnSwaps)	//checks if column exchange is possible
{
	int i;
	for (i = pos+1; i < col; ++i)
	{
		double tmp_pivot = a[pos][i];
		if(!compare(tmp_pivot,0)) //nonzero
		{
			columnExchange(pos,i,row,col,a,columnSwaps);
			return 1;
		}
	}
	int flag = shiftRow(pos,row,col,a,columnSwaps);
	return flag;
}

int checkRow(int pos,int row,int col,double** a,int* columnSwaps) //checks if row is zero
{
	int c=0,i;
	for (i = 0; i < col; ++i)
	{
		if(compare(a[pos][i],0)) //zero
		{
			c++;
		}
	}
	if(c==col)
		return 1;
	else
		return 0;
}

void swaprow(int r1,int r2,int row,int col,double** a,int* columnSwaps)  //swaps row r1 and r2
{
	double* tmp;
	tmp = (double*)calloc(col,sizeof(double));
	int i;
	for (i = 0; i < col; ++i)
	{
		tmp[i] = a[r1][i];
	}
	for (i = 0; i < col; ++i)
	{
		a[r1][i] = a[r2][i];
	}
	for (i = 0; i < col; ++i)
	{
		a[r2][i] = tmp[i];
	}
	free(tmp);

}

int shiftRow(int pos,int row,int col,double** a,int* columnSwaps)  //shifts all zero rows to the bottom of matrix
{
	int i;
	for (i = pos+1; i < row; ++i)
	{
		if(!checkRow(i,row,col,a,columnSwaps))
		{
			swaprow(pos,i,row,col,a,columnSwaps);
			return 0;
		}
	}
	return -1;

}

void columnExchange(int col1,int col2,int row,int col,double** a,int* columnSwaps) //exchange col1,col2
{
	double* temp;
	temp = (double*)calloc(row,sizeof(double));
	int i;
	for (i = 0; i < row; ++i)
	{
		temp[i] = a[i][col1];
	}
	for (i = 0; i < row; ++i)
	{
		a[i][col1] = a[i][col2];
	}
	for (i = 0; i < row; ++i)
	{
		a[i][col2] = temp[i];
	}
	//printMatrixx(row,col,a);
	pt++;
	columnSwaps[pt] = col1;
	
	pt++;
	columnSwaps[pt] = col2;
	
	free(temp);
}


void extractF(int rank,int row,int col,double** a,int* columnSwaps)  //Forms F from rref form of matrix [I F]
{
	printMatrixx(row,col,a);
	int c = col-rank;
	double** F;
	F = (double**)calloc(col,sizeof(double));
	
	int i,j;
	for (i = 0; i < col; ++i)
	{
		F[i] = (double*)calloc(c,sizeof(double));
	}
	
	for (i = 0; i < col; ++i)
	{
		if(i<rank){
			for (j = 0; j < c; ++j)
			{
				F[i][j] = a[i][rank+j];
				if(F[i][j]!=0)
				{
					F[i][j] = -F[i][j];
					
				}
					
			}
		}
		else{
			for (j = 0; j < c; ++j)
			{
				if(j==(i-rank)){
					F[i][j] = 1;
				}
				else{
					F[i][j] = 0;
				}
			}
		}
	}
	
	if(columnSwaps[0]!=-1)	//columns were exchanged
	{
		printf("\nSwapping rows as per previous operations\n");
		while(pt>0)
		{
			swaprow(columnSwaps[pt],columnSwaps[pt-1],col,c,F,columnSwaps);
			pt = pt - 2;
		}
	}
	printf("\n---------------Nullspace--------------- \n");
	printMatrixx(col,c,F);
	free(F);
}

void rref(int row,int col,double** a,int* columnSwaps)		//forms rref form of matrix
{
	if(isZero(row,col,a))
	{
		printf("\n--------------------------------------\n");
		printf("	Nullspace belongs to whole R%d\n",col );
		printf("--------------------------------------\n");
		return;
	}
	int rank=0,i,j;
	int out = 0;
	int nullity=0;
	int maxPossiblePivots = 0;
	if(row<col)
		maxPossiblePivots = row;
	else
		maxPossiblePivots = col;
	for (i = 0; i < maxPossiblePivots; ++i)
	{
		double pivot = a[i][i];
		if(pivot==0)
		{
			int flag = checkColumnExchange(i,row,col,a,columnSwaps);
			if(flag==-1)	// exchange possible
			{
				out = i;
				break;
			}
			else{			// pivot exchanged
				i--;
			}
		}
		else{		//elimination to form rref form
			rank++;
			double div = a[i][i];
			for (j = 0; j < col; ++j)
			{
				if(a[i][j]!=0){
					a[i][j]/=div;
				}
			}
			for (j = 0; j < row; ++j)
			{
				if(j!=i){
				 double k = a[j][i]/a[i][i];
			 		for (int z = i; z < col; ++z)
			 		{
			 			a[j][z] = a[j][z] - a[i][z]*k;			 				 	
			 		}
			 
			 		//printMatrixx(row,col,a);
			 	}			 
			}
		}
	}
	printf("\nrank is %d\n",rank);
	nullity = col-rank;
	if(nullity==0)
	{
		printMatrixx(row,col,a);
		printf("\nNull Space consists of zero vector only\n");
	}
	else if(rank==maxPossiblePivots)
	{
		extractF(rank,row,col,a,columnSwaps);
	}
	else{
		
		if(out>0)
		{
			extractF(out,row,col,a,columnSwaps);
		}
	}
}



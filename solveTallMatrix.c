#include<stdio.h>
#include<stdlib.h>
void backSubstitution(int c,int r,double** a,double* b);
void findParticularSolution(int r,int c,double** a,double* b);
int compare(double x,double y);

int solveTallMatrix(int r,int c,double** a,double* b)
{
	int zeros=0,rows=0,j,i;
	for (j = r-1; j >= 0; --j)
	{	
		for (i = 0; i < c; ++i)  
		{
			if(compare(a[j][i],0))
			{
				zeros++;
			}
		}
		if(zeros==c && compare(b[j],0))   //Check if row is zero and equation is consistent 
		{
			rows++;
		}
		else if(zeros==c && !compare(b[j],0))	//Equation is inconsistent
		{
			return -1;
		}
	}
	if(r-rows==c)			//If number of pivots are equal to number of columns
	{
		backSubstitution(c,r,a,b);
		return 1;
	}
	else{					//Consistent equation but pivots are less
		printf("\n---------------------------------\n");
		printf("	Infinite Solutions!!\n");
		printf("---------------------------------\n");
		findParticularSolution(r,c,a,b);
		return 0;	
	}
	
}
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int compare(double x,double y);

int r_exchange(int r,int c,int r1,int c1,double** a,double* b){
	int flag=0,i,j;
	for(i=r1+1;i<r;i++)
	{
		double tt = a[i][c1];
		if(!compare(tt,0)) // Row below r1 have non zero value at position c1
		{
			double* temp;
			temp = (double*)calloc(c,sizeof(double));
			for(j=0;j<c;j++)
			{
				temp[j] = a[r1][j];
			}
			for(j=0;j<c;j++)
			{
				a[r1][j] = a[i][j];
			}
			for (j = 0; j < c; ++j)
			{
				a[i][j] = temp[j];
			}
			flag=1;
			float t = b[r1];
			b[r1] = b[i];
			b[i] = t;
			break;
		}
	}				
	if(flag==0)			//No row exchange possible
	{
		return c1;
	}
	else{				//Row exchanged successfully
		
		return -1;
	}
}
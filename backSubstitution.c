#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int compare(double x,double y);
void findParticularSolution(int r,int c,double** a,double* b);

void backSubstitution(int c,int r,double** a,double* b)
{

	double value[c],sum=0;
		int i,j;
		for(i=0;i<c;i++)
		{
			value[i] = 0;
		}
		int count=0;
		for(i=r-1;i>=0;i--)
		{
			if(compare(a[i][i],0) && r-count>=c && compare(b[i],0))	//Sufficient pivots uptil now
			{
				count++;
				i++;
			}
			else if(a[i][i]==0 && r-count<=c)	//Less pivots than number of columns
			{
				printf("\n---------------------------------\n");
				printf("	Infinite Solutions\n");
				printf("---------------------------------\n");
				findParticularSolution(r,c,a,b);
				return;
			}
			else{
					int inf=0;
					for(j=0;j<c;j++)
					{
						if(!compare(a[i][j],0))
						{
							inf++;
						}
					}
					if(inf>(r-i))	//Number of nonzero values in one equation are greater than r-i
					{
						printf("\n---------------------------------\n");
						printf("	Infinite Solutions\n");
						printf("---------------------------------\n");
						findParticularSolution(r,c,a,b);
						return;
					}
					sum=0;
					int zz = i;
					int ff = c-1;
					value[i] = a[i][ff];
					for(j=r-1;j>zz;j--){
						sum+=value[j]*a[i][j];    
					}
					value[zz] = (b[zz]-sum)/a[i][i];
			}
		}
		printf("\n---------------------------------\n");
		printf("	Unique Solution Exists\n");
		printf("---------------------------------\n");
		int alpha = 1;
		printf("-----------------------------------------\n");
		for(i=0;i<c;i++)
		{
			printf("x%d = %f\n",alpha,value[i]);
			alpha++;
		}
		printf("-----------------------------------------\n");
}
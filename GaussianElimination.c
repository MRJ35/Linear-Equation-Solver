#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int compare(double x,double y);
void findParticularSolution(int row,int col,double** a,double* b);
int checkPossibleColumns(int r1,int c1,int row,int col,double** a,double* b);
int r_exchange(int r,int c,int r1,int c1,double** a,double* b);
void printMatrix(int r,int c,double** a,double* b);
void solve(int r,int c,double** a,double* b);
int solveTallMatrix(int r,int c,double** a,double* b);
void backSubstitution(int c,int r,double** a,double* b);
int nullSpaceFinder(int row,int col,double** a);


int compare(double x,double y){       // compares whether x=y precisely
    double epsilon = 1e-8;
    if(fabs(x-y)<epsilon) return 1;
    else return 0;
}

void findParticularSolution(int r,int c,double** a,double* b)     // finds particular solution of system of equations
{
	double value[c],sum=0;
	int i,j;
	for(i=r-1;i>=0;i--)
	{
		
		for(j=0;j<c;j++)
		{
			value[j] = 0;
		}
		sum=0;
		int zz = i;
		for(j=c-1;j>=0;j--){
			if(i!=j && !compare(a[i][j],0) && !compare(b[i],0))
			{
				value[j] = 1;
			}
			sum+=value[j]*a[i][j];    
		}
		value[zz] = (b[zz]-sum)/a[i][i];
	}	
			int alpha = 97;
			printf("\nParticular Solution is \n");
			printf("-----------------------------------------\n");
			for(i=0;i<c;i++)
			{
				printf("%c = %f\n",alpha,value[i]);
				alpha++;
			}
			printf("-----------------------------------------\n");
		
}


int checkPossibleColumns(int r1,int c1,int row,int col,double** a,double* b)  //checks whether there are any other columns such that its value for row1 is nonzero.(Possibility of pivot for FAT matrix)
{
	int i;
	for (i = c1+1; i < col; ++i)
	{
		if(a[r1][i]!=0)
		{
			return 1;
		}
	}
	return -1;
}


int matrixType(int row,int col,double** a,double* b)  //Returns type of matrix.
{
	if(row>col)
		return 1;		//Tall Matrix
	else if(row<col)
		return -1;		//Fat Matrix
	else
		return 0;		//Square Matrix
}

void GaussianElimination(int r,int c,double** a,double* b){ //Performs Gaussian Elimination
	int type = matrixType(r,c,a,b);

	if(type==1) //Tall Matrix
	{
		printf("\n--------Tall Matrix-------\n");
		int i,j,k,z;
		for (i = 0,j=0; j < c; ++i,++j)
		{
			if(!compare(a[i][j],0))  //Pivot is nonzero
			{
				for (k = i+1; k < r; ++k)  //Makes all elements below pivot element zero
				{

					double ratio = (!compare(a[k][j],0)) ? a[k][j]/a[i][j] : 0;
					for(z = 0;z<c;z++)
					{
						a[k][z] -= a[i][z]*ratio;
					}
					b[k] = b[k] - b[i]*ratio;
				}
				//printMatrix(r,c,a,b);
			}
			else{  //Pivot is zero
				int flag = r_exchange(r,c,i,j,a,b); //Check if row exchange is possible

				if(flag>=0 && compare(b[i],0))	//Row exchange not possible and consistent equation  
				{
					printMatrix(r,c,a,b);
					printf("\n------------------------------\n");
					printf("	Infinite Solutions \n");
					printf("------------------------------\n");
					findParticularSolution(r,c,a,b);
					return;
				} 							//Inconsistent equation
				else if(flag>=0){
					printMatrix(r,c,a,b);
					printf("\n------------------------------\n");
					printf("	NO Solutions\n");
					printf("------------------------------\n");
					return;
				}
				else{  				//Rows exchanged successfully.
					i--;
					j--;
					//printMatrix(r,c,a,b);
				}
			}
		}
		printMatrix(r,c,a,b);
		int flag = solveTallMatrix(r,c,a,b);  
		if(flag==-1)	
		{
			printMatrix(r,c,a,b);
			printf("\n-----------------------------------\n");
			printf("	No solutions exists\n");
			printf("-----------------------------------\n");
			return;
		}
	}
	else if(type==-1)	//Fat Matrix
	{
		printf("\n-------FAT matrix-------\n");
		//printMatrix(r,c,a,b);
		int i,j;
		for (i = 0,j=0; i < r; ++i,++j)
		{
			double pivotValue = a[i][j];
			if(compare(a[i][j],0))	//Pivot is zero
			{
				int flag = r_exchange(r,c,i,j,a,b); //Check for row exchange
			
				if(flag>=0)
				{
					
					if(j!=(c-1))	//Row exchange not possible but next column is nonzero
					{
						i--;
					}
					else{
						printMatrix(r,c,a,b);
						if(b[i]==0){		//Consistent equation but pivot is zero
							printf("\n-----------------------------------\n");
							printf("	Infinite Solutions possible\n");
							printf("-----------------------------------\n");
							findParticularSolution(r,c,a,b);
							return;
						}
						else{			//Inconsistent equation
							printf("\n-----------------------------------\n");
							printf("	No solutions possible\n");
							printf("-----------------------------------\n");
							return;
						}
					}
				}
					//printMatrix(r,c,a,b);		
			}
			else{							//Nonzero pivot
				int z;
				for (z = 0; z < c; ++z)
				{
					if(!compare(a[i][z],0))
					{
						a[i][z]/=pivotValue;
					}
				}
				if(!compare(b[i],0)){ 
					b[i]/=pivotValue;
				}
				//printMatrix(r,c,a,b);
				int k;
				for (k = i+1; k < r; ++k)	//Make all elements below the pivot as 0
				{
			 		double ratio = (!compare(a[k][j],0)) ? (double)a[k][j]/a[i][j] : 0;
			 		for (z = 0; z < c; ++z)
			 		{
			 			if(!compare(a[i][z],0) && !compare(ratio,0)){
			 				a[k][z] -= a[i][z]*ratio;
			 			}
			 			if(compare(a[k][z],0))
			 			{
			 				a[k][z] = 0;
			 			}			 				 	
			 		}
			 		b[k] -= b[i]*ratio;
			 		if(compare(b[k],0))
			 			b[k] = 0;

			 		//printMatrix(r,c,a,b);			 
				}
				
			}

		}
		printMatrix(r,c,a,b);
		backSubstitution(c,r,a,b);
	}
	else{  //for square matrix 
		int i,k,z,j;
		for (i = 0,j=0; i < r; ++i,++j)        
		{
			if(!compare(a[i][j],0)){	//Nonzero pivot
				for (k = i+1; k < r; ++k)		//Elimination
				{
			 		double ratio = (double)a[k][j]/a[i][j];
			 		for (z = j; z < c; ++z)
			 		{
			 			if(!compare(a[i][z],0) && !compare(ratio,0)){
			 				a[k][z] = a[k][z] - a[i][z]*ratio;
			 			}			 				 	
			 		}
			 		if(!compare(b[i],0) && !compare(ratio,0)){
			 			b[k] = b[k] - b[i]*ratio;
			 		}
			 		//printMatrix(r,c,a,b);			 
				}
			}
			else{			//Pivot is zero
				int flag = r_exchange(r,c,i,j,a,b);  //Check if row exchange is possible
				if(flag>=0 && compare(b[i],0))	//Consistent equation but pivot is zero
				{
					printMatrix(r,c,a,b);
					printf("\n-------------------------------\n");
					printf("	Infinte Solutions \n");
					printf("-------------------------------\n");
					findParticularSolution(r,c,a,b);
					return;					
				}
				else if(flag>=0)			//Inconsistent Equations
				{
					printf("\n-------------------------------\n");
					printf("	NO Solutions \n");
					printf("-------------------------------\n");
					return;	
				}
				else{			//Row exchanged successfully
					i--;
					j--;
				}
				
			}
		}
		printMatrix(r,c,a,b);
		backSubstitution(c,r,a,b);
	}
	
	
}



int main()
{
	printf("\n\n--------------------Linear Equation Solver--------------------\n\n");
	printf("\nEnter rows and columns of matrix A : ");
	int r,c;
	scanf("%d %d",&r,&c);
	printf("Enter %dx%d augmented matrix A \n", r,c+1);
	double** a;
	a = (double**)calloc(r,sizeof(double));
	double** nullSpaceMatrix;
	nullSpaceMatrix = (double**)calloc(r,sizeof(double));
	int i;
	for (i = 0; i < r; ++i)
	{
		a[i] = (double*)calloc(c,sizeof(double));
		nullSpaceMatrix[i] = (double*)calloc(c,sizeof(double));
	}
	double *b;
	b = (double*)calloc(r,sizeof(double));
	int j;
	for (i = 0; i < r; ++i)
	{
		for ( j = 0; j < c; ++j)
		{
			scanf("%lf",&a[i][j]);
			nullSpaceMatrix[i][j] = a[i][j];
		}
		scanf("%lf",&b[i]);
	}
	GaussianElimination(r,c,a,b);

	printf("\n\n------------Now calculating nullspace----------\n");
	nullSpaceFinder(r,c,nullSpaceMatrix);
	free(a);
	free(b);
	free(nullSpaceMatrix);
	return 0;
}
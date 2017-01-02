/*Gaussian Elimination using openmp without pivoting and using back substitution*/ 
/*Refer readme.txt file for compiling and creating object file*/ 
/*Solving Equation of Xy=Z type where X and y are input solution and Z is output solution*/


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>//import this lib to perform openmp implementation

void initializeMat(int N, double** X, double* y);
void backSubstitution(int N, double** X, double* y, double* Z, int n_threads);
void *evaluate(void *);
void displayMat(int N, double** X, double* y);
void printAnswer(double* Z, int N);
void gaussian_pthread();

/* Structure of each thread*/
struct Thread_Data {
	double** X;
	double* y;
	int N;
	int numThreads;
} t_info;

/* Structure of Arguments on which each thread will work*/
struct Func_Arg
{
	int row_index;//index to be worked on by thread
	int t_index;//thread index
}f_arg;

/* Prints the solution vector Z of N size*/
void printAnswer(double* Z,int N)
{
	int i;
	printf("\nSolution Vector (x):\n\n");
	for (i=0;i<N;i++)
	{
		printf("|%lf|\n", Z[i]);
	}
	
}
/*Initialize X and y matrix with random values and scaling by 50000*/
void initializeMat(int N,double** X,double* y)
{
	int i, j;
	for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++)
		{
						
			X[i][j]=rand()/50000.0;
			
			//scanf("%lf",&X[i][j]);
				
		}
	y[i] = rand()/50000.0;

	}
}

/*Display X and y matrix*/
void displayMat(int N, double** X, double* y)
{	
	int i,j;
	printf("Displaying Initial Matrix.\n");
	for (i=0;i<N;i++)
	{
		printf("| ");
		for(j=0;j<N; j++)
		{	
			printf("%lf ",X[i][j]);
		}
		printf("| | %lf |\n",y[i]);
	}
}

/*Performing back substitution on matrix Z to computes solution values*/
void backSubstitution(int N, double** X, double* y, double* Z, int n_Threads)
{
	int i,j;
	for (i=N-1;i>=0;i--)
	{
		int count=0;
		Z[i] = y[i];
		for (j=i+1;j<N;j++)
		{
			Z[i]-=X[i][j]*Z[j];
			
		}
		Z[i] = Z[i]/X[i][i];
	}
}
/* this routine performs openmp implementation with n number of threads*/
void gaussian_openmp(int N, double** X,double* y,int n_threads)
{
	int i,j,k,mp;
	for (i=0;i<N-1;i++)
	{	
		/*performing loop in parallel with n number of threads. N, X,y,i variables are shared and j, k, mp variables are private to each thread*/
		#pragma omp parallel default(none) num_threads(n_threads) shared(N,X,y,i) private(j,k,mp)
		{
			/* setting scheduling dynamic for execution*/
			#pragma omp for schedule(dynamic)
			for(j=i+1;j<N;j++)
			{
				mp=X[j][i]/X[i][i];
				for(k=i;k<N;k++)
				{
					X[j][k]-=X[i][k]*mp;
				}	
				y[j]-=y[i]*mp;
			}
		}	
	}

}
/*
void serial_exec(int N, double** X,double* y)
{
	int i,j,k,mp;
	for (i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			mp=X[j][i]/X[i][i];
			for(k=i;k<N;k++)
			{
				X[j][k]-=X[i][k]*mp;
			}	
			y[j]-=y[i]*mp;
		}
	}
}
*/
void main()
{
	int i,j,k,n_thread,N,q;
	
	printf("Enter Number of Variables:");
	scanf("%d",&N);
	t_info.N=N;//Order of Matrix

	printf("Enter Number of Threads:");
	scanf("%d",&n_thread);
	t_info.numThreads=n_thread;//No of threads
	
	struct timeval start_time, end_elimination, end_substitution;
	double  substitution_time, total_time,elimination_time;
	
	/*Create X y Z matrixes dynamically*/
	double **X = (double **)calloc(N,sizeof(double*));
	for (q=0; q < N; q++)
	X[q] = (double*)calloc(N,sizeof(double*));
	
	double* y = (double*) malloc(sizeof(double)*N);
	double* Z = (double*) malloc(sizeof(double)*N);

	initializeMat(N,X,y);
	/*updating structure of thread*/
	t_info.X=X;
	t_info.y=y;
	
	printf("Openmp Execution:\n");
	gettimeofday(&start_time, NULL);//starting the clock time.
	gaussian_openmp(N,X,y,n_thread);//calling openmp function
	gettimeofday(&end_elimination, NULL);//terminating the clock time for gaussian elimiantion.
	backSubstitution(N,X,y,Z,n_thread);	
	gettimeofday(&end_substitution, NULL);//terminating the clock time for back substitution.
	elimination_time = ((end_elimination.tv_sec - start_time.tv_sec) * 1000000u + end_elimination.tv_usec - start_time.tv_usec) / 1.e6;
	substitution_time = ((end_substitution.tv_sec - end_elimination.tv_sec) * 1000000u + end_substitution.tv_usec - end_elimination.tv_usec) / 1.e6;
	total_time = ((end_substitution.tv_sec - start_time.tv_sec) * 1000000u + end_substitution.tv_usec - start_time.tv_usec) / 1.e6;

	printAnswer(Z,N);
	
	printf("Substitution execution time: %.3f seconds.\n", elimination_time);
	printf("Substitution execution time: %.3f seconds.\n", substitution_time);
	printf("Total execution: \n%.3f seconds elapsed with %d threads used.\n\n", total_time, n_thread);

	
			

}

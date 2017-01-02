/*Gaussian Elimination using pthread without pivoting and using back substitution*/ 
/*Refer readme.txt file for compiling and creating object file*/ 
/*Solving Equation of Xy=Z type where X and y are input solution and Z is output solution*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h>//import this lib to create pthreads for implementation

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

/* this routine is called by each thread to work on independent index.This is function where pthread efficiency is acheived*/
/* each thread works on same row because there is dependencies between rows of matrix*/
void gaussian_pthread()
{
	pthread_t threads[t_info.numThreads];//initialize pthread variables
	
	int i,j,k;
	float mp;
	
	
	for(i=0;i<t_info.N-1;i++)
	{
		
		for (k=0;k<t_info.numThreads;k++) 
		{	
			struct Func_Arg *f_arg;
			f_arg = (struct Func_Arg *) malloc(sizeof(struct Func_Arg));
			f_arg->t_index=k;
			f_arg->row_index=i;
			pthread_create(&threads[k],NULL,evaluate,(void *)f_arg);//creates thread with unique id and calls evaluate routine with unique index
		}
		for(j=0;j<t_info.numThreads;j++)
		{
			pthread_join(threads[j],NULL);//join the threads to update each row.
		}
	}	
	
		
}

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

/*This routine is called by each thread with unique index to perform computation on same row*/
void *evaluate(void *f_arg) 
{
		int norm, row, col;
		float mp;
		struct Func_Arg *n;
		n = (struct Func_Arg *)f_arg;
		norm=n->row_index;//row index on which each thread performs computation.
		for (row = norm+(n->t_index)+1; row < t_info.N; row+=t_info.numThreads) 
		{	
			mp=t_info.X[row][norm]/t_info.X[norm][norm];
			for (col=norm;col< t_info.N;col++) 
			{
				t_info.X[row][col]-=t_info.X[norm][col]*mp;
			}
			t_info.y[row]-= t_info.y[norm]*mp;
		}
		pthread_exit(NULL);
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
	
	printf("Pthread Execution:\n");
	gettimeofday(&start_time, NULL);//starting the clock time.
	gaussian_pthread();//calling pthread function
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

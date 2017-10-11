#include<stdio.h>
#include<omp.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define NOW()(omp_get_wtime())
typedef double mytime_t;

double get_min (const double *x, size_t n)
{
	double min_val = 0.0;
	int i;
	
	#pragma omp parallel for reduction(min:min_val)
	for (i=0; i<n;i++)
	{
		if (x[i]<min_val)
		{
			min_val = x[i];
		}
	}

	return min_val;
}

double get_max (const double *x, size_t n)
{
	double max_val = 0.0;
	int i;
	
	#pragma omp parallel for reduction(max:max_val)
	for (i=0;i<n;i++)
	{
		if (x[i]>max_val)
		{
			max_val = x[i];
		}
	}

	return max_val;
}

void histogram (size_t n, const double *x, int num_bins, double **bin_bdry, size_t **bin_count)
{
	double min = get_min (x, n);
	double max = get_max (x, n);
	double interval = (max-min)/num_bins;

	//printf("1b min %lf max %lf interval %lf\n", min, max, interval);
	int len_bin_bdry = num_bins+1;

	double  *bb1 = new double[len_bin_bdry];
	//size_t  *bc1 = new size_t[num_bins];
	size_t  *bc1 = new size_t[num_bins]();
	
	

	//*bin_bdry = (double *)malloc(sizeof(double) * len_bin_bdry);
	bb1[0] = min;
	//bb1[len_bin_bdry-1] = (max+1) * 1e-14;	
	bb1[len_bin_bdry-1] = max;	
	//*bin_count = new size_t[num_bins];
	int k;
	int l;

#pragma omp parallel
	{
		size_t *bc1_loc = new size_t[num_bins]();
		#pragma omp for 
		for (k=1; k<len_bin_bdry; k++)
		{
			bb1[k] = (min) + (interval*(k));
		}
		
		#pragma omp barrier
	  
		#pragma omp for
		for	(l=0; l<n; l++)
		{
			double bin_index = (x[l] - min)/interval;
		
	
			int bi = (int)bin_index;
			if (bi == num_bins)
			{
				bc1_loc[bi-1]++;
				continue;
			}
			bc1_loc[bi]++;
		}


		#pragma omp barrier
        #pragma omp critical
		for (l=0;l <num_bins; l++)
		{
			bc1[l] += bc1_loc[l];
		}
	}
#if 0
	for(l=0;l<num_bins;l++)
	{
		printf("BC %ld\n", bc1[l]);
	}
#endif
	*bin_bdry = bb1;
	*bin_count = bc1;
//}
	}
int main(int argc, char *argv[])
{
	double outlier;
	int outlier_required = atoi(argv[4]);
	

		int j =0;
		int i=0;
		FILE *fp;
		char *line = NULL;
		size_t len_word = 0;
		ssize_t read;
		size_t num_of_elements;
		size_t total_num_of_elements;
		int nb;
			
		num_of_elements = atol(argv[1]);
		nb = atoi(argv[3]);

		total_num_of_elements = num_of_elements;
		if (outlier_required == 1)
		{
			outlier = nb;
			total_num_of_elements += 1;
		}
		

		double *x = new double[total_num_of_elements];
		fp = fopen(argv[2], "r");
		if (fp == NULL)
			exit(0);


		while(((read = getline(&line, &len_word, fp)) != -1) && (i<num_of_elements))
		{
			x[i++] = atof(line);
		}	
		
	
		if (outlier_required == 1)
		{
			//x[num_of_elements] = outlier;
			x[i] = outlier;
		}	

		fclose(fp);
		if (line)
			free(line);


#if 0	
	double min_num = 0;
	double max_num = 1;

	srand(5);

	for (size_t i=0; i<n; i++)
	{
		double f = (double)rand()/RAND_MAX;
		x[i] = min_num + f * (max_num - min_num);
		printf("%lf\n", x[i]);
	}
#endif
#if 0
	for (size_t i=0; i<total_num_of_elements; i++)
	{
		printf("%ld %lf\n", i, x[i]);
	}

#endif

	double *bb = NULL;
	size_t *bc = NULL;

	//printf("TOTAL NUM %ld\n", total_num_of_elements);
	mytime_t start = NOW();
	histogram (total_num_of_elements, x, nb, &bb, &bc);
	mytime_t end = NOW();
	
	printf("ELAPSED TIME %lf\n", end - start);

#if 0
	for (j=0; j<nb; j++)
	{
		printf("BC %ld\n", bc[j]);
	}

	for (j=0; j<nb+1; j++)
	{
		printf("BB %lf\n", bb[j]);
	}	
#endif

}


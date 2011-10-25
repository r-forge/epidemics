#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h> // for matrix & vectors operations
#include <gsl/gsl_permutation.h> //for LU decomposition
#include <gsl/gsl_sf_gamma.h> // for combinations calculation
#include <gsl/gsl_eigen.h>

#include <time.h>

#include "common.h"

/******************************************************************************/
/* DECLARATION OF GLOBAL VARIABLES                                            */
/******************************************************************************/

gsl_rng * rng; /* to be included as an external variable eslewhere */

/******************************************************************************/
/* MAIN                                                                       */
/******************************************************************************/

#define OPTIONS "awrh"

void readParameters(parameters * param) 
{
    int k,i,j;
    double a;
    int tau;
	FILE * paramInit, *tauInit, *KInit, *GTInit;
	char val[30];
	
	//////////* Parameters *//////////
	
	if ((paramInit=fopen("param.txt","r"))==NULL) 
	{
		printf("Cannot read param.txt");
		exit(2);
	}
	
	// NbGroups
    fscanf(paramInit,"%s %d",val,&(param->NbGroups));
    printf("%s %d\n",val,param->NbGroups);
	
	// T
	fscanf(paramInit,"%s %d",val,&(param->T));
    printf("%s %d\n",val,param->T);
    
    // S
	fscanf(paramInit,"%s %d",val,&(param->S));
    printf("%s %d\n",val,param->S);
    
    // p
	fscanf(paramInit,"%s %d",val,&(param->p));
    printf("%s %d\n",val,param->p);
    
    //////////* Changing dates *//////////   
    
    // Tau
    if(param->p>0)
    {
    	if ((tauInit=fopen("ChangingDates.txt","r"))==NULL) 
		{
			printf("Cannot read ChangingDates.txt");
			exit(2);
		}
	
    	param->tau = gsl_vector_alloc(param->p);
    	for (k=0 ; k<param->p ; k++)
		{
			fscanf(tauInit,"%d",&tau);
    		printf("tau[%d] %d\n",k,tau);
        	gsl_vector_set(param->tau,k,tau); 
    	}
    }
   
   	//////////* Next generation matrix *//////////   
    
    // K
    if ((KInit=fopen("NextGenerationMatrix.txt","r"))==NULL) 
	{
		printf("Cannot read NextGenerationMatrix.txt");
		exit(2);
	}
    for (i=0 ; i<param->p+1 ; i++)
    {
        (param->K)[i] = gsl_matrix_alloc(param->NbGroups,param->NbGroups);
    }
	for (k=0 ; k<param->p+1 ; k++)
	{
        for(i=0 ; i<param->NbGroups ; i++)
        {
             for(j=0 ; j<param->NbGroups ; j++)
             {
                fscanf(KInit,"%lg",&a);
    			printf("K_%d[%d,%d] %lg\n",k,i,j,a);
                gsl_matrix_set((param->K)[k],i,j,a);
             }
        }
	}
    
    //////////* Generation time distribution *//////////   
        		    	
    if ((GTInit=fopen("GTdistr.txt","r"))==NULL) 
	{
		printf("Cannot read GTdistr.txt");
		exit(2);
	}
  	param->GTdistr = gsl_matrix_alloc(param->NbGroups,param->S);
  	
  	for(i=0 ; i<param->NbGroups ; i++)
    {
        for(j=0 ; j<param->S ; j++)
        {
           	fscanf(GTInit,"%lg",&a);
    		printf("GT[%d,%d] %lg\n",i,j,a);
            gsl_matrix_set(param->GTdistr,i,j,a);
        }
    }
    
    fclose(paramInit);
    if(param->p>0){fclose(tauInit);}
    fclose(KInit);
    fclose(GTInit);
    
}
  		
void freeParam(parameters * param) 
{
	/************************************************************************/
	/*      Input : structure parameters                                    */
	/*      Output : strucutre parametres                                   */
    /*               - 1 gsl_vector freed - (p+2) matrices freed            */
	/************************************************************************/
	
	int i;
	
	gsl_matrix_free(param->GTdistr);
	
    gsl_vector_free(param->tau);
    
	for (i=0 ; i<param->p+1 ; i++)
    {
        gsl_matrix_free(param->K[i]);
    }
    
}

int CalcEpidemicSize(parameters *param,gsl_matrix *Incid)
{
	int i, j;
	int prov=0;
	for (i=0 ; i<param->NbGroups ; i++)
	{
		for (j=0 ; j<param->T ; j++)
		{
			prov+=gsl_matrix_get(Incid,i,j);
			
		}
	}	
	return(prov);
}


void SimulEpidPoiss(parameters *param, gsl_matrix *Incid, gsl_vector *Incid0)
{
	int k, t, g, i, s;
	int tMin, tMax;
	double lambda;
	double prov;
	int poiss;
	
	for (k=0 ; k<param->NbGroups ; k++)
	{
		gsl_matrix_set(Incid,k,0,gsl_vector_get(Incid0,k));	
	}

	for (i=0 ; i<param->p+1 ; i++)
	{
		if(i==0)
		{
			tMin=1;
		}else
		{
			tMin=gsl_vector_get(param->tau,i-1);
		}
		if(i==param->p)
		{
			tMax=param->T;
		}else
		{
			tMax=gsl_vector_get(param->tau,i);
		}
		for (t=tMin ; t<tMax ; t++)
		{
			//printf("t %d\n",t);
			//fflush(stdout);
			for(k=0 ; k<param->NbGroups ; k++)
			{
				//printf("k %d\n",k);
				//fflush(stdout);
				lambda=0;
				for(g=0 ; g<param->NbGroups ; g++)	
				{
					//printf("g %d\n",g);
					//fflush(stdout);
					prov=0;
					for(s=1 ; s<=GSL_MIN(t,param->S) ; s++)
					{
						
						prov+=gsl_matrix_get(Incid,g,t-s)*gsl_matrix_get(param->GTdistr,g,s-1);
						//printf("s %d %lg\n",s,prov);
						//fflush(stdout);
					}
					prov*=gsl_matrix_get(param->K[i],g,k);
					//printf("prov %lg\n",prov);
					//fflush(stdout);
        			lambda+=prov;
        			//printf("lambda %lg\n",lambda);
					//fflush(stdout);
        			
				}
				poiss=gsl_ran_poisson(rng,lambda);
				//printf("lambda %lg inci %d\n",lambda,poiss);
				//fflush(stdout);
				gsl_matrix_set(Incid,k,t,poiss);
			}
		}
			
	}
			
}


void SimulEpidPoissPerGroup(parameters *param, gsl_matrix *Incid, gsl_matrix *IncidPerGroup, gsl_vector *Incid0PerGroup)
{
	int k, t, g, i, s;
	int tMin, tMax;
	double lambda;
	double prov;
	gsl_vector * vectProv = gsl_vector_calloc(param->NbGroups);
	int poiss;

	for (k=0 ; k<param->NbGroups*param->NbGroups ; k++)
	{
		gsl_matrix_set(IncidPerGroup,k,0,gsl_vector_get(Incid0PerGroup,k));
	}
	for (k=0 ; k<param->NbGroups ; k++)
	{
		prov=0;
		for (g=0 ; g<param->NbGroups ; g++)
		{
			prov+=gsl_matrix_get(IncidPerGroup,g*param->NbGroups+k,0);
		}
		gsl_matrix_set(Incid,k,0,prov);
	}

	for (i=0 ; i<param->p+1 ; i++)
	{
		if(i==0)
		{
			tMin=1;
		}else
		{
			tMin=gsl_vector_get(param->tau,i-1);
		}
		if(i==param->p)
		{
			tMax=param->T;
		}else
		{
			tMax=gsl_vector_get(param->tau,i);
		}
		for (t=tMin ; t<tMax ; t++)
		{
			//printf("t %d\n",t);
			//fflush(stdout);
			for(k=0 ; k<param->NbGroups ; k++)
			{
				//printf("k %d\n",k);
				//fflush(stdout);
				lambda=0;
				for(g=0 ; g<param->NbGroups ; g++)
				{
					//printf("g %d\n",g);
					//fflush(stdout);
					prov=0;
					for(s=1 ; s<=GSL_MIN(t,param->S) ; s++)
					{

						prov+=gsl_matrix_get(Incid,g,t-s)*gsl_matrix_get(param->GTdistr,g,s-1);
						//printf("s %d %lg\n",s,prov);
						//fflush(stdout);
					}
					prov*=gsl_matrix_get(param->K[i],g,k);
					gsl_vector_set(vectProv,g,prov);
					//printf("prov %lg\n",prov);
					//fflush(stdout);
        			lambda+=prov;
        			//printf("lambda %lg\n",lambda);
					//fflush(stdout);

				}
				poiss=gsl_ran_poisson(rng,lambda);
				//printf("lambda %lg inci %d\n",lambda,poiss);
				//fflush(stdout);
				gsl_matrix_set(Incid,k,t,poiss);
				for(g=0 ; g<param->NbGroups ; g++)
				{
					gsl_matrix_set(IncidPerGroup,g*param->NbGroups+k,t,(double)poiss*gsl_vector_get(vectProv,g)/lambda);
				}
			}
		}

	}

}

void SimulEpidNegBin(double VarDivMean, parameters *param, gsl_matrix *Incid, gsl_vector *Incid0)
{
	// VarDivMean = Variance/Mean of the negative binomial considered
	int k, t, g, i, s;
	int tMin, tMax;
	double lambda;
	double prov;
	int poiss;
	
	for (k=0 ; k<param->NbGroups ; k++)
	{
		gsl_matrix_set(Incid,k,0,gsl_vector_get(Incid0,k));	
	}

	for (i=0 ; i<param->p+1 ; i++)
	{
		if(i==0)
		{
			tMin=1;
		}else
		{
			tMin=gsl_vector_get(param->tau,i-1);
		}
		if(i==param->p)
		{
			tMax=param->T;
		}else
		{
			tMax=gsl_vector_get(param->tau,i);
		}
		for (t=tMin ; t<tMax ; t++)
		{
			//printf("t %d\n",t);
			//fflush(stdout);
			for(k=0 ; k<param->NbGroups ; k++)
			{
				//printf("k %d\n",k);
				//fflush(stdout);
				lambda=0;
				for(g=0 ; g<param->NbGroups ; g++)	
				{
					//printf("g %d\n",g);
					//fflush(stdout);
					prov=0;
					for(s=1 ; s<=GSL_MIN(t,param->S) ; s++)
					{
						
						prov+=gsl_matrix_get(Incid,g,t-s)*gsl_matrix_get(param->GTdistr,g,s-1);
						//printf("s %d %lg\n",s,prov);
						//fflush(stdout);
					}
					prov*=gsl_matrix_get(param->K[i],g,k);
					//printf("prov %lg\n",prov);
					//fflush(stdout);
        			lambda+=prov;
        			//printf("lambda %lg\n",lambda);
					//fflush(stdout);
        			
				}
				poiss=gsl_ran_negative_binomial(rng,1/VarDivMean,lambda/(VarDivMean-1));
				//printf("lambda %lg inci %d\n",lambda,poiss);
				//fflush(stdout);
				gsl_matrix_set(Incid,k,t,poiss);
			}
		}
			
	}
			
}




int main(int argc, char *argv[]) 
{
	int i,j;
	int inc;
	double inc2;
	FILE *Incidence;
	FILE *IncidencePerGroup;
	int EpidemicSize=0;
	
	/*parameters */
	parameters param;
	
	/*****************************************************/
	/***     DEFINITION OF A RANDOM GENERATOR          ***/
	/*****************************************************/
		
    time_t t;
    t = time(NULL); // time in seconds, used to change the seed of the random generator
    const gsl_rng_type *typ;
    gsl_rng_env_setup();    
    typ=gsl_rng_default;
    rng=gsl_rng_alloc(typ);
    gsl_rng_set(rng,t); // changes the seed of the random generator

    /*****************************************************/

			
    /* reading simulation parameters */
   	readParameters(&param);
    
    /* opening and preparing output files */
    
    /* Simulation algorithm */
   	gsl_matrix *IncidPerGroup = gsl_matrix_calloc(param.NbGroups*param.NbGroups,param.T);
   	gsl_vector *Incid0PerGroup = gsl_vector_calloc(param.NbGroups*param.NbGroups);
   	gsl_matrix *Incid = gsl_matrix_calloc(param.NbGroups,param.T);

   	
   	FILE * IncidInit;
	if ((IncidInit=fopen("Incid0.txt","r"))==NULL) 
	{
		printf("Cannot open Incid0.txt");
		exit(2);
	}
	for (i=0 ; i<param.NbGroups*param.NbGroups ; i++)
	{
    	fscanf(IncidInit,"%d",&j);
    	printf("Incid0[%d] %d\n",i,j);
    	gsl_vector_set(Incid0PerGroup,i,j);
	}
   	
   	while(EpidemicSize<50*param.NbGroups*(param.p+1))//20)
   	{
   		//SimulEpidNegBin(4, &param, Incid, Incid0);
   		SimulEpidPoissPerGroup(&param, Incid, IncidPerGroup, Incid0PerGroup);
   		EpidemicSize=CalcEpidemicSize(&param,Incid);
   		/*printf("%u\n",EpidemicSize);
   		fflush(stdout);*/
   	}
   	
   printf("EpidemicSize %d\t",EpidemicSize);
   fflush(stdout);
   
   	
   	
   	if ((Incidence=fopen("Incidence.txt","w"))==NULL) 
	{
		printf("Cannot open Incidence.txt");
		exit(2);
	}
	for (i=0 ; i<param.NbGroups ; i++)
   	{
   		for (j=0 ; j<param.T ; j++)
   		{
   			inc=gsl_matrix_get(Incid,i,j);
   			//printf("%d\t",inc);
   			//fflush(stdout);
   			fprintf(Incidence,"%d\t",inc);
   			fflush(Incidence);	
   		}	
   		fprintf(Incidence,"\n");
   		fflush(Incidence);
   	}

	if ((IncidencePerGroup=fopen("IncidencePerGroup.txt","w"))==NULL)
		{
			printf("Cannot open IncidencePerGroup.txt");
			exit(2);
		}
		for (i=0 ; i<param.NbGroups*param.NbGroups ; i++)
	   	{
	   		for (j=0 ; j<param.T ; j++)
	   		{
	   			inc2=gsl_matrix_get(IncidPerGroup,i,j);
	   			//printf("%d\t",inc);
	   			//fflush(stdout);
	   			fprintf(IncidencePerGroup,"%lg\t",inc2);
	   			fflush(IncidencePerGroup);
	   		}
	   		fprintf(IncidencePerGroup,"\n");
	   		fflush(IncidencePerGroup);
	   	}
   
    /* Closing files and freeing memory */
    freeParam(&param);
    fclose(Incidence);
    
	//getchar();
    return 0;

}

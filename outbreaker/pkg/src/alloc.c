#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_eigen.h>

#include <omp.h>

#include "common.h"
#include "init.h"
#include "InputOutput.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "alloc.h"
#include "tuneVariances.h"

extern gsl_rng * rng;

/**************** nb_data ****************/

nb_data * createNbData()
{
	nb_data *nb = (nb_data *) malloc(sizeof(nb_data));

	nb->NbAdmissions = (int *) calloc(NbPatients, sizeof(int));
	nb->NbPosSwabs = (int *) calloc(NbPatients, sizeof(int));
	nb->NbNegSwabs = (int *) calloc(NbPatients, sizeof(int));

	return nb;
}

void freeNbData(nb_data *nb)
{
	free(nb->NbAdmissions);
	free(nb->NbPosSwabs);
	free(nb->NbNegSwabs);

	free(nb);
}
/******************************************/

/**************** raw_data ****************/

raw_data *createRawData(nb_data *nb)
{
	raw_data *data = (raw_data *) malloc(sizeof(raw_data));

	int i;

	data->ward = (int *) calloc(NbPatients, sizeof(int));
	data->PatientIndex = (int *) calloc(NbPatients, sizeof(int));
	data->timeSeq = (double *) calloc(NbPatients, sizeof(double));

	for(i=0 ; i<NbPatients ; i++)
	{
		data->A[i] = gsl_vector_calloc(nb->NbAdmissions[i]);
		data->D[i] = gsl_vector_calloc(nb->NbAdmissions[i]);
		data->P[i] = gsl_vector_calloc(nb->NbPosSwabs[i]);
		data->N[i] = gsl_vector_calloc(nb->NbNegSwabs[i]);
		data->IsInHosp[i] = gsl_vector_calloc(T);
	}

	return data;

}

void freeRawData(raw_data *data)
{
	int i;

	free(data->ward);
	free(data->timeSeq);

	for(i=0 ; i<NbPatients ; i++)
	{
		gsl_vector_free(data->A[i]);
		gsl_vector_free(data->D[i]);
		gsl_vector_free(data->P[i]);
		gsl_vector_free(data->N[i]);
		gsl_vector_free(data->IsInHosp[i]);
	}

	free(data);
}
/******************************************/

/**************** aug_data ****************/

aug_data *createAugData()
{
	aug_data *augData = (aug_data *) malloc(sizeof(aug_data));

	augData->C = (int *) calloc(NbPatients, sizeof(int));
	augData->E = (int *) calloc(NbPatients, sizeof(int));

	augData->I0 = (int *) calloc(NbPatients, sizeof(int));
	augData->I1 = (int *) calloc(NbPatients, sizeof(int));

	return augData;
}

void freeAugData(aug_data *augData)
{
	free(augData->C);
	free(augData->E);

	free(augData->I0);
	free(augData->I1);

	free(augData);
}

void copyAugData(aug_data *augDataDest, aug_data *augDataSource)
{
	memcpy(augDataDest->C, augDataSource->C, NbPatients*sizeof(double));
	memcpy(augDataDest->E, augDataSource->E, NbPatients*sizeof(double));
	memcpy(augDataDest->I0, augDataSource->I0, NbPatients*sizeof(double));
	memcpy(augDataDest->I1, augDataSource->I1, NbPatients*sizeof(double));
}
/******************************************/

/***************** param ******************/

parameters *createParam()
{
	parameters *param = (parameters *) malloc(sizeof(parameters));

	param->beta = gsl_matrix_calloc(2,2);

	return param;
}

void freeParam(parameters *param)
{
	gsl_matrix_free(param->beta);

	free(param);
}

void copyParam(parameters * paramDest, parameters * paramSource)
{
	gsl_matrix_memcpy (paramDest->beta,paramSource->beta);
	paramDest->betaWardOut=paramSource->betaWardOut;
	paramDest->betaOutOut=paramSource->betaOutOut;

	paramDest->Sp=paramSource->Sp;
	paramDest->Se=paramSource->Se;

	paramDest->Pi=paramSource->Pi;

	paramDest->mu=paramSource->mu;
	paramDest->sigma=paramSource->sigma;

	paramDest->nu1=paramSource->nu1;
	paramDest->nu2=paramSource->nu2;

	paramDest->tau=paramSource->tau;
	paramDest->alpha=paramSource->alpha;
}

/******************************************/

/************ MCMC internals **************/

mcmcInternals *createMcmcInternals()
{
	mcmcInternals *MCMCSettings = (mcmcInternals *) malloc(sizeof(mcmcInternals));

	MCMCSettings->Sigma_beta = gsl_matrix_calloc(2,2);

	return MCMCSettings;
}

void printStdProp(mcmcInternals *MCMCSettings)
{
	int i,j;

	for (i=0 ;  i<2 ; i++)
    {
        for (j=0 ; j<2 ; j++)
        {
        	printf("Std proposal for beta_%d,%d: %lg\n",i,j,gsl_matrix_get(MCMCSettings->Sigma_beta,i,j));
        }
    }
	printf("Std proposal for betaWardOut: %lg\n",MCMCSettings->Sigma_betaWardOut);
	printf("Std proposal for betaOutOut: %lg\n",MCMCSettings->Sigma_betaOutOut);
	printf("Std proposal for mu: %lg\n",MCMCSettings->Sigma_mu);
	printf("Std proposal for sigma: %lg\n",MCMCSettings->Sigma_sigma);
	printf("Std proposal for nu1: %lg\n",MCMCSettings->Sigma_nu1);
	printf("Std proposal for nu2: %lg\n",MCMCSettings->Sigma_nu2);
	printf("Std proposal for tau: %lg\n",MCMCSettings->Sigma_tau);
	printf("Std proposal for alpha: %lg\n",MCMCSettings->Sigma_alpha);

	fflush(stdout);
}

void freeMcmcInternals(mcmcInternals *MCMCSettings)
{
	gsl_matrix_free(MCMCSettings->Sigma_beta);

	free(MCMCSettings);
}

/***************************************/

/************ Acceptance ***************/

acceptance *createAcceptance()
{
	acceptance *accept = (acceptance *) malloc(sizeof(acceptance));

	accept->PourcAcc_beta = gsl_matrix_calloc(2,2);

	accept->PourcAcc_betaWardOut=0;
	accept->PourcAcc_betaOutOut=0;
	accept->PourcAcc_mu=0;
	accept->PourcAcc_sigma=0;
	accept->PourcAcc_nu1=0;
	accept->PourcAcc_nu2=0;
	accept->PourcAcc_tau=0;
	accept->PourcAcc_alpha=0;

	return accept;
}

void reInitiateAcceptance(acceptance *accept)
{
	int i,j;
	for(i=0 ; i<2 ; i++)
	{
		for(j=0 ; j<2 ; j++)
		{
			gsl_matrix_set(accept->PourcAcc_beta,i,j,0);
		}
	}

	accept->PourcAcc_betaWardOut=0;
	accept->PourcAcc_betaOutOut=0;
	accept->PourcAcc_mu=0;
	accept->PourcAcc_sigma=0;
	accept->PourcAcc_nu1=0;
	accept->PourcAcc_nu2=0;
	accept->PourcAcc_tau=0;
	accept->PourcAcc_alpha=0;

}

void printAcceptance(acceptance *accept, NbProposals *NbProp)
{
	int i,j;

	for(i=0 ; i<2 ; i++)
	{
		for(j=0 ; j<2 ; j++)
		{
			printf("Prob accept beta_%d,%d\t%lg\n",i,j,gsl_matrix_get(accept->PourcAcc_beta,i,j)/gsl_matrix_get(NbProp->NbProp_beta,i,j));
			fflush(stdout);
		}
	}

	printf("Prob accept betaWardOut\t%lg\n",accept->PourcAcc_betaWardOut/NbProp->NbProp_betaWardOut);
	fflush(stdout);
	printf("Prob accept betaOutOut\t%lg\n",accept->PourcAcc_betaOutOut/NbProp->NbProp_betaOutOut);
	fflush(stdout);
	printf("Prob accept mu\t%lg\n",accept->PourcAcc_mu/NbProp->NbProp_mu);
	fflush(stdout);
	printf("Prob accept sigma\t%lg\n",accept->PourcAcc_sigma/NbProp->NbProp_sigma);
	fflush(stdout);
	printf("Prob accept nu1\t%lg\n",accept->PourcAcc_nu1/NbProp->NbProp_nu1);
	fflush(stdout);
	printf("Prob accept nu2\t%lg\n",accept->PourcAcc_nu2/NbProp->NbProp_nu2);
	fflush(stdout);
	printf("Prob accept tau\t%lg\n",accept->PourcAcc_tau/NbProp->NbProp_tau);
	fflush(stdout);
	printf("Prob accept alpha\t%lg\n",accept->PourcAcc_alpha/NbProp->NbProp_alpha);
	fflush(stdout);

}

void freeAcceptance(acceptance *accept)
{
	gsl_matrix_free(accept->PourcAcc_beta);

	free(accept);
}

/********************************************/

/************* Is Acceptance OK *************/

isAcceptOK *createIsAcceptOK()
{
	isAcceptOK *acceptOK = (isAcceptOK *) malloc(sizeof(isAcceptOK));

	acceptOK->IsAccOK_beta = gsl_matrix_calloc(2,2);

	acceptOK->IsAccOK_betaWardOut=0;
	acceptOK->IsAccOK_betaOutOut=0;
	acceptOK->IsAccOK_mu=0;
	acceptOK->IsAccOK_sigma=0;
	acceptOK->IsAccOK_nu1=0;
	acceptOK->IsAccOK_nu2=0;
	acceptOK->IsAccOK_tau=0;
	acceptOK->IsAccOK_alpha=0;

	return acceptOK;
}

void freeIsAcceptOK(isAcceptOK *acceptOK)
{
	gsl_matrix_free(acceptOK->IsAccOK_beta);

	free(acceptOK);
}

/****************************************/

/************ NbProposals ***************/

NbProposals *createNbProposals()
{
	NbProposals *NbProp = (NbProposals *) malloc(sizeof(NbProposals));

	NbProp->NbProp_beta = gsl_matrix_calloc(2,2);

	NbProp->NbProp_betaWardOut=0;
	NbProp->NbProp_betaOutOut=0;
	NbProp->NbProp_mu=0;
	NbProp->NbProp_sigma=0;
	NbProp->NbProp_nu1=0;
	NbProp->NbProp_nu2=0;
	NbProp->NbProp_tau=0;
	NbProp->NbProp_alpha=0;

	return NbProp;
}

void reInitiateNbProp(NbProposals * NbProp)
{
	int i,j;

	for(i=0 ; i<2 ; i++)
	{
		for(j=0 ; j<2 ; j++)
		{
			gsl_matrix_set(NbProp->NbProp_beta,i,j,0);
		}
	}

	NbProp->NbProp_betaWardOut=0;
	NbProp->NbProp_betaOutOut=0;
	NbProp->NbProp_mu=0;
	NbProp->NbProp_sigma=0;
	NbProp->NbProp_nu1=0;
	NbProp->NbProp_nu2=0;
	NbProp->NbProp_tau=0;
	NbProp->NbProp_alpha=0;

}

void freeNbProposals(NbProposals *NbProp)
{
	gsl_matrix_free(NbProp->NbProp_beta);

	free(NbProp);
}

/******************************************/

/************** OUTPUT FILES **************/

output_files *createFILES(char *workspace)
{
	output_files *fich = (output_files *) malloc(sizeof(output_files));

	char fileName[300];


	strcpy(fileName, workspace);
	strcat(fileName,"LogL.txt");
	fich->LogL = fopen(fileName,"w");
	if ( fich->LogL == NULL )
	{
		printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
		exit(1);
	}

	strcpy(fileName, workspace);
	strcat(fileName,"ColonDates.txt");
	fich->ColonDates = fopen(fileName,"w");
	if ( fich->ColonDates == NULL )
	{
		printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
		exit(1);
	}

	strcpy(fileName, workspace);
	strcat(fileName,"EndColonDates.txt");
	fich->EndColonDates = fopen(fileName,"w");
	if ( fich->EndColonDates == NULL )
	{
		printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
		exit(1);
	}

	strcpy(fileName, workspace);
	strcat(fileName,"Parameters.txt");
	fich->Parameters = fopen(fileName,"w");
	if (fich->Parameters == NULL )
	{
		printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
		exit(1);
	}

	return fich;
}

void freeFILES(output_files *fich)
{
	fclose(fich->LogL);
	fclose(fich->ColonDates);
	fclose(fich->EndColonDates);
	fclose(fich->Parameters);

	free(fich);
}

/******************************************/


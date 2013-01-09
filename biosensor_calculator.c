#include <stdio.h>
#include "biosensor_calculator.h"
#include "explicit_calculator.h"
#include "implicit_calculator.h"

void calculate(struct BiosensorInformation *bio_info) 
{
  int i;	
	
  printf("explicitScheme = %i\n", bio_info->explicitScheme);
  printf("substrateInhibition = %i\n", bio_info->substrateInhibition);
  printf("productInhibition = %i\n", bio_info->productInhibition);    
  printf("k2 = %f\n", bio_info->k2);
  printf("kM = %f\n", bio_info->kM);
  printf("kS = %f\n", bio_info->kS);  
  printf("kP = %f\n", bio_info->kP);
  printf("timeStep = %f\n", bio_info->timeStep);
  printf("N = %i\n", bio_info->N);
  printf("responseTimeMethod = %i\n", bio_info->responseTimeMethod);
  printf("minTime = %f\n", bio_info->minTime);    
  printf("responseTime = %f\n", bio_info->responseTime);
  printf("outputFileName = %s\n", bio_info->outputFileName);
  printf("ne = %i\n", bio_info->ne);
  printf("s0 = %f\n", bio_info->s0);
  printf("p0 = %f\n", bio_info->p0);
  printf("noOfBiosensorLayers = %i\n", bio_info->noOfBiosensorLayers);
  
  for (i = 0; i < bio_info->noOfBiosensorLayers; i++) 
  {
	printf("biosensorLayers[%i].enzymeLayer = %i\n", i, bio_info->biosensorLayers[i].enzymeLayer);  
	printf("biosensorLayers[%i].Ds = %f\n", i, bio_info->biosensorLayers[i].Ds);  
	printf("biosensorLayers[%i].Dp = %f\n", i, bio_info->biosensorLayers[i].Dp);  
	printf("biosensorLayers[%i].d = %f\n", i, bio_info->biosensorLayers[i].d);  
	printf("biosensorLayers[%i].e0 = %f\n", i, bio_info->biosensorLayers[i].e0);  	
  }	  
  
  if (bio_info->explicitScheme)
	calculate_explicitly(bio_info);
  else
	calculate_implicitly(bio_info);  
}

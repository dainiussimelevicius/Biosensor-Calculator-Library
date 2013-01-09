#include <stdlib.h>
#include "biosensor_calculator.h"
#include "biosensor_information.h"

int main()
{
  struct BiosensorInformation *biosensorInformation;	
	
  biosensorInformation = (struct BiosensorInformation *) malloc(sizeof(struct BiosensorInformation));
  biosensorInformation->explicitScheme = 1;
  biosensorInformation->substrateInhibition = 1;
  biosensorInformation->productInhibition = 1;
  //[s^-1]
  biosensorInformation->k2 = 1;
  //[mol/l] -> [mol/cm^3]
  biosensorInformation->kM = 0.01 * 1e-3;
  //[mol/l] -> [mol/cm^3]
  biosensorInformation->kS = 0.001 * 1e-3;
  //[mol/l] -> [mol/cm^3]
  biosensorInformation->kP = 0.001 * 1e-3;
  //[s]
  biosensorInformation->timeStep = 1e-5;
  biosensorInformation->N = 200;
  biosensorInformation->responseTimeMethod = MIN_TIME;
  //[s]
  biosensorInformation->minTime = 200;
  //[s]
  biosensorInformation->responseTime = 0;
  biosensorInformation->outputFileName = "output.dat";
  biosensorInformation->ne = 1;
  //[mol/l] -> [mol/cm^3] 4e-5; 
  biosensorInformation->s0 = 0.04 * 1e-3;
  //[mol/l] -> [mol/cm^3]
  biosensorInformation->p0 = 0 * 1e-3;
  biosensorInformation->noOfBiosensorLayers = 2;
  biosensorInformation->biosensorLayers = (struct LayerInformation *) malloc(sizeof(struct LayerInformation) * 2);

  //Užpildoma sluoksnių informacija
  // 0
  biosensorInformation->biosensorLayers[0].enzymeLayer = 1;
  //[um^2/s] -> [cm^2/s]
  biosensorInformation->biosensorLayers[0].Ds = 100 * 1e-8;
  //[um^2/s] -> [cm^2/s]
  biosensorInformation->biosensorLayers[0].Dp = 100 * 1e-8;
  //[um] -> [cm]
  biosensorInformation->biosensorLayers[0].d = 10 * 1e-4;
  //[mol/l] -> [mol/cm^3]
  biosensorInformation->biosensorLayers[0].e0 = 0.01 * 1e-3;
  
  // 1  
  biosensorInformation->biosensorLayers[1].enzymeLayer = 0;
  //[um^2/s] -> [cm^2/s]
  biosensorInformation->biosensorLayers[1].Ds = 200 * 1e-8;
  //[um^2/s] -> [cm^2/s]
  biosensorInformation->biosensorLayers[1].Dp = 200 * 1e-8;
  //[um] -> [cm]
  biosensorInformation->biosensorLayers[1].d = 300 * 1e-4;
  //[mol/l] -> [mol/cm^3]
  biosensorInformation->biosensorLayers[1].e0 = 0 * 1e-3; 
  
  calculate(biosensorInformation);
    
  return 0;	
}


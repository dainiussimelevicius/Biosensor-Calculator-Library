#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "explicit_calculator.h"
#include "utils.h"
#include "constants.h"

void calculate_explicitly(struct BiosensorInformation *bio_info) 
{
  //„for“ ciklų kintamasis
  int a; 	
  
  //Srovės tankis
  double i, last_i = 0; 
  
  //Kintamasis rodo kaip pakito srovės tankis nuo praėjusios iteracijos
  double diff;
	
  //Medžiagų koncentracijų masyvai	
  double *current_s, *last_s;
  double *current_p, *last_p;
  
  //Žingsnių pagal erdvę masyvas 
  double *space_steps;
  
  //Fermentinių sluoksnių Vmax masyvas
  double *v_maxs; 
  
  //Tinklo taškų skaičius per visus biojutiklio sluoksnius
  int points;
  
  //Iteracija pagal laiką
  long long int t = 0; 
  
  //Simuliacijos laikas sekundėmis
  double execution_time;
  
  //Kintamasis nurodo ar jau pasiektas atsako laikas
  int response_time_reached;
  
  //Kintamasis nurodo ties kuriuo sluoksniu esame
  int layer;
  
  //Kintamasis nurodo, kad esame ties sluoksnio kraštu (ties sluoksnių sandūra)
  int is_boundary;
  
  //Kinetikos dedamoji
  double kinetics_part;
  
  //Rezultatų saugojimui skirtas failas
  FILE *output_file;
  
  //Sukuriamas rezultatų saugojimui skirtas failas
  output_file = fopen(bio_info->outputFileName, "w");
  fclose(output_file);  
  
  //Apskaičiuojamas tinklo taškų skaičius per visus biojutiklio sluoksnius  
  points = bio_info->noOfBiosensorLayers * bio_info->N + 1;  
  
  //Medžiagų koncentracijų masyvams išskiriama atmintis 
  last_s    = (double *) malloc(points * sizeof(double));
  current_s = (double *) malloc(points * sizeof(double));
  last_p    = (double *) malloc(points * sizeof(double));
  current_p = (double *) malloc(points * sizeof(double));
  
  //Priskiriamos pradinės ir kai kurios kraštinės sąlygos
  fill_array(last_s, points - 1, 0);
  fill_array(last_p, points - 1, 0);   
  last_s[points - 1] = bio_info->s0;
  last_p[points - 1] = bio_info->p0;  
  current_s[points - 1] = bio_info->s0;
  current_p[points - 1] = bio_info->p0;  
  current_p[0] = 0;  
    
  //Kiekvienam sluoksniui apskaičiuojami žingsniai pagal erdvę
  space_steps = (double *) malloc(bio_info->noOfBiosensorLayers * sizeof(double));
  for (a = 0; a < bio_info->noOfBiosensorLayers; a++)
    space_steps[a] = bio_info->biosensorLayers[a].d / bio_info->N;
    
  //Kiekvienam fermento sluoksniui apskaičiuojamas Vmax
  v_maxs = (double *) malloc(bio_info->noOfBiosensorLayers * sizeof(double));
  for (a = 0; a < bio_info->noOfBiosensorLayers; a++)
    if (bio_info->biosensorLayers[a].enzymeLayer)
      v_maxs[a] = bio_info->k2 * bio_info->biosensorLayers[a].e0;
      
  do
  {
    //Iteruojama per biojutiklio sluoksnius, skaičiuojamos medžiagų koncentracijos
    layer = 0;
    for (a = 1; a < points - 1; a++)
    {
	  //Nustatome ar tai nėra sluoksnių sandūra
	  is_boundary = !(a % bio_info->N);
	  
	  //Reikšmės sluoksnių sandūrose bus skaičiuojamos vėliau pagal derinimo sąlygas
	  if (is_boundary)
	  {
		//Nustatome kuriame sluoksnyje esame
        layer++;        
      }  	  
	  else
	  {
		//Įskaičiuojama difuzijos įtaka  
        current_s[a] = bio_info->biosensorLayers[layer].Ds * bio_info->timeStep * \
          (last_s[a + 1] - 2 * last_s[a] + last_s[a - 1]) / (space_steps[layer] * space_steps[layer]) + \
          last_s[a];
          
        current_p[a] = bio_info->biosensorLayers[layer].Dp * bio_info->timeStep * \
          (last_p[a + 1] - 2 * last_p[a] + last_p[a - 1]) / (space_steps[layer] * space_steps[layer]) + \
          last_p[a];          
          
        //Jeigu sluoksnis yra fermentinis, tuomet prisideda ir kinetikos dalis
        if (bio_info->biosensorLayers[layer].enzymeLayer)
        {
          kinetics_part = bio_info->timeStep * v_maxs[layer] * last_s[a] / \
            (bio_info->kM * (bio_info->productInhibition?(1 + last_p[a] / bio_info->kP):1) + \
            last_s[a] * (bio_info->substrateInhibition?(1 + last_s[a] / bio_info->kS):1));          
          current_s[a] -= kinetics_part;
          current_p[a] += kinetics_part;
        }
      }
    }    
    
    //Sluoksnių sandūroms pritaikomos derinimo sąlygos
    for (layer = 0; layer < bio_info->noOfBiosensorLayers - 1; layer++)    
    {
	  //Apskaičiuojame kuriame taške yra layer ir layer + 1 sluoksnių sandūra
      a = bio_info->N * (layer + 1);
      
      current_s[a] = (bio_info->biosensorLayers[layer + 1].Ds * space_steps[layer] * current_s[a + 1] + \
        bio_info->biosensorLayers[layer].Ds * space_steps[layer + 1] * current_s[a - 1]) / \
        (bio_info->biosensorLayers[layer + 1].Ds * space_steps[layer] + \
        bio_info->biosensorLayers[layer].Ds * space_steps[layer + 1]);
      
      current_p[a] = (bio_info->biosensorLayers[layer + 1].Dp * space_steps[layer] * current_p[a + 1] + \
        bio_info->biosensorLayers[layer].Dp * space_steps[layer + 1] * current_p[a - 1]) / \
        (bio_info->biosensorLayers[layer + 1].Dp * space_steps[layer] + \
        bio_info->biosensorLayers[layer].Dp * space_steps[layer + 1]);      
    }
	  
    //Kraštinė substrato nepratekėjimo sąlyga
    current_s[0] = current_s[1];  
	  
    //Skaičiuojamas srovės tankis
    i = bio_info->ne * F * bio_info->biosensorLayers[0].Dp * (current_p[1] - current_p[0]) / space_steps[0];
    diff = fabs(i - last_i);
    last_i = i;	     
    
    //Masyvai sukeičiami vietomis
    swap_arrays(&current_s, &last_s);
    swap_arrays(&current_p, &last_p);     
	  
	//Apskaičiuojamas laikas  
    t++;    
    execution_time = t * bio_info->timeStep;
    
    //Spausdinami rezultatai
    if ((t % INTERVAL) == 0)    
    {
      output_file = fopen(bio_info->outputFileName, "a");
      fprintf(output_file, "%e %e\n", i, execution_time);
      fclose(output_file);
    }            
        
    //Nustatoma ar tęsti simuliaciją
    switch(bio_info->responseTimeMethod) 
    {
      case MIN_TIME:
        if (execution_time < bio_info->minTime)
        {  
		  response_time_reached = 0;	
          break;		
	    }
      case DEFAULT_TIME:
        if (i > 1e-30)
          response_time_reached = ((execution_time / i) * (diff / bio_info->timeStep) <= EPSILON);
        else
          response_time_reached = 0;
        break;
      case FIXED_TIME:
        response_time_reached = (execution_time >= bio_info->responseTime);
        break;
    }    
        
  } while (!response_time_reached);
  
  //Atspausdinamas paskutinis taškas
  output_file = fopen(bio_info->outputFileName, "a");
  fprintf(output_file, "%e %e\n", i, execution_time);
  fclose(output_file);  
  
  //Atlaisvinama atmintis
  free(current_s);
  free(last_s);
  free(current_p);
  free(last_p);   	
  free(space_steps);  
  free(v_maxs);
}	

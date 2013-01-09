#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "explicit_calculator.h"
#include "utils.h"
#include "constants.h"

void calculate_explicitly(struct BiosensorInformation *bio_info) 
{
  int a; 	
    
  //Srovės tankis
  double i, last_i = 0; 
  
  //Kintamasis rodo kaip pakito srovės tankis nuo praėjusios iteracijos
  double di;
	
  //Medžiagų koncentracijų masyvai	
  double *current_s, *last_s;
  double *current_p, *last_p;
  
  //Žingsnių pagal erdvę masyvas 
  double *space_steps;
    
  //Tinklo taškų skaičius per visus biojutiklio sluoksnius
  int point_count;
  
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
  
  //Sukuriami lokalūs kintamieji dėl optimizavimo
  int subs_inh                 = bio_info->substrateInhibition;
  int prod_inh                 = bio_info->productInhibition;
  double k2                    = bio_info->k2;
  double km                    = bio_info->kM;
  double ks                    = bio_info->kS;
  double kp                    = bio_info->kP;
  double dt                    = bio_info->timeStep;
  int n                        = bio_info->N;
  enum RESP_METHOD resp_t_meth = bio_info->responseTimeMethod; 
  double min_t                 = bio_info->minTime;
  double resp_t                = bio_info->responseTime;
  char *out_file_name          = bio_info->outputFileName;
  int ne                       = bio_info->ne;
  double s0                    = bio_info->s0;
  double p0                    = bio_info->p0;
  int layer_count              = bio_info->noOfBiosensorLayers;
  int enz_layer;
  double Ds, Ds0, Ds1;
  double Dp, Dp0, Dp1;
  double dx, dx0, dx1;  
  double v_max;
    
  //Sukuriamas rezultatų saugojimui skirtas failas
  output_file = fopen(out_file_name, "w");
  fclose(output_file);  
  
  //Apskaičiuojamas tinklo taškų skaičius per visus biojutiklio sluoksnius  
  point_count = layer_count * n + 1;  
  
  //Medžiagų koncentracijų masyvams išskiriama atmintis 
  last_s    = (double *) malloc(point_count * sizeof(double));
  current_s = (double *) malloc(point_count * sizeof(double));
  last_p    = (double *) malloc(point_count * sizeof(double));
  current_p = (double *) malloc(point_count * sizeof(double));
  
  //Priskiriamos pradinės ir kai kurios kraštinės sąlygos
  fill_array(last_s, point_count - 1, 0);
  fill_array(last_p, point_count - 1, 0);   
  last_s[point_count - 1] = s0;
  last_p[point_count - 1] = p0;  
  current_s[point_count - 1] = s0;
  current_p[point_count - 1] = p0;  
  current_p[0] = 0;  
    
  //Kiekvienam sluoksniui apskaičiuojami žingsniai pagal erdvę
  space_steps = (double *) malloc(layer_count * sizeof(double));
  for (a = 0; a < layer_count; a++)
    space_steps[a] = bio_info->biosensorLayers[a].d / n;
          
  do
  {
    //Iteruojama per biojutiklio sluoksnius, skaičiuojamos medžiagų koncentracijos
    layer = 0;
    //Surenkami pirmojo sluoksnio parametrai       
    enz_layer = bio_info->biosensorLayers[layer].enzymeLayer;
    Ds = bio_info->biosensorLayers[layer].Ds;
    Dp = bio_info->biosensorLayers[layer].Dp;
    dx = space_steps[layer];            
    if (enz_layer)    
      v_max = k2 * bio_info->biosensorLayers[layer].e0;
    for (a = 1; a < point_count - 1; a++)
    {
	  //Nustatome ar tai nėra sluoksnių sandūra
	  is_boundary = !(a % n);
	  
	  //Reikšmės sluoksnių sandūrose bus skaičiuojamos vėliau pagal derinimo sąlygas
	  if (is_boundary)
	  {
		//Nustatome kuriame sluoksnyje esame
        layer++;
        //Surenkami kito sluoksnio parametrai       
        enz_layer = bio_info->biosensorLayers[layer].enzymeLayer;
        Ds = bio_info->biosensorLayers[layer].Ds;
        Dp = bio_info->biosensorLayers[layer].Dp;
        dx = space_steps[layer];            
        if (enz_layer)    
          v_max = k2 * bio_info->biosensorLayers[layer].e0;
      }  	  
	  else
	  {
		//Įskaičiuojama difuzijos įtaka  
        current_s[a] = dt * Ds * \
          (last_s[a + 1] - 2 * last_s[a] + last_s[a - 1]) / (dx * dx) + \
          last_s[a];
          
        current_p[a] = dt * Dp * \
          (last_p[a + 1] - 2 * last_p[a] + last_p[a - 1]) / (dx * dx) + \
          last_p[a];          
          
        //Jeigu sluoksnis yra fermentinis, tuomet prisideda ir kinetikos dalis
        if (enz_layer)
        {
          kinetics_part = dt * v_max * last_s[a] / \
            (km * (prod_inh?(1 + last_p[a] / kp):1) + \
            last_s[a] * (subs_inh?(1 + last_s[a] / ks):1));          
          current_s[a] -= kinetics_part;
          current_p[a] += kinetics_part;
        }
      }
    }    
    
    //Sluoksnių sandūroms pritaikomos derinimo sąlygos
    for (layer = 0; layer < layer_count - 1; layer++)    
    {
	  //Apskaičiuojame kuriame taške yra layer ir layer + 1 sluoksnių sandūra
      a = n * (layer + 1);
      
      Ds0 = bio_info->biosensorLayers[layer].Ds;
      Dp0 = bio_info->biosensorLayers[layer].Dp;
      dx0 = space_steps[layer]; 
           
      Ds1 = bio_info->biosensorLayers[layer + 1].Ds;
      Dp1 = bio_info->biosensorLayers[layer + 1].Dp;
      dx1 = space_steps[layer + 1];
      
      current_s[a] = (Ds1 * dx0 * current_s[a + 1] + Ds0 * dx1 * current_s[a - 1]) / \
                     (Ds1 * dx0 + Ds0 * dx1);      
      current_p[a] = (Dp1 * dx0 * current_p[a + 1] + Dp0 * dx1 * current_p[a - 1]) / \
                     (Dp1 * dx0 + Dp0 * dx1);      
    }
	  
    //Kraštinė substrato nepratekėjimo sąlyga
    current_s[0] = current_s[1];  
	  
    //Skaičiuojamas srovės tankis
    i = bio_info->ne * F * bio_info->biosensorLayers[0].Dp * (current_p[1] - current_p[0]) / space_steps[0];
    di = fabs(i - last_i);
    last_i = i;	     
    
    //Masyvai sukeičiami vietomis
    swap_arrays(&current_s, &last_s);
    swap_arrays(&current_p, &last_p);     
	  
	//Apskaičiuojamas laikas  
    t++;    
    execution_time = t * dt;
    
    //Spausdinami rezultatai
    if ((t % INTERVAL) == 0)    
    {
      output_file = fopen(out_file_name, "a");
      fprintf(output_file, "%e %e\n", i, execution_time);
      fclose(output_file);
    }            
        
    //Nustatoma ar tęsti simuliaciją
    switch(resp_t_meth) 
    {
      case MIN_TIME:
        if (execution_time < min_t)
        {  
		  response_time_reached = 0;	
          break;		
	    }
      case DEFAULT_TIME:
        if (i > 1e-30)
          response_time_reached = ((execution_time / i) * (di / dt) <= EPSILON);
        else
          response_time_reached = 0;
        break;
      case FIXED_TIME:
        response_time_reached = (execution_time >= resp_t);
        break;
    }    
        
  } while (!response_time_reached);
  
  //Atspausdinamas paskutinis taškas
  output_file = fopen(out_file_name, "a");
  fprintf(output_file, "%e %e\n", i, execution_time);
  fclose(output_file);  
  
  //Atlaisvinama atmintis
  free(current_s);
  free(last_s);
  free(current_p);
  free(last_p);   	
  free(space_steps);  
}	

#include <stdio.h>
#include <stdlib.h>
#include "biosensor_calculator.h"
#include "biosensor_information.h"

void callback_crunched(void *ptr, int time)
{
	printf("%ds simulated\n", time);
}

int main()
{
	struct bio_params *bio_info;

	bio_info = malloc(sizeof(*bio_info));
	bio_info->explicit_scheme = 0;
	bio_info->subs_inh = 1;
	bio_info->prod_inh = 1;
	//[s^-1]
	bio_info->k2 = 1;
	//[mol/l] -> [mol/cm^3]
	bio_info->km = 0.01 * 1e-3;
	//[mol/l] -> [mol/cm^3]
	bio_info->ks = 0.001 * 1e-3;
	//[mol/l] -> [mol/cm^3]
	bio_info->kp = 0.001 * 1e-3;
	//[s]
	bio_info->dt = 1e-3;
	bio_info->n = 200;
	bio_info->resp_t_meth = MIN_TIME;
	//[s]
	bio_info->min_t = 200;
	//[s]
	bio_info->resp_t = 0;
	bio_info->out_file_name = "output.dat";
	bio_info->ne = 1;
	//[mol/l] -> [mol/cm^3] 4e-5;
	bio_info->s0 = 0.04 * 1e-3;
	//[mol/l] -> [mol/cm^3]
	bio_info->p0 = 0 * 1e-3;
	bio_info->layer_count = 2;
	bio_info->layers = malloc(sizeof(*(bio_info->layers)) * bio_info->layer_count);

	//Užpildoma sluoksnių informacija
	// 0
	bio_info->layers[0].enz_layer = 1;
	//[um^2/s] -> [cm^2/s]
	bio_info->layers[0].Ds = 100 * 1e-8;
	//[um^2/s] -> [cm^2/s]
	bio_info->layers[0].Dp = 100 * 1e-8;
	//[um] -> [cm]
	bio_info->layers[0].d = 10 * 1e-4;
	//[mol/l] -> [mol/cm^3]
	bio_info->layers[0].e0 = 0.01 * 1e-3;

	// 1
	bio_info->layers[1].enz_layer = 0;
	//[um^2/s] -> [cm^2/s]
	bio_info->layers[1].Ds = 200 * 1e-8;
	//[um^2/s] -> [cm^2/s]
	bio_info->layers[1].Dp = 200 * 1e-8;
	//[um] -> [cm]
	bio_info->layers[1].d = 300 * 1e-4;
	//[mol/l] -> [mol/cm^3]
	bio_info->layers[1].e0 = 0 * 1e-3;

	calculate(bio_info, NULL, &callback_crunched);

	free(bio_info->layers);
	free(bio_info);

	return 0;
}

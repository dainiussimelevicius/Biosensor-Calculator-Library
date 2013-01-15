#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "implicit_calculator.h"
#include "utils.h"
#include "constants.h"

void calculate_implicitly(struct bio_params *bio_info, void *ptr, void (*callback_crunched)(void *, int))
{
	int a;

	//Srovės tankis
	double i, last_i = 0;

	//Kintamasis rodo kaip pakito srovės tankis nuo praėjusios iteracijos
	double di;

	//Medžiagų koncentracijų masyvai
	double *current_s, *last_s;
	double *current_p, *last_p;

	//Tarpiniai kintamieji naudojami sprendžiant atitinkamai substrato ir produkto lygčių sistemas
	double *pS, *qS, *pP, *qP;
	//Pagalbinis kintamasis laikysiantis reikšmę dt/(dx)^2
	double alfa;

	//Koeficientai substrato lygčių sistemoje
	double A1, A2, A3;
	//Koeficientai produkto lygčių sistemoje
	double B1, B2, B3;
	//Derinimo salygų lygties koeficientai substrato lygčių sistemoje
	double A1m, A2m, A3m;
	//Derinimo salygų lygties koeficientai produkto lygčių sistemoje
	double B1m, B2m, B3m;
	//Laisvieji nariai atitinkamai substrato ir produkto lygčių sistemose
	double fS, fP;

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
	int subs_inh                 = bio_info->subs_inh;
	int prod_inh                 = bio_info->prod_inh;
	double k2                    = bio_info->k2;
	double km                    = bio_info->km;
	double ks                    = bio_info->ks;
	double kp                    = bio_info->kp;
	double dt                    = bio_info->dt;
	int n                        = bio_info->n;
	enum resp_method resp_t_meth = bio_info->resp_t_meth;
	double min_t                 = bio_info->min_t;
	double resp_t                = bio_info->resp_t;
	char *out_file_name          = bio_info->out_file_name;
	int ne                       = bio_info->ne;
	double s0                    = bio_info->s0;
	double p0                    = bio_info->p0;
	int layer_count              = bio_info->layer_count;
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
	last_s    = malloc(point_count * sizeof(*last_s));
	current_s = malloc(point_count * sizeof(*current_s));
	last_p    = malloc(point_count * sizeof(*last_p));
	current_p = malloc(point_count * sizeof(*current_p));

	pS = malloc((point_count - 1) * sizeof(*pS));
	qS = malloc((point_count - 1) * sizeof(*qS));
	pP = malloc((point_count - 1) * sizeof(*pP));
	qP = malloc((point_count - 1) * sizeof(*qP));

	//Priskiriamos pradinės ir kraštinės sąlygos
	fill_array(last_s, point_count - 1, 0);
	fill_array(last_p, point_count - 1, 0);
	last_s[point_count - 1] = s0;
	last_p[point_count - 1] = p0;
	current_s[point_count - 1] = s0;
	current_p[point_count - 1] = p0;
	pS[0] = 1;
	qS[0] = 0;
	pP[0] = 0;
	qP[0] = 0;

	//Kiekvienam sluoksniui apskaičiuojami žingsniai pagal erdvę
	space_steps = malloc(layer_count * sizeof(*space_steps));
	for (a = 0; a < layer_count; a++)
		space_steps[a] = bio_info->layers[a].d / n;

	do {
		//Iteruojama per biojutiklio sluoksnius, skaičiuojamos medžiagų koncentracijos
		layer = 0;

		//Surenkami pirmojo sluoksnio parametrai
		enz_layer = bio_info->layers[layer].enz_layer;
		Ds = bio_info->layers[layer].Ds;
		Dp = bio_info->layers[layer].Dp;
		dx = space_steps[layer];
		if (enz_layer)
			v_max = k2 * bio_info->layers[layer].e0;
		alfa = dt / (dx * dx);
		A1 = A3 = alfa * Ds;
		A2 = -(2 * alfa * Ds + 1);
		B1 = B3 = alfa * Dp;
		B2 = -(2 * alfa * Dp + 1);
		for (a = 1; a < point_count - 1; a++) {
			//Nustatome ar tai nėra sluoksnių sandūra
			is_boundary = !(a % n);

			//Jeigu tai sluoksnių sandūra, tuomet pritaikome derinimo sąlygą
			if (!is_boundary) {
				//Jeigu tai fermento sluoksnis, tuomet įskaičiuojama kinetikos įtaka
				kinetics_part = 0;
				if (enz_layer)
					kinetics_part = dt * ((v_max * last_s[a]) / \
						(km * (prod_inh ? (1 + last_p[a] / kp) : 1) + \
						last_s[a] * (subs_inh ? (1 + last_s[a] / ks) : 1)));

				//Substrato lygčių sistema
				pS[a] = -A3 / (A1 * pS[a-1] + A2);
				fS = kinetics_part - last_s[a];
				qS[a] = (fS - A1 * qS[a - 1]) / (A1 * pS[a - 1] + A2);

				//Produkto lygčių sistema
				pP[a] = -B3 / (B1 * pP[a - 1] + B2);
				fP = -kinetics_part - last_p[a];
				qP[a] = (fP - B1 * qP[a - 1]) / (B1 * pP[a - 1] + B2);
			} else {
				//Suskaičiuojami parametrai derinimo sąlygoms
				Ds0 = bio_info->layers[layer].Ds;
				Dp0 = bio_info->layers[layer].Dp;
				dx0 = space_steps[layer];

				Ds1 = bio_info->layers[layer + 1].Ds;
				Dp1 = bio_info->layers[layer + 1].Dp;
				dx1 = space_steps[layer + 1];

				//Suskaičiuojami koeficientai derinimo sąlygoms
				A1m = Ds0 * dx1;
				A3m = Ds1 * dx0;
				A2m = -(A1m + A3m);

				B1m = Dp0 * dx1;
				B3m = Dp1 * dx0;
				B2m = -(B1m + B3m);

				//Pritaikomos derinimo sąlygos
				pS[a] = -A3m / (A1m * pS[a - 1] + A2m);
				fS = 0;
				qS[a] = (fS - A1m * qS[a - 1]) / (A1m * pS[a - 1] + A2m);

				pP[a] = -B3m / (B1m * pP[a - 1] + B2m);
				fP = 0;
				qP[a] = (fP - B1m * qP[a - 1]) / (B1m * pP[a - 1] + B2m);

				//Pereiname į kitą sluoksnį
				layer++;

				//Surenkami kito sluoksnio parametrai
				enz_layer = bio_info->layers[layer].enz_layer;
				Ds = bio_info->layers[layer].Ds;
				Dp = bio_info->layers[layer].Dp;
				dx = space_steps[layer];
				if (enz_layer)
					v_max = k2 * bio_info->layers[layer].e0;
				alfa = dt / (dx * dx);
				A1 = A3 = alfa * Ds;
				A2 = -(2 * alfa * Ds + 1);
				B1 = B3 = alfa * Dp;
				B2 = -(2 * alfa * Dp + 1);
			}
		}

		//Suskaičiuojami S ir P vektoriai
		for (a = point_count - 2; a >= 0; a--) {
			current_s[a] = pS[a] * current_s[a + 1] + qS[a];
			current_p[a] = pP[a] * current_p[a + 1] + qP[a];
		}

		//Skaičiuojamas srovės tankis
		i = bio_info->ne * F * bio_info->layers[0].Dp * \
			(current_p[1] - current_p[0]) / space_steps[0];
		di = fabs(i - last_i);
		last_i = i;

		//Masyvai sukeičiami vietomis
		swap_arrays(&current_s, &last_s);
		swap_arrays(&current_p, &last_p);

		//Apskaičiuojamas laikas
		t++;
		execution_time = t * dt;

		//Spausdinami rezultatai
		if ((t % INTERVAL) == 0) {
			output_file = fopen(out_file_name, "a");
			fprintf(output_file, "%e %e\n", i, execution_time);
			fclose(output_file);
			if (callback_crunched != NULL)
				callback_crunched(ptr, (int) execution_time);
		}

		//Nustatoma ar tęsti simuliaciją
		switch (resp_t_meth) {
		case MIN_TIME:
			if (execution_time < min_t) {
				response_time_reached = 0;
				break;
			}
			//Jeigu jau pasiekė minimalų laiką, tuomet tikrinama pagal DEFAULT_TIME sąlygas
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
	if (callback_crunched != NULL)
		callback_crunched(ptr, (int) execution_time);

	//Atlaisvinama atmintis
	free(current_s);
	free(last_s);
	free(current_p);
	free(last_p);
	free(pS);
	free(qS);
	free(pP);
	free(qP);
	free(space_steps);
}

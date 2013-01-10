#include "biosensor_calculator.h"
#include "explicit_calculator.h"
#include "implicit_calculator.h"

void calculate(struct bio_params *bio_info)
{
	if (bio_info->explicit_scheme)
		calculate_explicitly(bio_info);
	else
		calculate_implicitly(bio_info);
}

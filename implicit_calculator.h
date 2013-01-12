#ifndef IMPLICITCALCULATOR_H
#define IMPLICITCALCULATOR_H

#include "biosensor_information.h"

void calculate_implicitly(struct bio_params *bio_info, void *ptr, void (*callback_crunched)(void *, int));

#endif

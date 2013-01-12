#ifndef BIOSENSORCALCULATOR_H
#define BIOSENSORCALCULATOR_H

#include "biosensor_information.h"

void calculate(struct bio_params *bio_info, void *ptr, void (*callback_crunched)(void *, int));

#endif

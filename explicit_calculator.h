#ifndef EXPLICITCALCULATOR_H
#define EXPLICITCALCULATOR_H

#include "biosensor_information.h"

void calculate_explicitly(struct bio_params *bio_info, void *ptr, void (*callback_crunched)(void *, int));

#endif

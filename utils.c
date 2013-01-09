#include "utils.h"

void swap_arrays(double **array1, double **array2)
{
  double *temp;

  temp = *array1;
  *array1 = *array2;
  *array2 = temp;
}

void fill_array(double *array, int length, double value)
{
  int a;	 

  for (a = 0; a < length; a++)
    array[a] = value;
}

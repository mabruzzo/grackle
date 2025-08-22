/***********************************************************************
/
/ Declares the function to initialize the cloudy data
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#ifndef INITIALIZE_CLOUDY_DATA_H
#define INITIALIZE_CLOUDY_DATA_H

#include "grackle.h"

///Initializes an empty #cloudy_data struct with zeros and NULLs.
void initialize_empty_cloudy_data_struct(cloudy_data *my_cloudy);

// initialize cloudy cooling data
int initialize_cloudy_data(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           cloudy_data *my_cloudy, char *group_name,
                           code_units *my_units, int read_data);

int _free_cloudy_data(cloudy_data *my_cloudy,
                      chemistry_data *my_chemistry,
                      int primordial);

#endif /* INITIALIZE_CLOUDY_DATA_H */

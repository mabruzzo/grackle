/***********************************************************************
/
/ Declares the function to initialize the primordial rates
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#ifndef INITIALIZE_RATES_H
#define INITIALIZE_RATES_H

#include "grackle.h"

int initialize_rates(chemistry_data *my_chemistry,
                     chemistry_data_storage *my_rates,
                     code_units *my_units,
                     double co_length_unit,
                     double co_density_unit);

#endif /* INITIALIZE_RATES_H */

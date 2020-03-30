/***********************************************************************
/
/ Calculate temperature field
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* Set the mean molecular mass. */
 
#define MU_METAL 16.0
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
/* function prototypes */ 

extern void FORTRAN_NAME(calc_temp_cloudy_g)(
        gr_float *d, gr_float *e, gr_float *metal, gr_float *temperature,
	int *in, int *jn, int *kn, int *iexpand, int *imetal,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
	double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh,
        long long *priGridRank, long long *priGridDim,
        double *priPar1, double *priPar2, double *priPar3, 
        long long *priDataSize, double *priMMW, double *priForceMMW);

int local_calculate_pressure(chemistry_data *my_chemistry,
                             chemistry_data_storage *my_rates,
                             code_units *my_units,
                             grackle_field_data *my_fields,
                             gr_float *pressure);

int local_calculate_temperature_table(chemistry_data *my_chemistry,
                                      chemistry_data_storage *my_rates,
                                      code_units *my_units,
                                      grackle_field_data *my_fields,
                                      gr_float *temperature);
 
int local_calculate_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                grackle_field_data *my_fields,
                                gr_float *temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Compute the pressure first. */
 
  if (my_chemistry->primordial_chemistry > 0) {
    if (local_calculate_pressure(my_chemistry, my_rates, my_units,
                                 my_fields, temperature) == FAIL) {
      fprintf(stderr, "Error in calculate_pressure.\n");
      return FAIL;
    }
  }
 
  /* Compute the size of the fields. */
 
  int i, dim, size = 1;
  for (dim = 0; dim < my_fields->grid_rank; dim++)
    size *= my_fields->grid_dimension[dim];

  /* Calculate temperature units. */

  double temperature_units =  mh * POW(my_units->velocity_units, 2) / kboltz;

  double number_density, tiny_number = 1.-20;
  double inv_metal_mol = 1.0 / MU_METAL;
  
  if (my_chemistry->primordial_chemistry == 0) {
    if (local_calculate_temperature_table(my_chemistry, my_rates, my_units,
                                          my_fields, temperature) == FAIL) {
      fprintf(stderr, "Error in local_calculcate_temperature_table.\n");
      return FAIL;
    }
    return SUCCESS;
  }

 /* Compute temperature with mu calculated directly. */
 
# ifdef _OPENMP
# pragma omp parallel for schedule( runtime ) private( i, number_density )
# endif
  for (i = 0; i < size; i++) {
 
    if (my_chemistry->primordial_chemistry > 0) {
      number_density =
        0.25 * (my_fields->HeI_density[i] + my_fields->HeII_density[i] +
                my_fields->HeIII_density[i]) +
        my_fields->HI_density[i] + my_fields->HII_density[i] +
        my_fields->e_density[i];
    }

    /* Add in H2. */
 
    if (my_chemistry->primordial_chemistry > 1) {
      number_density += my_fields->HM_density[i] +
        0.5 * (my_fields->H2I_density[i] + my_fields->H2II_density[i]);
    }

    if (my_fields->metal_density != NULL) {
      number_density += my_fields->metal_density[i] * inv_metal_mol;
    }
 
    /* Ignore deuterium. */
 
    temperature[i] *= temperature_units / max(number_density, tiny_number);
    temperature[i] = max(temperature[i], MINIMUM_TEMPERATURE);
  }
 
  return SUCCESS;
}

int local_calculate_temperature_table(chemistry_data *my_chemistry,
                                      chemistry_data_storage *my_rates,
                                      code_units *my_units,
                                      grackle_field_data *my_fields,
                                      gr_float *temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  if (my_chemistry->primordial_chemistry != 0) {
    fprintf(stderr, "ERROR: this function requires primordial_chemistry set to 0.\n");
    return FAIL;
  }
 
  /* Compute the size of the fields. */
 
  int i, dim, size = 1;
  for (dim = 0; dim < my_fields->grid_rank; dim++)
    size *= my_fields->grid_dimension[dim];

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (my_fields->metal_density == NULL)
    metal_field_present = FALSE;

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }

  /* Calculate temperature units. */

  double temperature_units = mh * POW(my_units->velocity_units, 2) / kboltz;

  FORTRAN_NAME(calc_temp_cloudy_g)(
        my_fields->density,
        my_fields->internal_energy,
        my_fields->metal_density,
        temperature,
        my_fields->grid_dimension,
        my_fields->grid_dimension+1,
        my_fields->grid_dimension+2,
        &my_units->comoving_coordinates,
        &metal_field_present,
        my_fields->grid_start,
        my_fields->grid_start+1,
        my_fields->grid_start+2,
        my_fields->grid_end,
        my_fields->grid_end+1,
        my_fields->grid_end+2,
        &my_units->a_value,
        &my_chemistry->TemperatureStart,
        &my_chemistry->TemperatureEnd,
        &temperature_units,
        &co_length_units,
        &my_units->a_units,
        &co_density_units,
        &my_units->time_units,
        &my_chemistry->Gamma,
        &my_chemistry->HydrogenFractionByMass,
        &my_rates->cloudy_primordial.grid_rank,
        my_rates->cloudy_primordial.grid_dimension,
        my_rates->cloudy_primordial.grid_parameters[0],
        my_rates->cloudy_primordial.grid_parameters[1],
        my_rates->cloudy_primordial.grid_parameters[2],
        &my_rates->cloudy_primordial.data_size,
        my_rates->cloudy_primordial.mmw_data,
        &my_chemistry->forced_primordial_mmw);

  return SUCCESS;
}

int _calculate_temperature(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           code_units *my_units,
                           int grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                           gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                           gr_float *H2I_density, gr_float *H2II_density,
                           gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                           gr_float *e_density, gr_float *metal_density,
                           gr_float *temperature)
{

  grackle_field_data my_fields;
  my_fields.grid_rank                = grid_rank;
  my_fields.grid_dimension           = grid_dimension;
  my_fields.grid_start               = grid_start;
  my_fields.grid_end                 = grid_end;
  my_fields.density                  = density;
  my_fields.internal_energy          = internal_energy;
  my_fields.HI_density               = HI_density;
  my_fields.HII_density              = HII_density;
  my_fields.HM_density               = HM_density;
  my_fields.HeI_density              = HeI_density;
  my_fields.HeII_density             = HeII_density;
  my_fields.HeIII_density            = HeIII_density;
  my_fields.H2I_density              = H2I_density;
  my_fields.H2II_density             = H2II_density;
  my_fields.DI_density               = DI_density;
  my_fields.DII_density              = DII_density;
  my_fields.HDI_density              = HDI_density;
  my_fields.e_density                = e_density;
  my_fields.metal_density            = metal_density;

  if (local_calculate_temperature(my_chemistry, my_rates, my_units,
                                  &my_fields, temperature) == FAIL) {
    fprintf(stderr, "Error in local_calculate_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_temperature(code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *temperature)
{
  if (local_calculate_temperature(grackle_data, &grackle_rates, my_units,
                                  my_fields, temperature) == FAIL) {
    fprintf(stderr, "Error in local_calculate_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}

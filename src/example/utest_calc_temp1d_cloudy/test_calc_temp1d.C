#include <stdlib.h> // abort
#include <stdio.h>  // stderr, vfprintf
#include <stdarg.h> // va_list, va_start, va_end
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <string>
#include <vector>

#define mh     1.67262171e-24

extern "C" {
  #include <grackle.h>
  #include "../../clib/interop/calc_temp1d_cloudy_g.h"
  #include "../../clib/grackle_macros.h"

  // legacy version of the function
  extern void FORTRAN_NAME(calc_temp1d_cloudy_g)(
        const gr_float* d, const gr_float* metal, // 3D array
        const gr_float* e, // 3D array
        const double* rhoH, // 1D array
        const int* in, const int* jn, const int* kn,
        const int* is, const int* ie, const int* j, const int* k,
        double* tgas, double*mmw,
        const double* dom, const double* zr,
        const double* temstart, const double* temend, const double* gamma,
        const double* utem, const int* imetal,
        const long long* clGridRank,
        const long long* clGridDim,
        const double* clPar1,
        const double* clPar2,
        const double* clPar3,
        const long long* clDataSize,
        const double* clMMW,
        const int32_t* itmask);
}

[[noreturn]] void error(const char *fmt, ...){

  // parse argument list
  va_list arg_list;
  va_start(arg_list, fmt);         // access variable arguments after fmt
  vfprintf(stderr, fmt, arg_list); // print to stderr
  va_end(arg_list);                // cleanup variable arguments
  abort();                         // exit program with nonzero exit code
}

std::string vec_to_string(const std::vector<double>& vec) {
  std::string out = "{";

  std::size_t len = vec.size();

  std::size_t pause_start;
  std::size_t pause_stop;

  if (len > 30){
    pause_start = 3;
    pause_stop = len - 3;
  } else {
    pause_start = len *2;
    pause_stop = pause_start;
  }

  for (std::size_t i = 0; i < len; i++) {
    if ((i > pause_start) && (i < pause_stop)) { continue; }

    if (i == pause_stop) {
      out += ", ... ";
    } else if (i != 0) {
      out += ", ";
    }

    char buf[30];
    sprintf(buf, "%g", vec[i]);
    out += buf;
  }
  return out + "}";
}

void compare_values(const std::vector<double>& actual,
                    const std::vector<double>& desired,
                    double rtol = 0.0, double atol = 0.0,
                    std::string err_msg = "")
{
  if (actual.size() != desired.size()){
    error("the compared arrays have different lengths\n");
  }

  std::size_t num_mismatches = 0;
  double max_absDiff = 0.0;
  std::size_t max_absDiff_ind = 0;
  double max_relDiff = 0.0;
  std::size_t max_relDiff_ind = 0;

  for (std::size_t i = 0; i < actual.size(); i++) {
    double cur_absDiff = fabs(actual[i]-desired[i]);

    if ( cur_absDiff > (atol + rtol * fabs(desired[i])) ) {
      num_mismatches++;
      if (cur_absDiff > max_absDiff){
        max_absDiff = cur_absDiff;
        max_absDiff_ind = i;
      }

      if ( cur_absDiff > (max_relDiff * fabs(desired[i])) ) { // no div-by-zero
        max_relDiff = cur_absDiff / fabs(desired[i]);
        max_relDiff_ind = i;
      }
    }
  }

  if (num_mismatches == 0) { return; }

  std::string actual_vec_str = vec_to_string(actual);
  std::string ref_vec_str = vec_to_string(desired);

  error
    (("arrays are unequal for the tolerance: rtol = %g, atol = %g\n"
      "%s\n" // custom error message
      "Mismatched elements: %d / %d\n"
      "Max absolute difference: %g,   ind = %d, actual = %g, reference = %g\n"
      "Max relative difference: %g,   ind = %d, actual = %g, reference = %g\n"
      "actual: %s\n"
      "desired: %s\n"),
     rtol, atol, (int)num_mismatches, (int)actual.size(),
     err_msg.c_str(),
     max_absDiff, (int)max_absDiff_ind, actual[max_absDiff_ind],
     desired[max_absDiff_ind],
     max_relDiff, (int)max_relDiff_ind, actual[max_relDiff_ind],
     desired[max_relDiff_ind],
     actual_vec_str.c_str(), ref_vec_str.c_str());
}

void cut_down_to_1D_table(cloudy_data* ptr){
  if (ptr->grid_rank != 2) { error("can currently only cut down 2D to 1D"); }

  // table over rho and Temperature
  // -> I confirmed from the mmw table that temperature axis is contiguous
  long long num_rho_vals = ptr->grid_dimension[0];
  long long num_T_vals = ptr->grid_dimension[1];

  cloudy_data newObj;
  newObj.grid_rank = 1;

  // Dimension of dataset
  newObj.grid_dimension
    = (long long *)malloc(sizeof(long long) * newObj.grid_rank);
  newObj.grid_dimension[0] = num_T_vals;
  GRACKLE_FREE(ptr->grid_dimension);

  // Dataset parameter values.
  newObj.grid_parameters
    = (double**)malloc(sizeof(double**) * 1);
  newObj.grid_parameters[0] = ptr->grid_parameters[1];
  GRACKLE_FREE(ptr->grid_parameters[0]);
  GRACKLE_FREE(ptr->grid_parameters);

  // since Temperature is the fast-index, we are going to just reuse the
  // original pointers (even though they will have a bunch of unused data)
  newObj.heating_data = ptr->heating_data;
  newObj.cooling_data = ptr->cooling_data;
  newObj.mmw_data = ptr->mmw_data;
  newObj.data_size = num_T_vals;

  (*ptr) = newObj;
}

struct DummyGrackleConfig{
  // the central purpose here is to hold grackle configuration options for
  // tabulated solver
  //
  // do NOT try to use this as a general purpose C++ interface

  chemistry_data chem_data;
  chemistry_data_storage chem_rates;
  code_units units;

  DummyGrackleConfig(int n_tab_dims, double radiation_redshit) {
    // radiation_redshift is meaningless when n_tab_dims isn't 3
    
    // setup units!
    code_units my_units;
    my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
    my_units.density_units = mh;                  // = mass_H/cm^3 (in g/cm^3)
    my_units.length_units = 1.0 / 3.24077929e-25; // = 1 Mpc (in cm)
    my_units.time_units = 3.15576e13;             // = 1 Myr (in seconds)
    my_units.a_units = 1.0; // units for the expansion factor
    // Set expansion factor to 1 for non-cosmological simulation.
    my_units.a_value = 1. / (1. + radiation_redshit) / my_units.a_units;
    set_velocity_units(&my_units);

    // setup runtime parameters
    chemistry_data my_chem = _set_default_chemistry_parameters();
    my_chem.use_grackle = 1;            // chemistry on
    my_chem.with_radiative_cooling = 1; // cooling on
    my_chem.primordial_chemistry = 0;
    my_chem.metal_cooling = 1;         // metal cooling on
    my_chem.UVbackground = 0;          // UV background on

    // can't directly assign a string literal to char*
    if ((n_tab_dims <= 0) || (n_tab_dims > 3)) {
      error("n_tab_dims must lie be 1, 2, or 3\n");
    } else if (n_tab_dims <= 2) {
      const char* tmp = "../../grackle_data_files/input/CloudyData_noUVB.h5";
      my_chem.grackle_data_file = new char[strlen(tmp)+1];
      strcpy(my_chem.grackle_data_file, tmp);
    } else {
      const char* tmp = "../../grackle_data_files/input/CloudyData_UVB=HM2012.h5";
      my_chem.grackle_data_file = new char[strlen(tmp)+1];
      strcpy(my_chem.grackle_data_file, tmp);
    }

    chemistry_data_storage my_rates;
    _initialize_chemistry_data(&my_chem, &my_rates, &my_units);

    if (n_tab_dims == 1) {
      cut_down_to_1D_table(&(my_rates.cloudy_metal));
      cut_down_to_1D_table(&(my_rates.cloudy_primordial));
      //error("Can't currently support n_tab_dims == 1\n");
    }

    this->chem_data = my_chem;
    this->chem_rates = my_rates;
    this->units = my_units;
  }

  ~DummyGrackleConfig() {
    _free_chemistry_data(&(this->chem_data), &(this->chem_rates));
    delete[] (this->chem_data).grackle_data_file;
  }

  DummyGrackleConfig(const DummyGrackleConfig&) = delete;
  DummyGrackleConfig(DummyGrackleConfig&&) = delete;
  DummyGrackleConfig& operator=(const DummyGrackleConfig&) = delete;
  DummyGrackleConfig& operator=(DummyGrackleConfig&&) = delete;

};


struct calc_temp_outputs{
  std::vector<double> tgas;
  std::vector<double> mmw;
};

calc_temp_outputs run_test(DummyGrackleConfig& config,
                           bool use_fortran, bool slc_from_3D_arr = true){

  /* Calculate temperature units. */
  const double temperature_units = get_temperature_units(&(config.units));

  // initialize physical quantities!
  const gr_float Temp0 = 1.0e3;   // background temperature (in K)

  const std::size_t length = 60;
  std::vector<gr_float> density(length);
  std::vector<gr_float> metal_density(length);
  std::vector<gr_float> eint(length);
  std::vector<double> rhoH(length); // double is NOT a typo
  for (std::size_t i = 0; i < length; i++) {
    density[i] = gr_float(i+1);
    metal_density[i] = gr_float(config.chem_data.SolarMetalFractionByMass *
                                density[i]);
    rhoH[i] = config.chem_data.HydrogenFractionByMass * density[i];
    eint[i] = gr_float(Temp0 / temperature_units);
  }

  // prepare some other arguments
  //    -> when running in 3D, we aren't
  const int in = (!slc_from_3D_arr) ? int(length) : 5;
  const int jn = (!slc_from_3D_arr) ? 1           : 4;
  const int kn = (!slc_from_3D_arr) ? 1           : 3;

  if (std::size_t(in * jn * kn) > length) { error("something is wrong\n"); }

  const int j = jn; // remember, these are currently 1-indexed
  const int k = kn;

  const int is = 0; // not a typo! (even though other indices are 1-indexed)
  const int ie = in - 1; // not a typo!

  const double aye = config.units.a_value; // expansion factor (in code units)
  const double dom = config.units.density_units * pow(aye,3) / mh;
  const double zr = 1.0/(aye * config.units.a_units) - 1.;
  const int imetal = 1;

  // prepare the output arguments
  std::vector<double> tgas(in);
  std::vector<double> mmw(in);

  const std::vector<int32_t> itmask(std::size_t(in), int32_t(1));

  if (use_fortran) {
    FORTRAN_NAME(calc_temp1d_cloudy_g)
      (density.data(), metal_density.data(), eint.data(),
       rhoH.data(),
       &in, &jn, &kn,  &is, &ie, &j, &k,
       tgas.data(), mmw.data(),
       &dom, &zr,
       &config.chem_data.TemperatureStart, &config.chem_data.TemperatureEnd,
       &config.chem_data.Gamma, &temperature_units, &imetal,
       &config.chem_rates.cloudy_primordial.grid_rank, // clGridRank
       config.chem_rates.cloudy_primordial.grid_dimension, // clGridDim
       config.chem_rates.cloudy_primordial.grid_parameters[0], // clPar1
       config.chem_rates.cloudy_primordial.grid_parameters[1], // clPar2
       config.chem_rates.cloudy_primordial.grid_parameters[2], // clPar3
       &config.chem_rates.cloudy_primordial.data_size, // clDataSize
       config.chem_rates.cloudy_primordial.mmw_data, // clMMW
       itmask.data());

  } else {
    calc_temp1d_cloudy_g
      (density.data(), metal_density.data(), eint.data(),
       rhoH.data(),
       in, jn, kn,  is, ie, j, k,
       tgas.data(), mmw.data(),
       dom, zr,
       config.chem_data.TemperatureStart, config.chem_data.TemperatureEnd,
       config.chem_data.Gamma, temperature_units, imetal,
       config.chem_rates.cloudy_primordial.grid_rank, // clGridRank
       config.chem_rates.cloudy_primordial.grid_dimension, // clGridDim
       config.chem_rates.cloudy_primordial.grid_parameters[0], // clPar1
       config.chem_rates.cloudy_primordial.grid_parameters[1], // clPar2
       config.chem_rates.cloudy_primordial.grid_parameters[2], // clPar3
       config.chem_rates.cloudy_primordial.data_size, // clDataSize
       config.chem_rates.cloudy_primordial.mmw_data, // clMMW
       itmask.data());
  }

  return {tgas, mmw};
}

int main(void){
  for (int n_tab_dims = 1; n_tab_dims < 4; n_tab_dims++) {
    printf("\nConsidering a %dD table of MMW values\n", n_tab_dims);

    std::vector<double> z_vals {0.0};
    if (n_tab_dims == 3) {
      z_vals.push_back(0.13242);
      z_vals.push_back(0.5);
      z_vals.push_back(10);

      DummyGrackleConfig conf(n_tab_dims,0.0);
      double* clPar2 = conf.chem_rates.cloudy_primordial.grid_parameters[1];
      long long* clGridDim = conf.chem_rates.cloudy_primordial.grid_dimension;
      z_vals.push_back(clPar2[clGridDim[1]-3]);
      z_vals.push_back(clPar2[clGridDim[1]-2]);
      z_vals.push_back(clPar2[clGridDim[1]-1]);
      z_vals.push_back(clPar2[clGridDim[1]-1] * 1.01);
    }

    for (const auto& z_val : z_vals) {
      DummyGrackleConfig config(n_tab_dims,z_val);

      for (bool slc_from_3D_arr : {false, true}) {
        const char* descr = (slc_from_3D_arr) ? "slice from 3D arr" : "1D arr";
        printf("-> z = %g, comparing %s\n", z_val, descr);

        //printf("using the c version:\n"); fflush(stdout);
        calc_temp_outputs actual = run_test(config, false, slc_from_3D_arr);

        //printf("using the fortran version:\n"); fflush(stdout);
        calc_temp_outputs reference = run_test(config, true, slc_from_3D_arr);

        compare_values(actual.mmw, reference.mmw, 0.0, 0.0,
                       "**Error during comparison of mmw**");
        compare_values(actual.tgas, reference.tgas, 0.0, 0.0,
                       "**Error during comparison of tgas**");
      }
    }
  }

  return 0;
}

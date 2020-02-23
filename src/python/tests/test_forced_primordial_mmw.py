########################################################################
#
# Tests behavior when forced_primordial_mmw is used
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np
import os

from pygrackle import \
    chemistry_data, \
    FluidContainer

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    boltzmann_constant_cgs

from pygrackle.utilities.testing import \
    random_logscale, \
    assert_rel_equal

def test_forced_primordial_mmw():
    # Set solver parameters
    grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    data_file_path = bytearray(os.sep.join(
        [grackle_dir, "input", "CloudyData_UVB=HM2012.h5"]), 'utf-8')

    my_rand_state = np.random.RandomState(7921)

    for i in range(3):
        chem = chemistry_data()
        chem.use_grackle = 1
        chem.primordial_chemistry = 0
        chem.UVbackground = 1
        chem.self_shielding_method = 0
        chem.H2_self_shielding = 0
        chem.grackle_data_file = data_file_path
        chem.use_specific_heating_rate = 0
        chem.use_volumetric_heating_rate = 0

        # initialize mmw in [0.6,1.0)
        chem.forced_primordial_mmw = 0.4 * my_rand_state.random_sample(1) + 0.6

        chem.comoving_coordinates = 0 # proper units
        chem.a_units = 1.0
        chem.a_value = 1.0
        chem.density_units = random_logscale(-26,-22,
                                             random_state = my_rand_state)
        chem.length_units = random_logscale(12, 18,
                                            random_state = my_rand_state)
        chem.time_units = random_logscale(7, 15,
                                          random_state = my_rand_state)
        chem.velocity_units = chem.length_units / chem.time_units

        # set up a fluid container.
        rval = chem.initialize()
        if rval == 0:
            raise RuntimeError("Failed to initialize chemistry_data.")

        n_points = 200
        temperature = random_logscale(1, 9, size = n_points,
                                      random_state = my_rand_state)

        fc = FluidContainer(chem, n_points)
        fc["density"][:] = 10*my_rand_state.random_sample(1) # between 0 and 10
        fc["energy"] = temperature / chem.temperature_units / \
                       chem.forced_primordial_mmw / (chem.Gamma - 1.0)
        fc["x-velocity"][:] = 0.0
        fc["y-velocity"][:] = 0.0
        fc["z-velocity"][:] = 0.0
        fc.calculate_temperature()
        mu = fc["temperature"] / \
             (fc["energy"] * chem.temperature_units *
              (chem.Gamma - 1) )

        # correct for the difference in values between the boltzmann constant
        # and hydrogen mass used in the C/Fortran routines and the values used
        # in the python layer (used to compute chem.temperature_units)
        corr = ( (mass_hydrogen_cgs / boltzmann_constant_cgs) /
                 (1.67262171e-24    / 1.3806504e-16         ) )

        assert_rel_equal(
            mu*corr, (np.zeros_like(mu) + chem.forced_primordial_mmw), 15,
            ('The mean molecular weight computed from energy and temperature '
             'disagrees with chem.forced_primordial_mmw'))

import pandas as pd
import matplotlib as plt
import numpy as np
from scipy import integrate

evtoerg = 0.16 * 10 ** -11  # SED points unit is (erg cm-2 s-1)

# FermiPy returns *.npy file where gathered all information about SED of the source
file = np.load("4fgl_j1454.0+4927_sed.npy", encoding="latin1", allow_pickle=True).flat[0]

"""prefactor = file["param_values"][0]
prefactor_err = file["param_errors"][0]

index = file["param_values"][1]
index_err = file["param_errors"][1]

scale = file["param_values"][2]
"""
# In *.npy file e_min has "e_min" key name and the same for e_max
e_min = file["e_min"]
e_ref = file["e_ref"]
e_max = file["e_max"]
print(e_max)
""""
From .fits file flux and flux_err should be got.
Flux unit in that .fits is Mev cm-2 s-1.
"""
e2dnde = file["e2dnde"]
e2dnde_err = file["e2dnde_err"]

sed_point_erg = e2dnde * evtoerg * 10 ** 6  # Convert Mev flux to Erg flux
sed_point_err_erg = e2dnde_err * evtoerg * 10 ** 6  # Convert Mev flux_err to Erg flux_err


"""
resmin.append(integrate.quad(lambda x, prefactor, scale, index:
                             (x**2) * powerlaw(x, prefact, scale, index),
                             range_list[i], range_list[i + 1],
                             args=((prefactor - prefactor_err),
                                   (-1) * index)))
                                   
resmin.append(integrate.quad(lambda x, prefact, scale, index:
                             (x**2) * powerlaw(x, prefact, scale, index),
                             range_list[i], range_list[i + 1],
                             args=((prefactor[i] - deltaprefactor[i]) * prefactor_scale[i], power_scale[i],
                                   (-1) * power_index[i])))
resmax.append(integrate.quad(lambda x, prefact, scale, index:
                             (x**2) * powerlaw(x, prefact, scale, index),
                             range_list[i], range_list[i + 1],
                             args=((prefactor[i] + deltaprefactor[i]) * prefactor_scale[i], power_scale[i],
                                   (-1) * power_index[i])))

resmid.append(integrate.quad(lambda x, prefact, scale, index:
                             (x**2) * powerlaw(x, prefact, scale, index),
                             range_list[i], range_list[i + 1],
                             args=(
                                 prefactor[i] * prefactor_scale[i], power_scale[i], (-1) * power_index[i])))
"""

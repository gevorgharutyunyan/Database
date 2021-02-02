import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

evtoerg = 0.16 * 10 ** -11  # SED points unit is (erg cm-2 s-1)
def plot_sed():
    # FermiPy returns *.npy file where gathered all information about SED of the source
    file = np.load("4fgl_j1454.0+4927_sed.npy", encoding="latin1", allow_pickle=True).flat[0]

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

    sed_points = {'e_min':e_min,
                'e_ref':e_ref,
                'e_max':e_max,
                'sed_point_erg':sed_point_erg,
                'sed_point_err_erg':sed_point_err_erg}

    sed_data= pd.DataFrame(sed_points)



    fig, ax = plt.subplots(nrows =1 ,ncols = 1,figsize=(18, 9))
    ax.errorbar(sed_data.e_ref, sed_data.sed_point_erg, xerr=[e_ref-e_min,e_max-e_ref ], yerr=sed_point_err_erg, fmt=".-", elinewidth=1.002, ls=":",
                lw=0,
                color="r",
                capsize=1.98)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

plot_sed()
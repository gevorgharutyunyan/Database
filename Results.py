"""
Welcome to Results.py
The Python version of SEDResults.py is 2.7
Catalog version is 4FGL
Put this code into your working directory to draw SED results and save them.


--Version 1.1.0--
***Authors
---Gevorg Harutyunyan
---David  Israyelyan
---Mher   Khachatryan
***
"""
import json
import os
import re
import smtplib  # for e-mailing

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from scipy import integrate

evtoerg = 0.16 * 10 ** -11  # to convert ev to erg for changing axis values
upper_limit_ts = 2  # upper limit bound


#####################################################################################
#                                                                                   #
#                                ***Functions***                                    #
#                                                                                   #
#####################################################################################

def range_counter(emin, emax, range_num):
    """Divide energies to ranges """
    global range_list
    range_list = []  # list to store endpoints of energy ranges
    for i in range(0, range_num + 1):
        range_endpoint = 10 ** (np.log10(emin) + i * (np.log10(emax) - np.log10(emin)) / range_num)
        range_endpoint = np.round(range_endpoint, decimals=2)
        range_list.append(range_endpoint)
    return range_list


def bin_counter(emin, emax, range_count):
    """Compute enumbins"""
    low_limit = np.log10(100)
    upper_limit = np.log10(10 ** 6)
    cube_num = 40
    step = (upper_limit - low_limit) / cube_num
    bins_list = []  # list of every bin from start to end,
    count = 0  # for counting bins number
    # compute endpoints of bins
    for i in range(cube_num + 1):
        bins_list.append(np.round(10 ** (low_limit + i * step), decimals=2))
    for j in range(len(range_counter(emin, emax, range_count)) - 1):
        for bins in bins_list:
            if range_counter(emin, emax, range_count)[j] <= bins < \
                    range_counter(emin, emax, range_count)[j + 1]:
                count += 1
    return count


def powerlaw(e, power_prefactor, scale, index):
    """formula of power law"""
    return power_prefactor * (e / scale) ** index


def lc_bin_counter(tmin, tmax, range_num):
    """Compute binning of light curve"""
    global tmin_list
    tmin_list = np.array([])
    val = tmin
    tstep = (float(tmax - tmin)) / float(range_num)

    if isinstance(tstep, float):
        tmax = np.floor(tstep) * range_num + tmin

    if tmin != tmax:
        while val in np.arange(tmin, tmax):
            tmin_list = np.append(tmin_list, val)
            val += int(tstep)

        tmin_list = np.append(tmin_list, tmax)

    else:
        tmin_list = np.append(tmin_list, tmin)


######################################################################################
#                                                                                    #
#                                  ***READ RESULTS***                                #
#                                                                                    #
######################################################################################

def power_law_upper(e, power_scale, power_index):
    """Power law function for upper limits"""
    return (e / power_scale) ** power_index


def upper_limit_calc(power_scale_var, power_index_var, e_min, e_max, upper_limit_flux):
    """Change upper limit to [ERG*CM**-2*S**-1]"""
    diff_flux = integrate.quad(power_law_upper, e_min, e_max, args=(power_scale_var, power_index_var))[0]
    pref = upper_limit_flux / diff_flux
    res = \
        integrate.quad(lambda x, prefactor: x * power_law_upper(x, power_scale_var, power_index_var) * prefactor,
                       e_min,
                       e_max, args=pref)
    return res


def read_sed_results():
    """Read Results from result files of sed"""
    global resmin, resmax, err1, err2, fluxes_axis, fluxes_error_axis, prefactor, power_index, power_deltaindex
    global power_scale, ts_list, fluxes, fluxeserr, deltaprefactor, prefactor_scale
    # creating some lists  for appending corresponding values
    prefactor = []  # for taking prefactors
    deltaprefactor = []  # for taking prefactors
    prefactor_scale = []  # for taking prefactors
    MODELS = []  # for taking all models
    power_index = []  # to get power indexes
    power_deltaindex = []  # to get delta index
    power_scale = []  # to get power scales
    fluxes = []  # to get fluxes in [PHOTON*CM^-2*S^-1]
    fluxeserr = []  # to get errors of fluxes in [PHOTON*CM^-2*S^-1]
    resmin = []  # creating empty lists to store integrated values to get flux and flux error
    resmax = []
    resmid = []
    err1 = []
    err2 = []
    fluxes_axis = []  # flux values in [erg*sm^-2*s^-1]
    fluxes_error_axis = []  # errors of flux in [erg*sm^-2*s^-1]
    ts_list = []  # store data of TS'
    for i in range(range_count):
        '''read fluxes from files'''
        with open("range{i}_sed_fit_result.txt".format(i=i + 1), 'r') as range0:
            lines = range0.readlines()
        lines = [x.strip() for x in lines]  # get rid of whitespaces or /n's

        for line in lines:

            match_object = re.match("Spectrum.*", line)  # get the line of spectrum
            if match_object:  # if in Model's first line
                spectrum_index = lines.index(match_object.group())  # get index of Model's first line
                # regex is used to prevent only using power law, it's for seeing Spectrum and what comes aftrer it

                for i in range(spectrum_index, len(lines)):  # to get Model lines including Spectrum type to FLUX line
                    MODELS.append(lines[i])

            if line == "FLUX [PHOTON*CM^-2*S^-1]":
                flux_index = lines.index("FLUX [PHOTON*CM^-2*S^-1]")
                flux = lines[flux_index + 1]
                fluxes.append(float(flux))
            elif line == "FLUX ERROR [PHOTON*CM^-2*S^-1]":
                fluxerr_index = lines.index("FLUX ERROR [PHOTON*CM^-2*S^-1]")
                fluxerr = lines[fluxerr_index + 1]
                fluxeserr.append(float(fluxerr))

            elif line == "TS":
                ts_index = lines.index("TS")
                ts = lines[ts_index + 1]
                ts_list.append(float(ts))
        # iterate in MODELS list which continues all models as list of string
    for model in MODELS:
        '''get other parameters from model'''
        match_prefactor = re.match(".*Prefactor.*", model)  # get all lines of Prefactors
        if match_prefactor:
            models_prefactor = model.split(" ")  # split string to get lists
            models_prefactor = list(filter(None, models_prefactor))  # get rid of empty strings in list
            prefactor.append(models_prefactor[2])  # get prefactors
            deltaprefactor.append(models_prefactor[3])  # get prefactors
            prefactor_scale.append(models_prefactor[7][:-1])  # get prefactors
        match_index = re.match(".*Index.*", model)  # get all lines of Index
        if match_index:
            models_index = model.split(" ")  # split string to get lists
            models_index = list(filter(None, models_index))  # get rid of empty strings in list
            power_index.append(models_index[2])  # get indexes
            power_deltaindex.append(models_index[3])  # get errors of indixes
        match_scale = match_index = re.match(".*Scale.*", model)  # get all lines of Scale
        if match_scale:
            models_scale = model.split(" ")  # split string to get lists
            models_scale = list(filter(None, models_scale))  # get rid of empty strings in list
            power_scale.append(models_scale[2])  # get scales

    # convert parameters taken from model of .txt files to numpy array
    prefactor = np.array(prefactor, dtype=float)
    deltaprefactor = np.array(deltaprefactor, dtype=float)
    prefactor_scale = np.array(prefactor_scale, dtype=float)
    power_index = np.array(power_index, dtype=float)
    power_deltaindex = np.array(power_deltaindex, dtype=float)
    power_scale = np.array(power_scale, dtype=float)
    ts_list = np.array(ts_list, dtype=float)
    fluxes = np.array(fluxes)
    fluxeserr = np.array(fluxeserr)

    # integrate fluxes to get in [ERG*SM^-2*S^-1]
    for i in range(len(prefactor)):
        if ts_list[i] < upper_limit_ts:
            var1 = upper_limit_calc(power_scale[i], power_index[i], range_list[i], range_list[i + 1], fluxes[i])
            resmid.append(var1)
            resmax.append(var1)
            resmin.append(var1)
        else:
            resmin.append(integrate.quad(lambda x, prefact, scale, index:
                                         x * powerlaw(x, prefact, scale, index),
                                         range_list[i], range_list[i + 1],
                                         args=((prefactor[i] - deltaprefactor[i]) * prefactor_scale[i], power_scale[i],
                                               (-1) * power_index[i])))
            resmax.append(integrate.quad(lambda x, prefact, scale, index:
                                         x * powerlaw(x, prefact, scale, index),
                                         range_list[i], range_list[i + 1],
                                         args=((prefactor[i] + deltaprefactor[i]) * prefactor_scale[i], power_scale[i],
                                               (-1) * power_index[i])))

            resmid.append(integrate.quad(lambda x, prefact, scale, index:
                                         x * powerlaw(x, prefact, scale, index),
                                         range_list[i], range_list[i + 1],
                                         args=(
                                             prefactor[i] * prefactor_scale[i], power_scale[i], (-1) * power_index[i])))

        err1.append(resmid[i][0] * evtoerg * 10 ** 6 - resmin[i][0] * evtoerg * 10 ** 6)
        err2.append(resmax[i][0] * evtoerg * 10 ** 6 - resmid[i][0] * evtoerg * 10 ** 6)
        fluxes_axis.append(resmid[i][0] * evtoerg * 10 ** 6)
        fluxes_error_axis.append((err1[i] + err2[i]) / 2)

    resmin = np.array(resmin)
    resmax = np.array(resmax)
    resmid = np.array(resmid)
    err1 = np.array(err1)
    err2 = np.array(err2)
    fluxes_axis = np.array(fluxes_axis)
    fluxes_error_axis = np.array(fluxes_error_axis)


def read_lc_results():
    """Read Results from result files of lc"""
    global fluxes_lc, fluxeserr_lc, power_index_lc, power_deltaindex_lc, ts_list_lc
    fluxes_lc = []
    fluxeserr_lc = []
    MODELS = []  # for taking all models
    power_index_lc = []  # to get power indexes
    power_deltaindex_lc = []  # to get delta index
    ts_list_lc = []  # store data of TS'

    for i in range(len(tmin_list) - 1):
        '''read fluxes from files'''
        with open("range{i}_lc_fit_result.txt".format(i=i + 1), 'r') as range_lc:
            lines = range_lc.readlines()
        lines = [x.strip() for x in lines]  # get rid of whitespaces or /n's

        for line in lines:

            match_object = re.match("Spectrum.*", line)  # get the line of spectrum
            if match_object:  # if in Model's first line
                spectrum_index = lines.index(match_object.group())  # get index of Model's first line
                # regex is used to prevent only using power law, it's for seeing Spectrum and what comes aftrer it

                for i in range(spectrum_index, len(lines)):  # to get Model lines including Spectrum type to FLUX line
                    MODELS.append(lines[i])

            if line == "FLUX [PHOTON*CM^-2*S^-1]":
                flux_index = lines.index("FLUX [PHOTON*CM^-2*S^-1]")
                flux_lc = lines[flux_index + 1]
                fluxes_lc.append(float(flux_lc))
            elif line == "FLUX ERROR [PHOTON*CM^-2*S^-1]":
                fluxerr_index = lines.index("FLUX ERROR [PHOTON*CM^-2*S^-1]")
                fluxerr_lc = lines[fluxerr_index + 1]
                fluxeserr_lc.append(float(fluxerr_lc))

            elif line == "TS":
                ts_index = lines.index("TS")
                ts = lines[ts_index + 1]
                ts_list_lc.append(float(ts))

    for model in MODELS:
        '''get other parameters from model'''

        match_index = re.match(".*Index.*", model)  # get all lines of Index
        if match_index:
            models_index = model.split(" ")  # split string to get lists
            models_index = list(filter(None, models_index))  # get rid of empty strings in list
            power_index_lc.append(models_index[2])  # get indexes
            power_deltaindex_lc.append(models_index[3])  # get errors of indixes
        match_scale = match_index = re.match(".*Scale.*", model)  # get all lines of Scale

    # convert parameters taken from model of .txt files to numpy array

    power_index_lc = np.array(power_index_lc, dtype=float)
    power_deltaindex_lc = np.array(power_deltaindex_lc, dtype=float)
    ts_list_lc = np.array(ts_list_lc, dtype=float)
    fluxes_lc = np.array(fluxes_lc, dtype=float)
    fluxeserr_lc = np.array(fluxeserr_lc, dtype=float)


######################################################################################
#                                                                                    #
#                                  ***Plotting***                                    #
#                                                                                    #
######################################################################################


sed_or_lc = raw_input("Type SED for plotting SED.\nType LC for plotting Light Curve: ")
while True:
    if sed_or_lc.lower() in ['sed', 's']:
        # read object properties created previously
        with open("object_properties_sed.txt", "r") as tx_r:
            properties = json.load(tx_r)

        # get object name and some properties to use them in code
        obj_name = properties.get("obj_name")
        tmin = properties.get("tmin")
        tmax = properties.get("tmax")
        emin = properties.get("emin")
        emax = properties.get("emax")
        range_count = properties.get("range_count")
        receiver_email = properties.get("user_email")

        # get number of bins
        enumbins = bin_counter(emin, emax, range_count)  # counts bins
        range_counter(emin, emax, range_count)  # returns list of ranges
        read_sed_results()

        # get energy axis

        energy = []  # to get energies in given ranges
        xminerr = []  # minimum values of energy
        xmaxerr = []  # maximum values of energy
        energyerr = []  # to store energy errors

        # get x axis (energy axis) values, xmin is lower limit of range
        # xmax is upper limit of range
        # xcenter is center part of energy between xmin and xmax
        for i in range(len(range_list) - 1):
            xmin = range_list[i]
            xmax = range_list[i + 1]
            xcenter = np.sqrt(xmin * xmax)
            energy.append(xcenter)
            xminerr.append(xcenter - xmin)
            xmaxerr.append(xmax - xcenter)

        energy = np.array(energy)  # array of energy
        energyerr.append(xminerr)  # append lower bound of energy error to list
        energyerr.append(xmaxerr)  # append upper bound of energy error to list

        # plot sed

        fig, ax = plt.subplots(figsize=(18, 9))
        fontdict = {'fontsize': 18, "fontweight": 'bold'}  # customize axis names
        ax.errorbar(energy, fluxes_axis, xerr=energyerr, yerr=fluxes_error_axis, fmt=".-", elinewidth=1.002, ls=":",
                    lw=0,
                    color="r",
                    capsize=1.98)  # plot an errorbar of SED
        ax.set_xscale("log")
        ax.set_yscale("log")
        title_name = obj_name.replace("_", " ")  # change  object name to readable one
        title_name = title_name.upper()
        ax.set_title("{title_name} SED".format(title_name=title_name), fontdict=fontdict, loc='center')
        ax.set_ylabel("FLUX [ERG*SM^-2*S^-1] ", fontdict=fontdict)
        ax.set_xlabel("Energy, [MeV]", fontdict=fontdict)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(16)

        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(16)

        plt.savefig("{obj_name}_SED.pdf".format(obj_name=obj_name))

        # save results to dataframe

        sed_df = pd.DataFrame()
        energymin = np.array([range_list[i] for i in range(len(range_list) - 1)])
        energymax = np.array([range_list[i + 1] for i in range(len(range_list) - 1)])
        energymid = np.sqrt(energymin * energymax)
        sed_df["Energy Min [MeV]"] = energymin
        sed_df["Energy Center [MeV]"] = energymid
        sed_df["Energy Max [MeV]"] = energymax
        sed_df["Flux [ERG*SM^-2*S^-1]"] = fluxes_axis
        sed_df["Flux Error"] = fluxes_error_axis
        sed_df["Index"] = power_index
        sed_df["Delta Index"] = power_deltaindex
        sed_df["TS"] = ts_list
        sed_df.index = range(1, range_count + 1)
        sed_df.columns.name = "{}".format(title_name)
        sed_df.index.name = "Ranges"
        sed_df.to_csv("{}_sed_data.csv".format(obj_name), sep="\t")
        break

    elif sed_or_lc.lower() in ["lc", "light curve"]:
        # read object properties created previously
        with open("object_properties_lc.txt", "r") as tx_lc:
            properties_lc = json.load(tx_lc)
        range_num = properties_lc.get("range_num")
        tmin_lc = properties_lc.get("tmin")
        tmax_lc = properties_lc.get("tmax")
        obj_name = properties_lc.get("obj_name")
        receiver_email = properties_lc.get("user_email")

        # get number of bins

        lc_bin_counter(tmin_lc, tmax_lc, range_num)
        read_lc_results()
        time_mean = []
        tmin_list = np.array(tmin_list)

        for i in range(len(tmin_list)-1):
            time_mean.append((tmin_list[i] + tmin_list[i + 1]) / 2)
        time_mean = np.array(time_mean)  # middle value of time axis
        delta_time = time_mean - tmin_list[:-1]
        time_mjd = np.round(time_mean/86400+51910, decimals = 1)
        delta_mjd = np.round(delta_time/86400, decimals = 1)
        # plot light curve

        fig, ax = plt.subplots(nrows =3 ,ncols = 1,figsize=(18, 9), sharex = True)
        fontdict = {'fontsize': 18, "fontweight": 'bold'}  # customize axis names
        ax[0].errorbar(time_mjd, fluxes_lc, xerr=delta_mjd, yerr=fluxeserr_lc,
                       fmt=".-", elinewidth=1.002, ls=":", lw=0, color="r", capsize=1.98)  # plot an errorbar of SED
        ax[1].errorbar(time_mjd,power_index_lc,xerr =delta_mjd ,yerr= power_deltaindex_lc,fmt=".-",
                       elinewidth=1.002, ls=":", lw=0,color="g",capsize=1.98)
        ax[2].errorbar(time_mjd, ts_list_lc, xerr=delta_mjd, fmt=".-",
                       elinewidth=1.002, ls=":", lw=0, color="b", capsize=1.98)

        #ax.set_xscale("log")
        #ax.set_yscale("log")
        title_name = obj_name.replace("_", " ")  # change  object name to readable one
        title_name = title_name.upper()
        ax[0].set_title("{title_name} Light Curve".format(title_name=title_name), fontdict=fontdict, loc='center')
        ax[0].set_ylabel("FLUX [PHOTON*SM^-2*S^-1] ", fontdict=fontdict)
        ax[1].set_ylabel("Index", fontdict=fontdict)
        ax[2].set_ylabel("TS", fontdict=fontdict)
        ax[2].set_xlabel("Time, [MJD]", fontdict=fontdict)

        for i in range(3):
            for tick in ax[i].xaxis.get_major_ticks():
                tick.label.set_fontsize(16)

            for tick in ax[i].yaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        fig.subplots_adjust(hspace=0)
        plt.savefig("{obj_name}_Light_Curve.pdf".format(obj_name=obj_name))


        # save results to dataframe
        lc_df = pd.DataFrame()
        lc_df["Time [MET]"] = time_mean
        lc_df["Delta Time"] = delta_time
        lc_df["Flux [ERG*SM^-2*S^-1]"] = fluxes_lc
        lc_df["Flux Error"] = fluxeserr_lc
        lc_df["Index"] = power_index_lc
        lc_df["Delta Index"] = power_deltaindex_lc
        lc_df["TS"] = ts_list_lc
        lc_df.index = range(1, len(tmin_list))
        lc_df.columns.name = "{}".format(title_name)
        lc_df.index.name = "Ranges"
        lc_df.to_csv("{obj_name}_lc_data.csv".format(obj_name=obj_name), sep="\t")
        break

######################################################################################
#                                                                                    #
#                                  ***EMAILING***                                    #
#                                                                                    #
######################################################################################

# e-mailing constants
sender_email = "andreasmessi@mail.ru"  # Enter your address

message = """Dear User, \n
Analysis on Object {} is finally done!""".format(obj_name)
error_message = """Dear User, \n
Analysis stoped, please check working progress!"""
login = "andreasmessi@mail.ru"
password = "Nt8e4f0b#%(!@"

server = smtplib.SMTP("smtp.mail.ru:587")  # establish connection to mail.ru server
server.ehlo()
server.starttls()
server.ehlo()

server.login(login, password)  # login to account of emailing server

# send e-mail
if os.path.exists("{obj_name}_SED.pdf".format(obj_name=obj_name)):
    try:
        # if plot is done send "DONE" email
        msg = MIMEMultipart()
        msg["From"] = sender_email
        msg["To"] = receiver_email
        msg["Subject"] = "SED Analysis is Done(FERMI)!"
        msg.attach(MIMEText(message))
        server.sendmail(sender_email, receiver_email, msg.as_string())
        server.quit()
    except:
        print("Connection to server unexpectedly closed.")
elif os.path.exists("{obj_name}_Light_Curve.pdf".format(obj_name=obj_name)):
    try:
        # if plot is done send "DONE" email
        msg = MIMEMultipart()
        msg["From"] = sender_email
        msg["To"] = receiver_email
        msg["Subject"] = "LightCurve Analysis is Done(FERMI)!"
        msg.attach(MIMEText(message))
        server.sendmail(sender_email, receiver_email, msg.as_string())
        server.quit()
    except:
        print("Connection to server unexpectedly closed.")
else:
    try:
        msg = MIMEMultipart()
        msg["From"] = sender_email
        msg["To"] = receiver_email
        msg["Subject"] = "SED or Light Curve Analysis is Aborted(FERMI)!"
        msg.attach(MIMEText(error_message))
        server.sendmail(sender_email, receiver_email, msg.as_string())
        server.quit()
    except:
        print("Connection to server unexpectedly closed.")

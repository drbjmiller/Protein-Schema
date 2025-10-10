#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============
Generate Additional Figures for Protein Rarity Schema

Written by Brian Miller
===============
"""

############################################################
# Import modules and functions

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as mtick
from mut_functions import apply_styles
from math import (comb, log10, sqrt)
import numpy as np
from scipy.optimize import curve_fit

PROTEIN = "HisA"
FIND_N = False


############################################################
# Linear Equation
def linear(x, a, b):
    return a*x + b

# Exponential decay function
def exp_model(x, a):
    return np.exp(-a*x)

# Hyper-exponential decay function
def hyper_exp(x, a, b):
    return np.exp(-a*x-b*x**2)

############################################################

# Computer number of mutations
if FIND_N == True:
    x = [1, 2, 3, 5, 10]
    muts = [1.16, 3.9, 5, 9.5, 20.7]
    std_muts = [1.4, 1.6, .7, 3.7, 2.8]
    par, cov = curve_fit(linear, x, muts, sigma=std_muts, absolute_sigma=True)  
    a = par[0]
    b = par[1]
    siga = sqrt(cov[0][0])
    sigb = sqrt(cov[1][1])
    param_str = "a = {:.3f} +/- {:.3f}   b = {:.3f} +/- {:.3f}".format(a, siga, b, sigb)
    print(param_str)
    
    mut_list = [4, 6, 7, 8, 9]
    for i in mut_list:
        n = .69*linear(i, a, b)
        print("{:d}: {:.2f}".format(i, n))


# Curve fit hyper-exponential function - exp(-alpha*n - beta*n^2)

if PROTEIN == "Beta":
    mut_num = [0, 1.11, 1.11, 1.11, 2.67, 2.67, 3.46, 3.46, 5.08, 6.55, 6.55, 6.55, 8.02, 8.02, 9.49, 9.49, 11.0, 11.0, 11.0, 12.4, 12.4, 14.3]
    functional = [1, .84, .797, .851, .723, .695, .393, .479, .309, .191, .138, .234, .16, .117, .064, .085, .043, .021, .053, .032, .021, .006]
    parameters, covariance = curve_fit(hyper_exp, mut_num, functional)
    parameters_exp, covariance_exp = curve_fit(exp_model, mut_num, functional)
    y_text = r"$\beta$-lactamase Functional %"
    fig_file = "Beta.png"
elif PROTEIN == "GFP":
    mut_num = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]      
    functional = [1, 0.905745063, 0.8761, 0.7302, 0.5213, 0.3159, 0.1699, 0.0801, 0.0413, 0.0154, .0057]  # Cutoff = 1000     
    std_per = [0.001, 0.008754119, 0.002888557, 0.003941268, 0.005056783, 0.005472745, 0.005511837, 0.005146389, 0.005092847, 0.004609212, 0.004006222]
    parameters, covariance = curve_fit(hyper_exp, mut_num, functional, sigma=std_per, absolute_sigma=True)  
    parameters_exp, covariance_exp = curve_fit(exp_model, mut_num, functional, sigma=std_per, absolute_sigma=True)  
    y_text = "GFP Functional %"
    fig_file = "GFP.png"
elif PROTEIN == "HisA":
    mut_num = [0, 1, 2, 3, 4, 5, 6, 7]
    functional = [1, 0.894, 0.679, 0.477, 0.390, 0.292, 0.222, .158]
    std_per = [.01, 0.0236, 0.0399, 0.0539, 0.0762, 0.0928, 0.139, .158]   # For P(7) = 0, using Clopper–Pearson interval of .453 for STD
    parameters, covariance = curve_fit(hyper_exp, mut_num, functional, sigma=std_per, absolute_sigma=True)  
    parameters_exp, covariance_exp = curve_fit(exp_model, mut_num, functional, sigma=std_per, absolute_sigma=True)  
    y_text = "HisA Functional %"
    fig_file = "HisA.png"
else:
    mut_num = [0, 1, 2, 3, 4, 5]
    functional = [1, 0.815, 0.524, 0.295, 0.292, 0.1620]
    std_per = [0.0001,0.0256, 0.0369, 0.0419, 0.0655, 0.0588]
    parameters, covariance = curve_fit(hyper_exp, mut_num, functional, sigma=std_per, absolute_sigma=True)  
    parameters_exp, covariance_exp = curve_fit(exp_model, mut_num, functional, sigma=std_per, absolute_sigma=True)  
    y_text = "HisA Activity"
    fig_file = "HisAs.png"

alpha = parameters[0]     
beta = parameters[1]
siga = sqrt(covariance[0][0])
sigb = sqrt(covariance[1][1])
param_str = "a = {:.3f} +/- {:.3f}   b = {:.3f} +/- {:.3f}".format(alpha, siga, beta, sigb)
print(param_str)

# Plot P(n)
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(5, 3)
apply_styles(ax, 0, 15, -.01, 1.05, 1, 1, 'n', y_text, True)
ax.yaxis.set_minor_locator(MultipleLocator(.1))

model_x = list(np.linspace(0, 15, num=150))
model_y = list(map(lambda x: hyper_exp(x, alpha, beta), model_x))
plt.plot(model_x, model_y, color='black', linestyle='dotted', linewidth=1, label="model")

# --- NEW: add error bars if available; otherwise plot points as before ---
if PROTEIN == "Beta":
    plt.plot(mut_num, functional, color='blue', marker=".", linestyle="None", label="experiment")    
else:
    yerr = np.array(std_per)  # use 1σ; for ~95% CI use 1.96*np.array(std_per)
    plt.errorbar(
        mut_num, functional,
        yerr=yerr,
        fmt='.',               # point marker only (no connecting line)
        linewidth=1,
        elinewidth=1,
        capsize=3,
        label='experiment (±1σ)'
    )

plt.legend(loc=1)
plt.savefig(fig_file, dpi=300)

# Compute variance explained (R^2) for hyper-exponential model
yhat = []
ss_res = 0
ss_tot = 0
func_mean = sum(functional)/len(functional)
for n in range(len(mut_num)):
    yhat = hyper_exp(mut_num[n], parameters[0], parameters[1])
    ss_res += (functional[n] - yhat)**2
    ss_tot += (functional[n] - func_mean)**2
R2 = 1 - ss_res/ss_tot
print(f"Variance explained (R^2) for hyper-exponential: {R2*100:.2f}%")

# Compute variance explained (R^2) for exponential model
yhat = []
ss_res = 0
ss_tot = 0
func_mean = sum(functional)/len(functional)
for n in range(len(mut_num)):
    yhat = exp_model(mut_num[n], parameters_exp[0])
    ss_res += (functional[n] - yhat)**2
    ss_tot += (functional[n] - func_mean)**2
R2 = 1 - ss_res/ss_tot
print(f"Variance explained (R^2) for exponential: {R2*100:.2f}%")




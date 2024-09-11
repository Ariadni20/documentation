'''It applies the spectrum broadening and prepares the IR and
Raman plots.'''

#AUTHOR: Ariadni Boziki

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import cauchy, norm


def spectrum_broadening(x, x0, y0, distribution, scale, 
			fit_points=True):

	"""It applies the spectrum broadening either using a guassian
	or a lorentzian."""

	if distribution == 'gaussian':
		distribution = norm
	elif distribution == 'lorentzian':
		distribution = cauchy

	s = np.sum([yp * distribution.pdf(x, xp, scale=scale)
			for xp, yp in zip(x0, y0)], axis=0)

	if fit_points:
		s_max = np.max(s)
		if s_max == 0.0:
			s_max = 1.0
		return s * np.max(y0) / s_max


def plot_spectrum_IR(freq, IRintensity, code):

    """It prepares the IR plot."""

    x = np.linspace(-55,4500, num=500, endpoint=True)
    y_IR = spectrum_broadening(x, freq, IRintensity, 'lorentzian', 10.0, fit_points=True)

    np.savetxt("IRspectrum.txt", np.column_stack([x, y_IR]))	

    plt.title("IR spectrum")
#   plt.xlabel(r'Frequency (THz)')
    plt.xlabel(r'Frequency (cm$^{\rm -1}$)')
    if (code == 'aims'):
        plt.ylabel(r'IR Intensity (D$^{\rm 2}$/${\rm \AA}$$^{\rm 2}$ amu)')
    if (code == 'dftb+'):
        plt.ylabel(r'IR Intensity (A.U.)')
    plt.plot(x, y_IR, color='r')
    plt.stem(freq, IRintensity, markerfmt=',', basefmt=" ") 
    plt.savefig("IR_spectrum.png")
    plt.show()


def plot_spectrum_Raman(freq, Ramanactivity, code):

    """It prepares the Raman plot."""

    x = np.linspace(-55,4500, num=500, endpoint=True)
    y_Raman = spectrum_broadening(x, freq, Ramanactivity, 'lorentzian', 10.0, fit_points=True)	

    np.savetxt("Ramanspectrum.txt", np.column_stack([x, y_Raman]))

    plt.title("Raman spectrum") 
#   plt.xlabel(r'Frequency (THz)') 
    plt.xlabel(r'Frequency (cm$^{\rm -1}$)')
    if (code == 'aims'):
        plt.ylabel(r'Raman Activity (${\rm \AA}$$^{\rm 4}$/amu)')
    if (code == 'dftb+'):
        plt.ylabel(r'Raman Activity (A.U.)')
    plt.plot(x, y_Raman, color='r')
    plt.stem(freq, Ramanactivity, markerfmt=',', basefmt=" ") 
    plt.savefig("Raman_spectrum.png")
    plt.show()

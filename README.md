# XPSinv.jl
  this package should contain tools for the estimation of the concentration profiles of a  chemical species across a microjet (i.e. a very small stream) probed with X-rays. The data are the spectra of the emitted photoelectrons from the core subshell in the ground state (i.e. 1s orbital).

The package is divided into two:

- XPSpack: model (deterministic and associated uncertainty)
- XPSinv:  some inversion methods (or these could be taken from another package, and only included in this package)

The repository can be found at [XPSinv](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv)

## Model: XPSpack
Two models are implemented in this module; both are modelling the photoelectric signal generated by photoelectron emitted from C1s in a microjet and collected at an angle by a kinetic energy analyzer.
### Single point spectrum
The PE point model describes each data point collected by a kinetic energy analyzer with solid angular aperture &alpha; and gain T. It takes into account the total photon current F(h&nu;)&sigma;(h&nu;,Ke) and the attenuation of the amount of electrons depending on the depth at which the electron is emitted through an integral term spanning  the depth of the target -- the photon attenuation is not accounted for since it is negligible in the relevant cases. It assumes that the background signal (e.g. fluorescence, inelastic electron scattering, etc) has already been extracted using a baseline removal algorithm (e.g. that implemented by the function baseline_removal in XPSutils.jl). The intensity measured for an photon of energy h&nu; and a electron emitted with kinetic energy Ke is given by [1]:

I(h&nu;,Ke) - I<sub>bg</sub>(h&nu;,Ke) = &alpha;T(&mu;<sub>Ke</sub>) F(h&nu;)&sigma;(h&nu;,Ke) &int;<sub>[0,&infin;[</sub> &rho;<sub>A</sub>(z) e<sup>-g(z)&frasl;&lambda;(Ke)</sup> dz +  &epsilon;(h&nu;,Ke)

where g is the function that gives the relative amount of matter the emitted electron has to go through depending on the depth -- asymptotically, it is equivalent to identity.
The term &epsilon;(h&nu;,Ke) is a stochastic term that describes the measurement noise due to uncontrolled parameters -- it can also include uncertainty due to the deterministic term of the model itself when inverting the data.

The deterministic part of the model is implemented by the functions whose names begin with &Psi; in the XPSmeas.jl file. In particular Ψ_lin_peaks_mean_and_std, implements the linear basis operator and some uncertainty for the deterministic part of the above equation.

In this model, the effect of the analyzer is reduced to two parameters &alpha; and T. As much as the angular aperture is a good model for the integration of the angular density of the isotropic orbital 1s, the factor T might not be as good of a model. Indeed, the mechanisms by which the kinetic energy spectrum is collected, see [2] for an overview of principles, are probably more complex than just the multiplication by the factor T; however, for fine enough resolution and discrimination power, the complexity of the analyzer may be well approximated by a single multiplicative factor.

### Peaks
The spectra collected during an experiment involving organic compounds other that single isolated atom of carbons is usually made up of several overlapping peaks which can be attributed different carbon atoms in the molecules. Several reasons can explain why the signal is not concentrated at a single point in the kinetic energy spectrum:

- chemical shift
- satellites lines
- exchange/relaxation

among others. The interesting feature is that the spread and shift in kinetic energy space can be attributed to the potential energy terms in the Hamiltonian that does not involve the 1s orbital. It is the interaction between the 1s electron and the field created by the other electrons (including the valence electrons, those "binding" with the neighboring atom) that determines the spectral shape of the ionization cross section [3].
NOTE: In this package, we do not aim at creating a predictive model for the cross section spectral shape, but it is estimated from the spectra.

The spectral shape and the different peaks (the number of peak must be known/assumed) can be computed with the functions

- EM_peaks
- cross_section_spread_function and cross_section_spread_function_sample

implemented in the file XPSutils.jl. Later on they will be moved to the inversion module of the package.
man! Am I not creative when it comes to name functions!

### Peak area
An often used quantity  is the peak area, by which it is meant that the spectrum created by a specific neighborhood around a site, e.g. 1s orbital of a C atom surrounded by another carbon, one oxygen and one hydrogen, is a quantity that depends on the concentration of said configuration. It is defined by:

A(h&nu;, &mu;<sub>Ke</sub>) = &int; I(h&nu;,Ke) - I<sub>bg</sub>(h&nu;,Ke) dKe

Note that it is less sensitive to the measurement noise than single point spectra (but there are considerably fewer data points for an algorithm to invert the data once integrated).

### Relative peak area
In some cases, some parameters are not known, e.g. the analyzer characteristics or the photon flux. Another useful quantity is then the ratio of two different compounds taken measured for the same kinetic energy ranges; in this case, the ratio will cancel out the factor &alpha;T.
One very interesting case is the ratio involving C1s and O1s signals. Indeed, the ionization cross section for these atoms are well defined and tabulated, e.g. see [4], and the bulk of the O1s signal can be approximated by the solvent signal (in case of dilute aqueous solutions) whose concentration w.r.t. depth can be approximated by a step function (for O1s, g(z) = z, and &rho;<sub>H2O 1s</sub>(z) = &rho;<sub>0</sub> for z>0 and 0 otherwise).

R<sub>C1s,O1s</sub>(h&nu;<sub>1</sub>, &mu;<sub>Ke<sub>1</sub></sub>,h&nu;<sub>2</sub>, &mu;<sub>Ke<sub>2</sub></sub>) = A<sub>C1s</sub>(h&nu;<sub>1</sub>, &mu;<sub>Ke<sub>1</sub></sub>)/A<sub>O1s</sub>(h&nu;<sub>2</sub>, &mu;<sub>Ke<sub>2</sub></sub>)

NOTE: the peak area and relative area are not implemented as such in the package, but rather derive from the PE signal model and and implemented in examples.

## Inversion: XPSinv

### inversion algorithms

- Null Space Optimization: seems to work for simple cases, but not for all cases
- Chambolle and Pock (CP) Primal-Dual optimization is somewhat more time consuming but the results are more consistent

### Sampling

The data and the measurement model can be sampled so that the uncertainty due to measurement noise and parameter uncertainty can be assessed in the final result, i.e. the concentration profile

## Examples

### data_generation_exp.jl

Test for:

- baseline removal
- peak fitting
- PE model
- data generation (simulation)
- data inversion: some test for PE signal using some version of sampling the data and the model

[link to file](test/data_generation_exp.jl)

### test_all.jl

Test focusing on the inversion of simulated PE signal with all the inversion methods currently implemented

[link to file](test/test_all.jl)

### peak_area_model.jl

Test focusing on the inversion of simulated peak area data

[link to file](test/peak_area_model.jl)

### data_gen_exp_5.jl

Test for:

- baseline removal
- cross section spectral shape (estimate and uncertainty)
- PE and peak area model (and uncertainty)
- simulation of some data

[link to file](test/data_gen_exp_5.jl)

## Refs

- [1] N. Ottosson et al., Photoelectron spectroscopy of liquid water and aqueous solution: Electron effective attenuation lengths and emission-angle anisotropy, Journal of Electron Spectroscopy and Related Phenomena, 2012, Vol. 177, p. 60 ([DOI: 10.1016/j.elspec.2009.08.007](https://www.doi.org/10.1016/j.elspec.2009.08.007))
- [2] D. Roy and D. Tremblay, Design of electron spectrometers, Reports on Progress in Physics, 1990, Vol. 53, No. 12, p. 1621 ([DOI: 10.1088/0034-4885/53/12/003](https://www.doi.org/10.1088/0034-4885/53/12/003))
- [3] S. Manson and J. Cooper, Photo-Ionization in the Soft x-Ray Range: 1 Z Dependence in a Central-Potential Model, Physical Review, 1968, Vol. 165, p. 126 ([DOI: 10.1103/PhysRev.165.126](https://www.doi.org/10.1103/PhysRev.165.126))
- [4] J. Yeh and I. Lindau, Atomic subshell photoionization cross sections and asymmetry parameters: 1⩽ Z⩽ 103, 1985, Vol. 32, No. 1, p. 1--155, ([DOI: 10.1016/0092-640X(85)90016-6](https://www.doi.org/10.1016/0092-640X\(85\)90016-6))


# Dependence

Interpolations
utilsFun
utilsFun

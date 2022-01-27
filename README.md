# SPT_FFTlog_simple
Simple Mathematica notebook to obtain the SPT matter real space power spectrum using FFTLog. 
(The tracers version including RSD will be uploaded soon.)

Alejandro Aviles
avilescervantes@gmail.com

#

Computes P22 and P13 in real space, using the method of  Simonovic et al https://arxiv.org/abs/1708.08130

Just run the file Ploop.m

The main functions are in modules.m. 

#

# Remove BAO from power spectrum

Obtain the non-wiggle linear power spectrum necessary to get the IR-resummed power spectrum:

call: module pknwM[kmin, kmax, Nk, inputPkLinear, hubble] in  file PkNW/pnw.m 

Based in the method of  J. Hamann, S. Hannestad, J. Lesgourgues, C. Rampf and Y. Y. Wong, Cosmological
parameters from large scale structure - geometric versus shape information, JCAP 07 (2010)
022, [arxiv:1003.3999].



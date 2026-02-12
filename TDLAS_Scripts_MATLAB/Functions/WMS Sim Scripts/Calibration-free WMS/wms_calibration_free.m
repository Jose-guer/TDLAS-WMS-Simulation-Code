function [X1f,Y1f,X2f,Y2f,X4f,Y4f,H] = wms_calibration_free(t,tau,i0,psi1,i2,psi2)
% This function returns the X an Y components for the WMS-1f, -2f, and -4f
% signals and Fourier coefficients Hk, using the calibration-free WMS model. 
%
% Inputs
% t - time vector
% tau - transmittance = exp(-absorption)
% i0, i2, - first and second order DC normalized intensity modulation
% amplitudes
% psi1, psi2, - first and second order relative phase-shifts
% 
% Rieker, Gregory B., Jay B. Jeffries, and Ronald K. Hanson. 
% "Calibration-free wavelength-modulation spectroscopy for measurements of 
% gas temperature and concentration in harsh environments." 
% Applied optics 48.29 (2009): 5546-5560.
%
% written by Jose Guerrero, University of Michigan - Aerospace Department
% joseguer@umich.edu


H0 = wms_Hk(t,tau,0);
H1 = wms_Hk(t,tau,1);
H2 = wms_Hk(t,tau,2);
H3 = wms_Hk(t,tau,3);
H4 = wms_Hk(t,tau,4);
H5 = wms_Hk(t,tau,5);
H6 = wms_Hk(t,tau,6);

H = [H1,H2,H3,H4,H5,H6];

X1f = H1 + i0*(H0+H2/2)*cos(psi1) + i2/2*(H1+H3)*cos(psi2);
Y1f = -i0*(H0-H2/2)*sin(psi1) + i2/2*(H1-H3)*sin(psi2) ;

X2f = H2 + i0/2*(H1+H3)*cos(psi1) + i2*(H0 + H4/2)*cos(psi2);
Y2f = -i0/2*(H1-H3)*sin(psi1) + i2*(H0 - H4/2)*sin(psi2);

X4f = H4 + i0/2*(H5+H3)*cos(psi1) + i2/2*(H6 + H2)*cos(psi2);
Y4f = i0/2*(H5-H3)*sin(psi1) + i2/2*(H6 - H2)*sin(psi2);

end
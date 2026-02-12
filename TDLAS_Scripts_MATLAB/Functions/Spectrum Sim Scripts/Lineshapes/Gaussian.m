function [y] = Gaussian(v,vo, wG )
%v  : wavenumber array
%vo : linecenter 
%wG : Doppler linewidth, FWHM

y = 2 / wG * sqrt(log(2)/pi) * exp( -4*log(2) * (v-vo).^2/ wG^2 );

end
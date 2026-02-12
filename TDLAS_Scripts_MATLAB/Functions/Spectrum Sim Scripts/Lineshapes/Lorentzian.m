function [y] = Lorentzian( v,vo, wL )
%v  : wavenumber array
%vo : linecenter 
%wL : Lorentzian linewidth, FWHM

y = 1/(2*pi) * wL ./ ( (v-vo).^2 + (wL/2)^2 );

end
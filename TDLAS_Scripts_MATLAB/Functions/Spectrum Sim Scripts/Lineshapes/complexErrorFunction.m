function [w] = complexErrorFunction(x,y)
% complexErrorFunction  Calculation of complex error function using dimentionless coordinates
% 
% [w] = complexErrorFunction(x,y)   Computes the complex error function
%   using the algorithm developed by Dr. F. Schreier and kindly presented
%   in Fortran. The function was rewriten to MATLAB by Dr. N. Cherkasov
%   For more details on algorithm see the publication:
%   F. Schreier: Optimized Implementations of Rational Approximations for the Voigt and Complex Error Function. 
%   J. Quant. Spectrosc. & Radiat. Transfer, 112(6), 1010–1025, doi 10.1016/j.jqsrt.2010.12.010, 2011. 
%
%   Briefly, the algorithm is compiled from two:
%       for    large x+y     J  Humlicek, JQSRT 27, 437, 1982
%       for small x+y:    J.A.C. Weideman,  SIAM J. Numer. Anal. 31 (1994) pp. 1497-1518,  equation (38.I) and table I
%
% INPUT ARGUMENTS are dimentioneless coordinates x and y
%   x - array 1*N, and y - single variable
%
% OUTPUT
%   w - complex array 1*N
%
% The function was used for the deconvolution of IR spectra
% see the publication
%
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com


half=0.5;
one=1; 
two=2; 
recSqrtPi=1/sqrt(pi);

lenX=length(x);
w=zeros(lenX,1);

%   "Weideman" constants
% n=24;
l=4.1195342878142354; % l=sqrt(n/sqrt(2.))  ! L = 2**(-1/4) * N**(1/2)
a=[-1.5137461654527820e-10,  4.9048215867870488e-09,  1.3310461806370372e-09, -3.0082822811202271e-08, ...
-1.9122258522976932e-08,  1.8738343486619108e-07,  2.5682641346701115e-07, -1.0856475790698251e-06, ...
-3.0388931839840047e-06,  4.1394617248575527e-06,  3.0471066083243790e-05,  2.4331415462641969e-05, ...
-2.0748431511424456e-04, -7.8166429956142650e-04, -4.9364269012806686e-04,  6.2150063629501763e-03, ...
3.3723366855316413e-02,  1.0838723484566792e-01,  2.6549639598807689e-01,  5.3611395357291292e-01, ...
9.2570871385886788e-01,  1.3948196733791203e+00,  1.8562864992055408e+00,  2.1978589365315417e+00];
%   humlicek prbFct region I bounds
s15=15e0;
% -------------------------------------------------------------------------
x12 = y - s15; %        left wing -- center
x21 = -x12;    % 15-y   center -- right wing

if (y>s15 || x(1)>x21 || x(lenX)<x12)
%       all points are in Humlicek's region I
    for ii=1:lenX
        t= y - x(ii)*1i;
        w(ii) = (recSqrtPi*t) / (half + t*t);
    end
else
    for ii=1:lenX
        s  = abs(x(ii)) + y;
        if (s>s15)
            t     = y-x(ii)*1i;
            w(ii) = (recSqrtPi*t) / (half + t*t);
        else
            recLmZ  = one / (l+y-x(ii)*1i);
            t       = (l-y+x(ii)*1i) * recLmZ;
            w(ii) =  recLmZ  *  (recSqrtPi + two*recLmZ*...
                (a(24)+(a(23)+(a(22)+(a(21)+(a(20)+(a(19)+(a(18)+(a(17)+(a(16)+(a(15)+(a(14)+(a(13)+(a(12)+(a(11)+(a(10)+(a(9)+...
                (a(8)+(a(7)+(a(6)+(a(5)+(a(4)+(a(3)+(a(2)+a(1)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t));
        end
    end
end
    
end


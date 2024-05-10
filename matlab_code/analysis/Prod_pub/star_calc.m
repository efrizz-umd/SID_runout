% star_calc.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [Estar,Rstar,mstar,Gstar] = star_calc(E,R,nu,m)

% ************************************************************************
% This function computes equivalent quantities. See associated work 1
%
% % ----------- output ----------- %
% - Estar - equivalent Elastic modulus (Pa)
% - Rstar - equivalent radius (m)
% - mstar - equivalent mass (kg)
% - Gstar - equivalent shear modulus (Pa)
% % ----------- intput ----------- %
% - E - particle Elastic modulus (Pa)
% - R - particle radius (m)
% - nu - poisson ratio
% - m - particle mass (kg)
%
% ************************************************************************

%   Compute reduced quantities for a two particle collision where the
%   particles are identical

Rstar = R.*R./(R+R);

Estar = ((1-nu^2)./E + (1-nu^2)./E ).^(-1);

mstar = (1./m + 1./m).^(-1);

Gstar = ( 2*(2+nu)*(1-nu)./E + 2*(2+nu)*(1-nu)./E ).^(-1);

end

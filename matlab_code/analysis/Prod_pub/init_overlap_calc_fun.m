% init_overlap_calc_fun.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [init_overlap,init_overlap_cohesion,Fm_withoverlap,Fm_withoutoverlap] =  ...
  init_overlap_calc_fun(rho_p,phi,z,R,g,E,kc,nu,v0)

% ************************************************************************
% This function calculates the predicted initial overlap, with and without cohesion
%
% % ----------- output ----------- %
% - init_overlap - initial overlap vs depth (m)
% - init_overlap_calc_cohesion - initial overlap vs depth with cohesion (m)
% - Fm_withoverlap - initial force on particle as a result of impact vs depth (m)
% - Fm_withoutoverlap - initial force on particle as a result of impact vs depth (m)
% % ----------- intput ----------- %
% - rho_p - particle density (kg/m3)
% - phi - packing fraction
% - z - a depth (m)
% - R - particle radius (m)
% - g - gravity (m/s2)
% - E - elastic modulus (Pa)
% - kc - cohesion energy density (Pa)
% - nu - poisson ratio
% - v0 - impact velocity (m/s)
%
% ************************************************************************

F0 = rho_p*phi*z*pi*g*R^2;

init_overlap = 2*(F0*3*(1-nu^2)/(4*E)).^(2/3) / R^(1/3);

vol = 4/3 * pi * R^3;
m = rho_p*vol;

[Estar,Rstar,~,~] = star_calc(E,R,nu,m);
kn = (4/3)*Estar*sqrt(Rstar);
%keyboard

L = length(z);
init_overlap_cohesion = zeros(1,L);

% for jj = 1:L
%     syms d0
%
% %     kn = (4/3)*Estar*sqrt(Rstar);
%     Ac = -pi/4 * (d0^2 - 4*R*d0);
%     eqn = F0(jj) - norm(-kn*d0^(3/2) + kc*Ac);
%
%     d0solve = vpasolve(eqn == 0,d0);
%
%     % sometimes the solver returns a very small imaginary portion (like
%     % 1E-30), we are just going to overwrite this
%     if imag(d0solve)
%
%         if abs(imag(d0solve)) < eps
%
%             d0solve = real(d0solve);
%         else
%             d0solve = NaN;
%         end
%     end
%
%     init_overlap_cohesion(1,jj) = d0solve;
% end


Fm_withoverlap = kn*((5/4)*m*(v0^2)/kn + init_overlap.^(5/2)).^(3/5);
Fm_withoutoverlap = kn*((5/4)*m*(v0^2)/kn).^(3/5);

end

% FSW_impulse_fun_v4.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [Fm_vs_depth,Fm_vs_depth_scaled,Fm_pal, ...
    Fm_avg, Fm_std, Fm_scaled_avg, Fm_scaled_std,Fm_vs_depth_cohesive,Fm_vs_depth_c_analytical] = ...
    FSW_impulse_fun_v4(v0,delta0,SF,E_p,r_p,rho_p,depthplot,phi_init,kc)

% ************************************************************************
% This function computes the predicted force in the wave front as a function of
% impact velocity, initial overlap, elastic modulus, particle size, and density.
% See associated work 2 for details.
%
% % ----------- output ----------- %
% Fm is the maximum force in the wave front at each depth. There are several
% different varieties, the one used in the output plot (average and standard dev.)
% is Fm_vs_depth
% % ----------- intput ----------- %
% - v0 - piston velocity  (m/s)
% - delta0 - initial overlap (m)
% - SF - scaling factor
% - E_p - particle elastic modulus (Pa)
% - r_p - particle radius (m)
% - rho_p - particle density (kg/m3)
% - depthplot - array of sensor positions
% - phi_init - initial packing fraction
% - kc - cohesion energy density (Pa)
% ************************************************************************


% particle volume
vol_p = (4/3)*pi*r_p^3;
% mass of particle
m_p= vol_p*rho_p;
% hardcoded, adapt as needed
nu = 0.2;

% compute equivalent quantities
[Estar,Rstar,~,~] = star_calc(E_p,r_p,nu,m_p);

% spring constant
kn = (4/3)*Estar*sqrt(Rstar);



%SF = 1/50; % scaling factor, determined from sim (peak compression force vs SW force)

%Fm = ((5/4)*m_p^3*vimp^6*kn^2 + kn^5*d0.^(15/2)).^(3/5);

%% vs depth
% my derivation without cohesion
C1 = 0.0073;
Fm_vs_depth = C1*kn*( (5/(4))*m_p*(v0.^2)/kn + delta0.^(5/2)).^(3/5);


% kn = kn
% m_p = m_p
% v0 = v0
% delta0_p = delta0(10)

Fm_vs_depth_scaled = Fm_vs_depth * SF;

%% scalar
% pal derivation
Fm_pal = 0.719*(m_p^3 * Estar^2 * Rstar * v0^6).^(1/5);

Fm_avg = mean(Fm_vs_depth);
Fm_std = std(Fm_vs_depth);

Fm_scaled_avg = mean(Fm_vs_depth_scaled);
Fm_scaled_std = std(Fm_vs_depth_scaled);

%% with cohesion
syms Fms

L = length(delta0);
Fm_vs_depth_cohesive = zeros(L,1);

kc = 1000; % hardcode for now, make variable later
k = kn;

% for jj = 1:L
%
%     d0 = delta0(jj);
%
%     eqns = (2/5)*k*((Fms/k).^(5/3) - d0^(5/2)) + ...
%         (pi*kc/4)*( (1/3)*((Fms/k).^(2/3)-d0).^3 - 2*((Fms/k).^(2/3)-d0).^2*r_p ...
%         +((Fms/k).^(2/3)-d0)*8*r_p^2  ) - (1/2)*m_p*v0.^2;
%
%     Fm_vs_depth_cohesive(jj,1) = vpasolve(eqns == 0);
% end


%% all analytical prediction
z = depthplot/100; % m
theta = 3*((1-nu^2)/(4*E_p));
gL = 1.625; % m/s2 lunar gravity
L_c = length(z);
Fm_vs_depth_c_analytical = zeros(L_c,1);

%kc = 1000; % hardcode for now, make variable later
k = kn;

rho_b = rho_p*phi_init;
d0_computed = 2*(theta*gL*rho_b*pi*r_p^2*z.^(2/3)) / r_p^(1/3);

% for jj = 1:L_c
%
%
%
%     d0 = d0_computed(jj);
%
%     eqns = (2/5)*k*((Fms/k).^(5/3) - d0^(5/2)) + ...
%         (pi*kc/4)*( (1/3)*((Fms/k).^(2/3)-d0).^3 - 2*((Fms/k).^(2/3)-d0).^2*r_p ...
%         +((Fms/k).^(2/3)-d0)*8*r_p^2  ) - (1/2)*m_p*v0.^2;
%
%     Fm_vs_depth_c_analytical(jj,1) = vpasolve(eqns == 0);
% end

end

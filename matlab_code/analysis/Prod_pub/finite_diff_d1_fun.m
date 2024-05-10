% finite_diff_d1_fun.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [deriv] = finite_diff_d1_fun(u,dx)

% ************************************************************************
% This function calclates a derivative using a finite difference method
%
% % ----------- output ----------- %
% - u - position (or whatever you take the derivative of)
% - dx - grid spacing
% % ----------- intput ----------- %
% - deriv - velocity (or whatever computed derivative you found)
%
% ************************************************************************


    deriv = zeros(1,length(u));

    % use a forward difference method on LHS (first order)
    deriv(1,1) = (u(2)-u(1))/dx;

    % use backward difference on the RHS (first order)
    deriv(1,end) = (u(end) - u(end-1))/dx;

    % use central difference (second order) on RHS - 1 and LHS + 1
    deriv(1,2) = (u(3) - u(1))/(2*dx);
    deriv(1,end-1) = (u(end) - u(end-2))/(2*dx);

    % use central difference (fourth order) on interior points
    deriv(1,3:(end-2)) = (-1*u(5:end) + 8*u(4:(end-1)) - 8*u(2:(end-3)) + u(1:(end-4)))/(12*dx);

end

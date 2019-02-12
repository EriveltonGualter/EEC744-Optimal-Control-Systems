% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
% Problem 1-1
%
% Erivelton Gualter, 01/22/2018

close all

A = [ -0.16 0 ; 0.16 -0.16];

% State Transition Matrix
syms s
STM = ilaplace(inv(s*eye(size(A))-A)); % State Transition Matrix
pretty(STM)

% q(t) and p(t)

X0 = [60; 0];
xt = STM*X0;

qt = xt(1);
pt = xt(2);

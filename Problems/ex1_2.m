% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
% Problem 1-2
%
% Erivelton Gualter, 01/22/2018

close all

R = 3;
L = 1;
C = 1/2;

A = [-R/L -1/L; 1/C 0];
B = [1/L; 0];

% State Transition Matrix
syms s
STMs = inv(s*eye(size(A))-A);
STMt = ilaplace(STMs)            % State Transition Matrix

% q(t) and p(t)

X0 = [0; 0];
U= 2/s;

xt = ilaplace(STMs*X0 + STMs*B*U)

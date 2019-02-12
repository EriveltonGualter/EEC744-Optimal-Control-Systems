clear all, clc

syms rtf1 rtf2 xtf1 xtf2 c11 c12 c21 c22 h11 h12 h21 h22

T = 10;
dt = 0.02;
t = 0:dt:T; 
r = [sin(t); t/2];

C = [c11 c12; c21 c22];
H = [h11 h12; h21 h22];
rtf = [r(1,end); r(2,end)];
xtf = [xtf1; xtf2];

f = transpose(C*xtf-rtf)*H*(C*xtf-rtf) - 10*xtf1^2;

pretty(f)

sol = solve(f,C)

h11 = 20;
h12 =  0;
h21 =  0;
h22 =  0;

sol
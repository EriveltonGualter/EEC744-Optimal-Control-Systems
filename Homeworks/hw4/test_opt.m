
clear all
clc

W0 = ones(1,7);
%HQRC

A = [];
b = [];

Aeq = [];
Beq = [];

LB = zeros(1,7);
UB = 100*ones(1,7);

FUN = @(W) opt(W);

[x,~,flag,~] = fmincon(FUN,W0,A,b,Aeq,Beq,LB,UB)
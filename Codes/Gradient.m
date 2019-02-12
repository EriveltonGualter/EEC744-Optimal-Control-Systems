function [u, SizeHuArr, JArr] = Gradient(q, r)

% Gradient descent optimization of a second-order Newtonian system.
% Analytical calculation of derivatives.
% Try with q = 1, and r = 0.01, 100, and 1000.
% Try with different values of eps.
% Try with different initial u guesses.

eps = 0.02;
eps = 0.002; % gradient descent step size
%eps = 0.0002;
if ~exist('q', 'var')
    q = 1;
end
if ~exist('r', 'var')
    r = 1;
end

tf = 10; % final time
dt = 0.01; % time step size
t = 0 : dt : tf;
N = size(t, 2); % number of time steps

% Define the dynamic system.
A = [0 1; 0 0];
B = [0; 1];
C = [1 0; 0 1];
D = [0; 0];
sys = ss(A, B, C, D);

% Guess an initial optimal control.
%u = 2*ones(size(t));
u = zeros(size(t));

NormHu = inf;
NormHuArr = [];
JArr = [];
for iter = 1 : inf
   [x, t] = lsim(sys, u, t); % system simulation
   % Compute the costate
   p = zeros(N, 2); 
   p(N, 1) = 2 * q * (x(N, 1) - 100);
   p(N, 2) = 0;
   for i = N-1 : -1 : 1
      pdot = [0 -p(i+1, 1)];
      p(i, :) = p(i+1, :) - dt * pdot;
   end
   % Compute the Hamiltonian
   for i = 1 : N
      H(i) = r * u(i)^2 + p(i, :) * [x(i, 2); u(i)];
   end
   % Compute the partial of the Hamiltonian with respect to the control
   for i = 1 : N
      Hu(i) = 2 * r * u(i) + p(i, 2);
   end
   % Compute the cost function.
   J = q * (x(N, 1) - 100)^2 + r * dt * trapz(u.^2);
   JArr = [JArr J];
   % Compute the norm of Hu. Quit when Hu gets close to zero.
   NormHuOld = NormHu;
   NormHu = dt * trapz(Hu.^2);
   disp(['Iteration # ',num2str(iter),', Hu = ',num2str(NormHu)]);
   if (NormHu > NormHuOld) break, end
   NormHuArr = [NormHuArr NormHu];
   if (NormHu < 0.01) break, end
   % Update the control using gradient descent.
   u = u - eps * Hu;
end
disp(['Cost = ',num2str(J)]);
% Plot the results.
close all;
figure;
plot(JArr);
title('Cost Function');
xlabel('Iteration Number');
figure;
plot(t, u);
title('Optimal Control');
xlabel('Time');
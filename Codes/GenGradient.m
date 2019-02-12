function [u, SizeHuArr, JArr] = GenGradient(q, r)

% Generalized gradient descent optimization of a second-order Newtonian system.
% Try with q = 1, and r = 0.01, 100, and 1000.
% Try with different values of eps.
% Try with different initial u guesses.

epsInit = 0.004; % initial gradient descent step size
if ~exist('q', 'var')
    q = 1;
end
if ~exist('r', 'var')
    r = 1;
end

tf = 10; % final time
dt = tf / 100; % time step size
t = 0 : dt : tf;
N = size(t, 2);

% Define the dynamic system.
A = [0 1; 0 0];
B = [0; 1];
C = [1 0; 0 1];
D = [0; 0];
sys = ss(A, B, C, D);

% Guess an initial optimal control.
% We have 3 control histories - one for u0, one for u1, and one for u2.
u = 2*ones(3, size(t, 2));
u = zeros(3, size(t, 2));

NormHu = inf;
NormHuArr = [];
JArr = [];
% We have 3 cost functions - one for u0, one for u1, and one for u2.
J = inf * [1 1 1];
for iter = 1 : inf
   % Simulate the system with the given control.
   [x, t] = lsim(sys, u(1, :), t);
   % Compute the costate.
   p = zeros(N, 2);
   p(N, 1) = 2 * q * (x(N, 1) - 100);
   p(N, 2) = 0;
   for i = N-1 : -1 : 1
      pDot = [0 -p(i+1, 1)];
      p(i, :) = p(i+1, :) - dt * pDot;
   end
   % Compute the partial of the Hamiltonian with respect to the control.
   for i = 1 : N
      Hu(i) = 2 * r * u(1, i) + p(i, 2);
   end
   NormHuArr = [NormHuArr trapz(Hu.^2)];
   % Compute the cost at u(1, :).
   JOld = J(1);
   J(1) = q * (x(N, 1) - 100)^2 + r * dt * trapz(u(1,:).^2);
   disp(['Iteration # ',num2str(iter),', Cost = ',num2str(J(1))]);
   JArr = [JArr J(1)];
   % Quit when the cost function stops decreasing.
   if abs((J(1) - JOld)/J(1)) < 1e-5, break, end
   % Update the control using gradient descent.
   J(2) = inf;
   eps = epsInit;
   while J(2) >= J(1)
        eps = eps / 2;
        u(2, :) = u(1, :) - eps * Hu;
        [x, t] = lsim(sys, u(2, :), t);
        J(2) = q * (x(N, 1) - 100)^2 + r * dt * trapz(u(2,:).^2);
   end
   % Now that eps has been computed, calculate the third control and the
   % corresponding cost.
   u(3, :) = u(1, :) - 2 * eps * Hu;
   [x, t] = lsim(sys, u(3, :), t);
   J(3) = q * (x(N, 1) - 100)^2 + r * dt * trapz(u(3,:).^2);
   % Find the "optimal" value of eps by fitting a quadratic to the 
   % three cost values.
   A = [0 0 1;
        eps^2 eps 1;
        4*eps^2 2*eps 1];
   K = inv(A) * J';
   eps = -K(2) / 2 / K(1);  
   % Update the control.
   u(1, :) = u(1, :) - eps * Hu;
end
% Plot the results.
close all;
figure;
plot(JArr);
title('Cost Function');
xlabel('Iteration Number');
figure;
plot(t, u(1,:));
title('Optimal Control');
xlabel('Time');
function [u, SizeHuArr, JArr] = VariationExtremals(q1, yf, q2, zf, r)

% Thrust direction programming solution using variation of extremals.
if ~exist('q1', 'var')
    q1 = 1;
end
if ~exist('yf', 'var')
    yf = 1;
end
if ~exist('q2', 'var')
    q2 = 1;
end
if ~exist('zf', 'var');
    zf = 1;
end
if ~exist('r', 'var')
    r = 5;
end

% initial state
x0 = [0; 0; 0; 0];
%x0 = [0; 1; 0; 0];
%x0 = [0; 0; 0; 1];

tf = 10; % final time
dt = 0.01; % integration step size
tArray = 0 : dt : tf;
N = size(tArray, 2); % number of time steps

p0 = [1; 1; 1; 1]; % initial guess for p(0)
for iter = 1 : inf
    % state and costate integration
    x = x0;
    p = p0;
    for t = 1 : N
        beta = atan2(p(4), p(2)); % obtained from Hu = 0
        Thrust = -(p(2) * cos(beta) + p(4) * sin(beta)) / 2 / r;
        xdot = [x(2); Thrust*cos(beta); x(4); Thrust*sin(beta)];
        pdot = [0; -p(1); 0; -p(3)];
        x = x + xdot * dt;
        p = p + pdot * dt;
        betaArray(t) = beta;
        ThrustArray(t) = Thrust;
        zArray(t) = x(1);
        yArray(t) = x(3);
    end
    % Compute the cost function.
    J = q1 * (x(3) - yf)^2 + q2 * (x(1) - zf)^2 + r * dt * trapz(ThrustArray.^2);
    JArr(iter) = J;
    % Check for convergence
    dhdx = [2*q2*(x(1)-zf); 0; 2*q1*(x(3)-yf); 0];
    ErrorNorm = norm(p-dhdx, 2);
    disp(['Cost = ', num2str(J), ', Error norm = ', num2str(ErrorNorm)]);
    if ErrorNorm < 0.01, break, end
    % integrate the Px and Pp matrices
    Px = zeros(4);
    Pp = eye(4);
    d2Hdx2 = zeros(4);
    d2Hdp2 = [0 0 0 0; 0 -1/2/r 0 0; 0 0 0 0; 0 0 0 -1/2/r];
    d2Hdpdx = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
    d2Hdxdp = d2Hdpdx';
    for t = 1 : N
        Pxdot = d2Hdpdx * Px + d2Hdp2 * Pp;
        Ppdot = -d2Hdx2 * Px - d2Hdxdp * Pp;
        Px = Px + Pxdot * dt;
        Pp = Pp + Ppdot * dt;
    end
    % Newton's method update
    d2hdx2 = [2*q2 0 0 0; 0 0 0 0; 0 0 2*q1 0; 0 0 0 0];
    p0 = p0 + inv(d2hdx2 * Px - Pp) * (p - dhdx); % updated p(0) guess
end
% Plot the results.
close all;
figure;
plot(JArr);
title('Cost Function');
xlabel('Iteration Number'); ylabel('Cost');

figure;
plot(tArray, ThrustArray);
title('Optimal thrust');
xlabel('Time'); ylabel('Thrust');

figure;
plot(tArray, betaArray);
title('Optimal thrust angle');
xlabel('Time'); ylabel('Angle (rad)');
figure;

plot(zArray, yArray);
title('Optimal trajectory');
xlabel('horizontal'); ylabel('vertical');
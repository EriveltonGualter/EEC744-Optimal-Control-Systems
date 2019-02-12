% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 03/26/2018
% Problem 6.34

clc;
clear all; 
close all;

% Simulation Parameters
tf = 0.78;
dt = 0.01; % time step
t = 0 : dt : tf;
N = length(t); % number of time steps
X0  = [0.05; 0];      

% Parameters
R = 0.1;
eps = 0.01;

u = zeros(size(t));

NormHu = inf;
JArr = [];

for iter = 1 : Inf
    % Simulatio system
    X = X0;
    for k=1 : length(t)-1
        x1 = X(1,k);
        x2 = X(2,k);

        xd1 = -2*(x1+0.25) + (x2+0.5)*exp(25*x1/(x1+2)) -(x1+0.25)*u(k);
        xd2 = 0.5 - x2 - (x2+0.5)*exp(25*x1/(x1+2));

        XDOT = [xd1; xd2];

        X(:,k+1) = X(:,k) + XDOT*dt;
    end

    % Compute the costate
    p = zeros(2, N); 
    p(1, N) = 0;
    p(2, N) = 0;
    
    for i = N-1 : -1 : 1
        x1 = X(1,i+1);
        x2 = X(2,i+1);
        p1 = p(1,i+1);
        p2 = p(2,i+1);

        pDot(1, i) = -2*x1 + 2*p1 - ...
                  p1*(x2+0.5)*(50/(x1+2)^2)*exp(25*x1/(x1+2)) + ...
                  p1*u(i+1) + p2*(x2+0.5)*(50/(x1+2)^2)*exp(25*x1/(x1+2));
        pDot(2, i) = -2*x2 - p1*exp(25*x1/(x1+2)) + p2*(1 + exp(25*x1/(x1+2)));

        p(:, i) = p(:, i+1) - dt * pDot(:, i);
    end

    % Compute cost
    J = dt*trapz(X(1,:).^2 + X(2,:).^2 + R*u.^2);
    JArr = [JArr J];
       
    % Compute the partial of the Hamiltonian with respect to the control
    for i=1:N
        Hu(i) = 2*R*u(i)-p(1,i)*(X(1,i)+0.25);
    end
    
    NormHu = sqrt(dt * trapz(Hu.^2));
    disp(['Iteration # ',num2str(iter),', Hu = ',num2str(NormHu)]);
    if NormHu < 0.01
        break;
    end
    u = u - eps * Hu;
end

%%
close all
f1 = figure;
subplot(211); plot(t, X(1,:), '-k','LineWidth',2)
    title('Steepest Descent Solution', 'Interpreter','Latex', 'FontSize',14);
    ylabel('$x_1$', 'Interpreter','Latex', 'FontSize',14);
subplot(212); plot(t, X(2,:),'-k','LineWidth',2)
    ylabel('$x_2$', 'Interpreter','Latex', 'FontSize',14);
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);

f2 = figure;
plot(JArr,'xk','LineWidth',2);
    ylabel('Cost', 'Interpreter','Latex', 'FontSize',14);
    xlabel('Interaction Number', 'Interpreter','Latex', 'FontSize',14);

f3 = figure;
plot(t,u, '-k','LineWidth',2);
    title('Optimal Control', 'Interpreter','Latex', 'FontSize',14);
    ylabel('$u(t)$', 'Interpreter','Latex', 'FontSize',14);
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    
saveFigureToPdf('fig1',f1);
saveFigureToPdf('fig2',f2);
saveFigureToPdf('fig3',f3);
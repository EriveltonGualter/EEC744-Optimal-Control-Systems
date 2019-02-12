% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 03/26/2018
% Problem 6.35

clc;
clear all; 
close all;

% Simulation Parameters
tf = 0.78;
dt = 0.01; % time step
t = 0 : dt : tf;
N = length(t); % number of time steps

% Parameters
R = 0.1;

P0 = [1; 0.5];
% P0=[1.0782; 0.1918];
X0  = [0.05; 0];

JArr = [];
count = 0;
for iter = 1 : Inf

    % Simulatio system
    X = X0;
    P = P0;
    for k=1 : N-1
        x1 = X(1,k);
        x2 = X(2,k);
        p1 = P(1,k);
        p2 = P(2,k);

        u(k) = p1*(x1+0.25)/(2*R);
        
        xd1 = -2*(x1+0.25) + (x2+0.5)*exp(25*x1/(x1+2)) -(x1+0.25)*u(k);
        xd2 = 0.5 - x2 - (x2+0.5)*exp(25*x1/(x1+2));        
        pd1 = -2*x1 + 2*p1 - ...
                  p1*(x2+0.5)*(50/(x1+2)^2)*exp(25*x1/(x1+2)) + ...
                  p1*u(k) + p2*(x2+0.5)*(50/(x1+2)^2)*exp(25*x1/(x1+2));
        pd2 = -2*x2 - p1*exp(25*x1/(x1+2)) + p2*(1 + exp(25*x1/(x1+2)));
        
        XDOT = [xd1; xd2];
        PDOT = [pd1; pd2]; 

        X(:,k+1) = X(:,k) + XDOT*dt;
        P(:,k+1) = P(:,k) + PDOT*dt;
    end

    J = dt*trapz(X(1,1:end-1).^2 + X(2,1:end-1).^2 + R*u.^2);    
    JArr = [JArr J];
count = count + 1
    p = P(:,end);
    ErrorNorm = norm(p);
    if ErrorNorm < 0.001
        break
    end
    if (norm(P(1,end)) + norm(P(2,end))) < 1e-5
        break;
    end
    disp(['Cost = ', num2str(J), ', Error norm = ', num2str(ErrorNorm)]);

    
    % integrate the Px and Pp matrices
    Px = zeros(2,2);
    Pp = eye(2);
        
    for i = 1 : N
        x1 = X(1,i);
        x2 = X(2,i);
        
        alpha = exp(25*x1/(x1+2));
        d2Hdpdx = [-2+50*(x2+0.5)*alpha/(x1+2)^2-(x1+0.25)*p1/R, alpha; ...
              -50*(x2+0.5)*alpha/(x1+2)^2, -1-alpha];
        d2Hdp2 = [-(x1+0.25)^2/(2*R), 0; 0, 0];
        d2Hdx2 = [-2+(p2-p1)*(100*(23-x1)*(x2+0.5)/(x1+2)^4)*alpha + p1^2/(2*R), 50*(p2-p1)*alpha/(x1+2)^2; ...
                50*(p2-p1)*alpha/(x1+2)^2,  -2];
        d2Hdxdp = d2Hdpdx';
    
        Pxdot = d2Hdpdx * Px + d2Hdp2 * Pp;
        Ppdot = d2Hdx2 * Px - d2Hdxdp * Pp;
        Px = Px + Pxdot * dt;
        Pp = Pp + Ppdot * dt;
    end
    P0 = P0 - inv(Pp) * p % updated p(0) guess
end

%%
close all
f4 = figure;
subplot(211); plot(t, X(1,:), '-k','LineWidth',2)
    title('Variation of Extremal Solution', 'Interpreter','Latex', 'FontSize',14);
    ylabel('$x_1$', 'Interpreter','Latex', 'FontSize',14);
subplot(212); plot(t, X(2,:),'-k','LineWidth',2)
    ylabel('$x_2$', 'Interpreter','Latex', 'FontSize',14);
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);

f5 = figure;
plot(JArr,'xk','LineWidth',2);
    ylabel('Cost', 'Interpreter','Latex', 'FontSize',14);
    xlabel('Interaction Number', 'Interpreter','Latex', 'FontSize',14);

f6 = figure;
plot(t(1:end-1),u, '-k','LineWidth',2);
    title('Optimal Control', 'Interpreter','Latex', 'FontSize',14);
    ylabel('$u(t)$', 'Interpreter','Latex', 'FontSize',14);
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    
saveFigureToPdf('fig4',f4);
saveFigureToPdf('fig5',f5);
saveFigureToPdf('fig6',f6);
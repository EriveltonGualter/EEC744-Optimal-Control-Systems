% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
% Problem 1-3
%
% Erivelton Gualter, 01/22/2018

close all

% Mechanical System Description
M = 1; % kg
K = 2; % N/m
B = 2; % N/m/s

A = [ 0 1; -K/M -B/M];
B = [0 ; 1/M];

% State Transition Matrix
syms s
STM = ilaplace(inv(s*eye(size(A))-A)) % State Transition Matrix

X0 = [0.2;0];	% Initial Condition 

% Numerically simulate the system
Ts = 1e-1;  % Sample Time
Tf = 10;    % Final time
solver = 'ode4';    % Runge-Kutta Method

Options = simset('Solver', solver, 'FixedStep', Ts);
sim('ex1_3sim.slx', [0 Tf], Options);

% Analytical solutions

Ts = 1e-3;  % Sample Time
t = 0:Ts:Tf;
xt = ilaplace(inv(s*eye(size(A))-A)*X0 + inv(s*eye(size(A))-A)*B*2/(s+2))
xta =  [exp(-2*t) - (4*exp(-t).*(cos(t) - (3*sin(t))/2))/5 ;
       2*exp(-t).*(cos(t) - sin(t)/5) - 2*exp(-2*t) ];

%% 
close all
subplot(311); plot(tsim, xsim,'LineWidth',2); 
    legend('x1 - Numerical','x2 - Numerical');
    title('Numerical Solution','Interpreter','latex','FontSize',14); 
subplot(312); plot(t, xta,'LineWidth',2); 
    legend('x1 - Analytical','x2 - Analytical'); 
    title('Analytical Solution','Interpreter','latex','FontSize',14); 
subplot(313); hold on; plot(t, xta,'LineWidth',2); plot(tsim, xsim,'LineWidth',1); 
    legend('x1 - Numerical','x2 - Numerical','x1 - Analytical','x2 - Analytical');
    title('Numerical and Analytical Solution','Interpreter','latex','FontSize',14); 
    xlabel('Time [s] ','Interpreter','latex','FontSize',14); 
    
%% Others Numerically simulations of the system
Ts = 1e-1;  % Sample Time
Tf = 10;    % Final time
close all
Options = simset('Solver', 'ode1', 'FixedStep', Ts);
sim('ex1_3sim.slx', [0 Tf], Options);
hold on;
ax1 = subplot(211); plot(ax1, tsim, xsim(:,1),'LineWidth',2); 
ax2 = subplot(212); plot(ax2, tsim, xsim(:,2),'LineWidth',2); 

Options = simset('Solver', 'ode4', 'FixedStep', Ts);
sim('ex1_3sim.slx', [0 Tf], Options);
hold(ax1,'on'); plot(ax1, tsim, xsim(:,1),'LineWidth',2); hold(ax1,'off')
hold(ax2,'on'); plot(ax2, tsim, xsim(:,2),'LineWidth',2); hold(ax2,'off')

Options = simset('Solver', 'ode45', 'FixedStep', Ts);
sim('ex1_3sim.slx', [0 Tf], Options);
hold(ax1,'on'); plot(ax1, tsim, xsim(:,1),'LineWidth',2); hold(ax1,'off')
hold(ax2,'on'); plot(ax2, tsim, xsim(:,2),'LineWidth',2); hold(ax2,'off')

hold(ax1,'on'); plot(ax1, t, xta(1,:),'LineWidth',2); hold(ax1,'off')
hold(ax2,'on'); plot(ax2, t, xta(2,:),'LineWidth',2); hold(ax2,'off')
    
legend(ax1, 'Eulers Method', ...
            'Runge-Kutta 4th Order Method',...
            'Runge-Kutta, Dormand-Prince (4,5) pair', ...
            'Analitical Method');
legend(ax2, 'Eulers Method', ...
            'Runge-Kutta 4th Order Method',...
            'Runge-Kutta, Dormand-Prince (4,5) pair', ...
            'Analitical Method');
        
title(ax1, 'State $x_1$','Interpreter','latex','FontSize',14); 
title(ax2, 'State $x_2$','Interpreter','latex','FontSize',14); 

xlabel(ax1, 'Time [s] ','Interpreter','latex','FontSize',14); 
xlabel(ax2, 'Time [s] ','Interpreter','latex','FontSize',14); 
ax = [ax1; ax2];
linkaxes(ax,'x');
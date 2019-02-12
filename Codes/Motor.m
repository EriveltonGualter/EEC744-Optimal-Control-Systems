function Motor(ControlChange)

% Two-phase step motor simulation.
% This file can be used to investigate the effects of linearization.
% If ControlChange = 0 then the nonlinear and linearized simulations should
% match exactly (in steady state). As ControlChange increases, the linearized simulation
% diverges more and more from the nonlinear simulation.

if ~exist('ControlChange', 'var')
    ControlChange = 0;
end

Ra = 1.9; % Winding resistance
L = 0.003; % Winding inductance
lambda = 0.1; % Motor constant
J = 0.00018; % Moment of inertia
B = 0.001; % Coefficient of viscous friction

dt = 0.0005; % Integration step size
tf = 4; % Simulation length 

x = [0; 0; 0; 0]; % Initial state
xlin = x; % Linearized approximation of state
dx = x - xlin; % Difference between true state and linearized state
w = 2 * pi; % Control input frequency

% Initialize arrays for plotting at the end of the program
N = round(tf/dt) + 1;
xArray = zeros(4, N);
xlinArray = zeros(4, N);

i = 0;
% Begin simulation loop
for t = 0 : dt : tf
    i = i + 1;
    % Save data for plotting
    xArray(:,i) = x;
    xlinArray(:, i) = xlin;
    % Nonlinear simulation
    ua0 = sin(w*t); % nominal winding A control input
    ub0 = cos(w*t); % nominal winding B control input
    ua = (1 + ControlChange) * sin(w*t); % true winding A control input
    ub = (1 + ControlChange) * cos(w*t); % true winding B control input
    xdot = [-Ra/L*x(1) + x(3)*lambda/L*sin(x(4)) + ua/L;
        -Ra/L*x(2) - x(3)*lambda/L*cos(x(4)) + ub/L;
        -3/2*lambda/J*x(1)*sin(x(4)) + 3/2*lambda/J*x(2)*cos(x(4)) - B/J*x(3);
        x(3)];
    x = x + xdot * dt;
    x(4) = mod(x(4), 2*pi);
    % Linearized simulation
    w0 = -2*pi; % nominal rotor speed
    theta0 = 3*pi/4 - 2*pi*t; % nominal rotor position
    ia0 = 0.3708 * cos(2*pi*(t-1.36)); % nominal winding a current
    ib0 = -0.3708 * sin(2*pi*(t-1.36)); % nominal winding b current
    du = [ua-ua0; ub-ub0];
    F = [-Ra/L 0 lambda/L*sin(theta0) w0*lambda/L*cos(theta0);
        0 -Ra/L -lambda/L*cos(theta0) w0*lambda/L*sin(theta0);
        -3/2*lambda/J*sin(theta0) 3/2*lambda/J*cos(theta0) -B/J -3/2*lambda/J*(ia0*cos(theta0)+ib0*sin(theta0));
        0 0 1 0];
    G = [1/L 0; 0 1/L; 0 0; 0 0];
    dxdot = F * dx + G * du;
    dx = dx + dxdot * dt;
    xlin = [ia0; ib0; w0; theta0] + dx;
    xlin(4) = mod(xlin(4), 2*pi);
end

% % Compare linear and nonlinear simulation.
close all;
tArray = 0 : dt : tf+dt/10;
figure;
set(gcf,'Color','White'); 

subplot(2,2,1); hold on;
plot(tArray,xArray(3,:),'b-','LineWidth',1.5);
plot(tArray,xlinArray(3,:),'r--','LineWidth',1.5); 
set(gca,'FontSize',12); set(gca,'Box','on');
ylabel('Speed (Rad / Sec)');
legend('Nonlinear', 'Linearized');

subplot(2,2,2); hold on;
plot(tArray,xArray(4,:),'b-','LineWidth',1.5);
plot(tArray,xlinArray(4,:),'r--','LineWidth',1.5); 
set(gca,'FontSize',12); set(gca,'Box','on');
ylabel('Position (Radians)');

subplot(2,2,3); hold on;
plot(tArray,xArray(1,:),'b-','LineWidth',1.5);
plot(tArray,xlinArray(1,:),'r--','LineWidth',1.5); 
set(gca,'FontSize',12); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Current A (Amps)');

subplot(2,2,4); hold on;
plot(tArray,xArray(2,:),'b-','LineWidth',1.5);
plot(tArray,xlinArray(2,:),'r--','LineWidth',1.5); 
set(gca,'FontSize',12); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Current B (Amps)');
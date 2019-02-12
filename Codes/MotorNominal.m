function MotorNominal

% Two-phase step motor simulation.

Ra = 1.9; % Winding resistance
L = 0.003; % Winding inductance
lambda = 0.1; % Motor constant
J = 0.00018; % Moment of inertia
B = 0.001; % Coefficient of viscous friction

dt = 0.0005; % Integration step size

tf = 4; % Simulation length 

x = [0; 0; 0; 0]; % Initial state
w = 2 * pi; % Control input frequency

% Initialize arrays for plotting at the end of the program
N = round(tf/dt) + 1;
xArray = zeros(4, N);

i = 0;
% Begin simulation loop
for t = 0 : dt : tf
    i = i + 1;
    xArray(:, i) = x;
    ua = sin(w*t); % nominal winding A control input
    ub = cos(w*t); % nominal winding B control input
    xdot = [-Ra/L*x(1) + x(3)*lambda/L*sin(x(4)) + ua/L;
        -Ra/L*x(2) - x(3)*lambda/L*cos(x(4)) + ub/L;
        -3/2*lambda/J*x(1)*sin(x(4)) + 3/2*lambda/J*x(2)*cos(x(4)) - B/J*x(3);
        x(3)];
    x = x + xdot * dt;
    x(4) = mod(x(4), 2*pi);
end

% Plot results
close all;
tArray = 0 : dt : tf+dt/10;

figure;
set(gcf,'Color','White'); 

subplot(2,2,1); 
plot(tArray,xArray(3,:),'b-','LineWidth',1.5);
set(gca,'FontSize',12); set(gca,'Box','on');
ylabel('Speed (Rad/Sec)');

subplot(2,2,2);
plot(tArray,xArray(4,:),'b-','LineWidth',1.5);
set(gca,'FontSize',12); set(gca,'Box','on');
ylabel('Position (Radians)');

subplot(2,2,3);
plot(tArray,xArray(1,:),'b-','LineWidth',1.5);
set(gca,'FontSize',12); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Current A (Amps)');

subplot(2,2,4);
plot(tArray,xArray(2,:),'b-','LineWidth',1.5);
set(gca,'FontSize',12); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Current B (Amps)');
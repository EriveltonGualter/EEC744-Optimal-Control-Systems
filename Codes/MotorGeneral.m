function MotorGeneral(ODEOption)

% Two-phase step motor simulation written in a general way so the programmer
% can implement any kind of integration desired.
% INPUT: ODEOption = 1 (Rectangular), 2 (Trapezoidal), 3 (Runge Kutta), 4 (ODE45)

if ~exist('ODEOption', 'var')
    ODEOption = 1;
end

dt = 0.0005; % Integration step size
tf = 4; % Simulation length 

x = [0; 0; 0; 0]; % Initial state

if ODEOption == 4
    % Matlab ODE solver
    [tArray, xArray] = ode45(@CalcXdot, [0, tf], x);
    xArray = xArray';
    xArray(4,:) = mod(xArray(4,:), 2*pi);
else
    i = 0;
    N = round(tf/dt) + 1;
    % Initialize arrays for plotting at the end of the program
    tArray = zeros(1, N);
    xArray = zeros(4, N);
    % Begin simulation loop
    for t = 0 : dt : tf
        i = i + 1;
        % Save data for plotting
        xArray(:,i) = x;
        tArray(i) = t;
        % Update x
        if ODEOption == 1
            % Rectangular integration
            xdot = CalcXdot(t, x);
            x = x + xdot * dt;
        elseif ODEOption == 2
            % Trapezoidal integration
            xdot1 = CalcXdot(t, x);
            xdot2 = CalcXdot(t+dt, x+xdot1*dt);
            x = x + (xdot1 + xdot2) * dt / 2;
        elseif ODEOption == 3
            xdot1 = CalcXdot(t, x);
            xdot2 = CalcXdot(t+dt/2, x+xdot1*dt/2);
            xdot3 = CalcXdot(t+dt/2, x+xdot2*dt/2);
            xdot4 = CalcXdot(t+dt, x+xdot3*dt);
            x = x + (xdot1/6 + xdot2/3 + xdot3/3 + xdot4/6) * dt;
        end
        x(4) = mod(x(4), 2*pi);
    end
end


% Compare linear and nonlinear simulation.
close all;
figure;
set(gcf,'Color','White'); 

subplot(2,2,1); hold on;
plot(tArray,xArray(3,:),'b-','LineWidth',1.5);
set(gca,'FontSize',12); set(gca,'Box','on');
ylabel('Speed (Rad / Sec)');

subplot(2,2,2); hold on;
plot(tArray,xArray(4,:),'b-','LineWidth',1.5);
set(gca,'FontSize',12); set(gca,'Box','on');
ylabel('Position (Radians)');

subplot(2,2,3); hold on;
plot(tArray,xArray(1,:),'b-','LineWidth',1.5);
set(gca,'FontSize',12); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Current A (Amps)');

subplot(2,2,4); hold on;
plot(tArray,xArray(2,:),'b-','LineWidth',1.5);
set(gca,'FontSize',12); set(gca,'Box','on');
xlabel('Seconds'); ylabel('Current B (Amps)');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = CalcXdot(t, x)
w = 2 * pi; % Control input frequency
Ra = 1.9; % Winding resistance
L = 0.003; % Winding inductance
lambda = 0.1; % Motor constant
J = 0.00018; % Moment of inertia
B = 0.001; % Coefficient of viscous friction
ua = sin(w*t); % nominal winding A control input
ub = cos(w*t); % nominal winding B control input
xdot = [-Ra/L*x(1) + x(3)*lambda/L*sin(x(4)) + ua/L;
    -Ra/L*x(2) - x(3)*lambda/L*cos(x(4)) + ub/L;
    -3/2*lambda/J*x(1)*sin(x(4)) + 3/2*lambda/J*x(2)*cos(x(4)) - B/J*x(3);
    x(3)];
return
function tdp1(theta0, tf)

% Thrust direction programming problem.
% Inputs: theta0 = initial thrust angle (degrees)
%         tf = final time (seconds)

theta0 = theta0 * pi / 180;
x = 0; xdot = 0;
y = 0; ydot = 0;
xArr = [x]; xdotArr = [xdot];
yArr = [y]; ydotArr = [ydot];
theta = theta0; thetaArr = [theta];
dt = 0.01;
for t = dt : dt : tf
   xdot = tf * log((tan(theta0) + sec(theta0)) / (tan(theta) + sec(theta))) / 2 / tan(theta0);
   x = x + xdot * dt;
   ydot = tf * (sec(theta0) - sec(theta)) / 2 / tan(theta0);
   y = y + ydot * dt;
   thetadot = -2 * tan(theta0) * cos(theta) * cos(theta) / tf;
   theta = theta + thetadot * dt;
   thetaArr = [thetaArr theta];
   xArr = [xArr x];
   xdotArr = [xdotArr xdot];
   yArr = [yArr y];
   ydotArr = [ydotArr ydot];
end
close all;
t = 0 : dt : tf;
figure;
subplot(2,2,1);
plot(t, (180/pi)*thetaArr);
title('Thrust Direction Programming');
xlabel('time'); ylabel('Thrust angle (degrees)');

subplot(2,2,2);
plot(t, xdotArr);
xlabel('time'); ylabel('x velocity'); 

subplot(2,2,3);
plot(xArr, yArr);
xlabel('Horizontal Position');
ylabel('Vertical Position');

subplot(2,2,4);
plot(t, ydotArr);
xlabel('time'); ylabel('y velocity');

disp(['xdot = ',num2str(xdot)]);
disp(['yf = ',num2str(y)]);
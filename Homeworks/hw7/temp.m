% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
% 
% Erivelton Gualter, 03/19/2018

clear all; clc; close all

% Boundary Conditions
xo = 3.66;  
yo = -1.86;
xf = 0;     
yf = 0;

ts = 0.01;          % Sampling time
q = 0:0.001:2*pi;   % angle guess array
t = linspace(0,10,length(q)); % Array time

Z = zeros(length(t),3); % Initialize states
minimum = Inf;          % Variable to find the best gues

% Brush the angle array to find the best angle guess
for j=1:length(q)
    Z(1,:) = [ xo yo q(j)]; % Initial state condition
    
    % Euler Integration
    for i=1:length(t)-1 

        % Load states
        z1 = Z(i,1);
        z2 = Z(i,2);
        z3 = Z(i,3);

        % Dot States
        z1d = cos(z3)-z2;
        z2d = sin(z3);
        z3d = cos(z3)^2;

        Zdot = [z1d z2d z3d];

        % Integration
        Z(i+1,:) = Z(i,:) + Zdot*ts; 

        % Euclidian distance of the final position to [xf yf]
        best = norm(Z(i+1,1:2)-[xf yf]);
        if best < minimum
            minimum = best; % Update minimum flag
            
            Zopt = Z;           % Final States
            topt = t(1:i);      % time array
            Zopt = Zopt(1:i,:); 
            tfinal = t(i); 
        end
    end
end

figure; hold on;
plot(Zopt(:,1), Zopt(:,2),'LineWidth',2);
plot(xo, yo, '*','LineWidth',2);
plot(xf, yf, '*','LineWidth',2);
axis equal

% Find the angles using the guess condition
syms theta thetaf
sol = vpasolve([yo == sec(theta) - sec(thetaf), ...
                xo == -(sec(thetaf)*(tan(thetaf)-tan(theta)) - ...
                tan(theta)*(sec(thetaf)-sec(theta)) + ...
                log((tan(thetaf)+sec(thetaf))/(tan(theta)+sec(theta))))/2] ...
                ,[theta thetaf],[Zopt(1,3) Zopt(end,3)]); %[2.6270 1.2935]);
            
Z(1,:) = [ xo yo double(sol.theta)]; % Initial condition
% Euler Integration Method
for i=1:length(t)-1

    % Load states
    z1 = Z(i,1);
    z2 = Z(i,2);
    z3 = Z(i,3);

    % Dot States
    z1d = cos(z3)-z2;
    z2d = sin(z3);
    z3d = cos(z3)^2;

    Zdot = [z1d z2d z3d];

    % Integration
    Z(i+1,:) = Z(i,:) + Zdot*ts; 
end
    
figure; hold on;
plot(Zopt(:,1), Zopt(:,2),'LineWidth',2);
plot(xo, yo, '*','LineWidth',2);
plot(xf, yf, '*','LineWidth',2);
axis equal

y = Zopt(:,2);
x = Zopt(:,1);

idx = 1:10:length(x);
x=x(idx);
y=y(idx);

dy = gradient(y);
dx = gradient(x);

figure; hold on;
quiver(x,y,dy,dx);
plot(Zopt(:,1), Zopt(:,2),'LineWidth',2);
axis equal

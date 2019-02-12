% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 02/21/2018

clear all; clc; close all

% parameter
theta0 = pi/2;
X0 = 0; 
Y0 = 0;
    
Xf = [1, 1, 3]; 
Yf = [1, 3, 1];

[X(1,:), Y(1,:)] = brachistochrone(X0, Y0, Xf(1), Yf(1), theta0);
[X(2,:), Y(2,:)] = brachistochrone(X0, Y0, Xf(2), Yf(2), theta0);
[X(3,:), Y(3,:)] = brachistochrone(X0, Y0, Xf(3), Yf(3), theta0);

Xf(4) = 0.380000000000000; 
Yf(4) = 0.623586752778062;
[X(4,:), Y(4,:)] = brachistochrone(X0, Y0, Xf(4), Yf(4), theta0);
    
f1 = figure; plot(X(1,:),Y(1,:), 'LineWidth', 2); 
set(gca, 'Ydir','reverse');
title_str = strcat('Brachistochrone Solution for $x_f: $',num2str(Xf(1)),' and $y_f: $',num2str(Yf(1)));
title(title_str, 'Interpreter','Latex', 'FontSize',14);
axis equal; axis([min(X(1,:)) max(X(1,:)) min(Y(1,:)) max(Y(1,:))]);
saveFigureToPdf('f1',f1);

f2 = figure; plot(X(2,:),Y(2,:), 'LineWidth', 2); 
set(gca, 'Ydir','reverse');
title_str = strcat('Brachistochrone Solution for $x_f: $',num2str(Xf(1)),' and $y_f: $',num2str(Yf(1)));
title(title_str, 'Interpreter','Latex', 'FontSize',14);
axis equal; axis([min(X(2,:)) max(X(2,:)) min(Y(2,:)) max(Y(2,:))]);
saveFigureToPdf('f2',f2);

f3 = figure; plot(X(3,:),Y(3,:), 'LineWidth', 2); 
set(gca, 'Ydir','reverse');
title_str = strcat('Brachistochrone Solution for $x_f: $',num2str(Xf(1)),' and $y_f: $',num2str(Yf(1)));
title(title_str, 'Interpreter','Latex', 'FontSize',14);
axis equal; axis([min(X(3,:)) max(X(3,:)) min(Y(3,:)) max(Y(3,:))]);
saveFigureToPdf('f3',f3);

f4 = figure; hold on; 
plot(X(1,:),Y(1,:), 'LineWidth', 2); plot(X(4,:),Y(4,:), '--', 'LineWidth', 3); 
set(gca, 'Ydir','reverse');
title_str = strcat('Brachistochrone Solution for $x_f: $',num2str(Xf(1)),' and $y_f: $',num2str(Yf(1)));
title(title_str, 'Interpreter','Latex', 'FontSize',14);
axis equal; axis([min(X(1,:)) max(X(1,:)) min(Y(1,:)) max(Y(1,:))]);
legend('[0,0]-[1,1]','[0,0]-[0.4,0.6]');
saveFigureToPdf('f4',f4);

function [X,Y] = brachistochrone(X0, Y0, xf, yf, theta0)

    N = 200;

    fun = @(thetaf) f1([X0; Y0], [xf; yf], theta0, thetaf);  % function of x alone
    tf = fzero(fun,0);

    X = X0:(xf-X0)/N:xf(1);
    for i=1:length(X)

        fun = @(theta) f2(X(i), [xf; yf], theta, tf);  % function of x alone
        theta_out = fzero(fun,0);

        Y(i) = yf*cos(theta_out)^2/(cos(tf)^2);

    end
    
    % Functions
    function minxy = f1(P0, Pf, theta, thetaf)

        x = P0(1);
        y = P0(2);
        xf = Pf(1);
        yf = Pf(2);

        minxy = x - (xf + yf/(2*cos(thetaf)^2) * (2*(thetaf-theta) + ...
            sin(2*thetaf)) - sin(2*theta));

        minxy = minxy + y - yf/(2*cos(thetaf)^2)*cos(theta)^2;
    end

    function minx = f2(x, Pf, theta, thetaf)

        xf = Pf(1);
        yf = Pf(2);

        minx = x - (xf + yf/(2*cos(thetaf)^2) * (2*(thetaf-theta) + sin(2*thetaf) ...
            - sin(2*theta)));
    end
end
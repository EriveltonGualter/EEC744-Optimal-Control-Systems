% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 02/28/2018

clear all; clc; close all

%% Question 3-a
a = 2;

myfun = @(L,h,a) L - 2 * h * sinh(a / 2/ h);
 
L = 2.05; fun1 = @(h) myfun(L,h,a); h1 = fzero(fun1,0.1); h1=abs(h1);
L = 2.90; fun2 = @(h) myfun(L,h,a); h2 = fzero(fun2,0.1); h2=abs(h2);

x = -a/2:0.01:a/2; 

y1 = h1 * (cosh(x / h1) - cosh(a / 2 / h1));
y2 = h2 * (cosh(x / h2) - cosh(a / 2 / h2));

f1 = figure; plot(x, y1, x, y2, 'LineWidth', 2); 

title('Catenary Problem', 'Interpreter','Latex', 'FontSize',14);
xlabel('x', 'Interpreter','Latex', 'FontSize',14);
ylabel('y', 'Interpreter','Latex', 'FontSize',14);
legend('L = 2.05', 'L = 2.90');

%% Question 3 - natural catenary
f4 = figure; hold on
ax1 = subplot(121); 
ax2 = subplot(122);
%% 
a = 2; syms h;
dy = 2; % 2, 2.5, 3, 3.5, 4, 4.5, 5
hn = vpasolve(h * cosh(a / 2 / h) == dy, h ); 
Ln =  2 * hn * sinh(a / 2 / hn)
yn = hn * (cosh(x / hn));
yn2 = hn * (cosh(x / hn))-dy;

hold(ax1,'on'); plot(ax1, x, yn, 'LineWidth', 2);  hold(ax1,'off')
    title(ax1, 'Natural Catenary for several B.C.', 'Interpreter','Latex', 'FontSize',14);
    xlabel(ax1, 'x', 'Interpreter','Latex', 'FontSize',14);
    ylabel(ax1, 'y', 'Interpreter','Latex', 'FontSize',14);
    xlim(ax1,[-1.5 1.5]);
    grid(ax1,'on')
    
hold(ax2,'on'); plot(ax2, x, yn2, 'LineWidth', 2);  hold(ax2,'off')
    title(ax2,'Natural Catenary Shift', 'Interpreter','Latex', 'FontSize',14);
    xlabel(ax1, 'x', 'Interpreter','Latex', 'FontSize',14);
    ylabel(ax1, 'y', 'Interpreter','Latex', 'FontSize',14);
    xlim(ax2,[-1.5 1.5]);
    grid(ax2,'on')

saveFigureToPdf('f4',f4);

%% Question 3-b
syms h
h3 = vpasolve(-h * (cosh(0 / h) - cosh(a / 2 / h)) == 2, h );

L3 = 2 * h3 * sinh(a / 2/ h3);

y3 = h3 * (cosh(x / h3) - cosh(a / 2 / h3));

f2 = figure; plot(x, y3, 'LineWidth', 2); 

title('Catenary Problem', 'Interpreter','Latex', 'FontSize',14);
xlabel('x', 'Interpreter','Latex', 'FontSize',14);
ylabel('y', 'Interpreter','Latex', 'FontSize',14);
legend(['L = ' num2str(double(L3))]);

saveFigureToPdf('f2',f2);

%% Question 4
close all
r = 6; a = 5; 

f3 = figure;
subplot(121);
y = -4+0.6823;
circle(0,y,r); 
semicrc = r.*[cos(acos(a/r):0.01:(pi-acos(a/r))); ...
              sin(acos(a/r):0.01:(pi-acos(a/r)))];
semicrc(2,:) = semicrc(2,:) + y;

area(semicrc(1,:), semicrc(2,:), 'FaceColor','b')
axis([-6 6 -15 15]); axis equal

subplot(122);
y = 3;
circle(0,y,r); 
semicrc = r.*[cos(acos(a/r):0.01:(pi-acos(a/r))); ...
              sin(acos(a/r):0.01:(pi-acos(a/r)))];
semicrc(2,:) = semicrc(2,:) + y;

area(semicrc(1,:), semicrc(2,:), 'FaceColor','b')
axis([-6 6 -15 15]); axis equal
saveFigureToPdf('f3',f3);





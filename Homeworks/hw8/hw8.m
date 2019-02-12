% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 03/26/2018

clear all; close all;

% Plant
A = [0 1; 0 -0.02];
B = [0 ; 0.02];

% Simulation Parametes
X0 = [0.2; 0];  % Initial States
tf = 0.4; %0.2; % Final time [s]
dt = 1e-3;      % Sampling time
N = tf/dt+3;
t = 0:dt:tf;

% LQR Parameters
H = zeros(size(A));
Q = [1 0; 0 0];
R = 1e-8;

%% A, B and C probelm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ad, Bd] = c2d(A,B, dt);
N = tf/dt+2;
P(:,:,1)= eye(2);    
for k=2:N-1 % Backwards loop to find P and K
    S = R + Bd'*P(:,:,k-1)*Bd;
    F(:,:,N-k) =  -( inv(S) * Bd' * P(:,:,k-1) * Ad );
    P(:,:,k) = (Ad + Bd*F(:,:,N-k))'*P(:,:,k-1)*(Ad + Bd*F(:,:,N-k)) + F(:,:,N-k)'*R*F(:,:,N-k) + Q;
end

X1 = X0;
X2 = X0;
for k=1:N-2  % Simulate systen using Euler Integration

    % Time varying 
    u1(k) = F(:,:,k)*X1(:,k);
    
    XD1 = A*X1(:,k) + B*u1(k);
    X1(:,k+1) = X1(:,k) + XD1*dt;
    
    % Steady-stae
    Pss = dare(Ad,Bd,Q,R);
    S = R + Bd'*Pss*Bd;
    Fss = -( inv(S) * Bd' * Pss * Ad );
    u2(k) = Fss*X2(:,k);   
    
    XD2 = A*X2(:,k) + B*u2(k);
    X2(:,k+1) = X2(:,k) + XD2*dt;
end

% Performance Measure (calc cost)
Jt1 = 0;
Jt2 = 0;
for k=1:N-2
    J1(k) = X1(:,end)'*H*X1(:,end)/2 + (X1(:,k)'*Q*X1(:,k) + u1(k)*R*u1(k))/2;
    Jt1 = Jt1 + J1(k);
    
    J2(k) = X2(:,end)'*H*X2(:,end)/2 + (X2(:,k)'*Q*X2(:,k) + u2(k)*R*u2(k))/2;
    Jt2 = Jt2 + J2(k);
end

%%
close all
f1 = figure;
ax1 = subplot(211); hold on; plot(t, X1(1,:),'LineWidth',2); plot(t, X2(1,:),'--','LineWidth',2);
ax2 = subplot(212); hold on; plot(t, X1(2,:),'LineWidth',2); plot(t, X2(2,:),'--','LineWidth',2);
title(ax1, 'State $x_1$','Interpreter','latex','FontSize',14); 
title(ax2, 'State $x_2$','Interpreter','latex','FontSize',14);   
legend(ax1, 'Time-varying LQR', 'Steady-state LQR');
legend(ax2, 'Time-varying LQR', 'Steady-state LQR');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig1',f1);
    
f2 = figure; hold on;
tu = t(1:end-1);
plot(tu, u1,'LineWidth',2); plot(tu, u2,'--','LineWidth',2);
legend('Time-varying LQR', 'Steady-state LQR');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
ylabel('Optimal Control Law', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig2',f2);

f3 = figure; hold on;
F_plot1(:,:) = F(1,:,:); 
F_plot2(1,:) = Fss(1)*ones(1,length(F));
F_plot2(2,:) = Fss(2)*ones(1,length(F));
plot(t(1:end-1), F_plot1(1,:), t(1:end-1), F_plot1(2,:),'LineWidth',2);
plot(t(1:end-1), F_plot2(1,:), t(1:end-1), F_plot2(2,:),'--','LineWidth',2);
legend('Time-varying LQR', 'Steady-state LQR');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
title('State Feedback Matrix', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig3',f3);

f4 = figure; hold on;
plot(J1,'LineWidth',2); plot(J2,'--','LineWidth',2);
xlim([1 length(J1)]);
title('Cost', 'Interpreter','Latex', 'FontSize',14);
legend('Time-varying LQR', 'Steady-state LQR');

saveFigureToPdf('fig4',f4);

%% D Probelm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(t);
X3 = X0;
for k=1:N-1

    K = care(A,B,Q,R);
    u3(k) = -inv(R)*B'*K*X3(:,k);
    
    XD3 = A*X3(:,k) + B*u3(k);
    
    X3(:,k+1) = X3(:,k) + XD3*dt;    
end

% Performance Measure (calc cost)
Jt3 = 0;
for k=1:N-2
    J3(k) = X3(:,end)'*H*X3(:,end)/2 + (X3(:,k)'*Q*X3(:,k) + u3(k)*R*u3(k))/2;
    Jt3 = Jt3 + J3(k);
end
%%
f5 = figure;
ax1 = subplot(211); hold on; plot(t, X2(1,:),'LineWidth',2); plot(t, X3(1,:),'--','LineWidth',2);
ax2 = subplot(212); hold on; plot(t, X2(2,:),'LineWidth',2); plot(t, X3(2,:),'--','LineWidth',2);
title(ax1, 'State $x_1$','Interpreter','latex','FontSize',14); 
title(ax2, 'State $x_2$','Interpreter','latex','FontSize',14);   
legend(ax1, 'Steady-state LQR',  'Hamiltonian matrix approach');
legend(ax2, 'Steady-state LQR',  'Hamiltonian matrix approach');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig5',f5);
    
f6 = figure; hold on;
tu = t(1:end-1);
plot(tu, u2,'LineWidth',2); plot(tu, u3,'--','LineWidth',2);
legend('Steady-state LQR',  'Hamiltonian matrix approach');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
ylabel('Optimal Control Law', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig6',f6);

f7 = figure; hold on;
gain = -inv(R)*B'*K;
F_plot3(1,:) = gain(1)*ones(1,length(F));
F_plot3(2,:) = gain(2)*ones(1,length(F));
plot(t(1:end-1), F_plot2(1,:), t(1:end-1), F_plot2(2,:),'LineWidth',2);
plot(t(1:end-1), F_plot3(1,:), t(1:end-1), F_plot3(2,:),'--','LineWidth',2);
legend('Steady-state LQR',  'Hamiltonian matrix approach');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
title('State Feedback Matrix', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig7',f7);

f8 = figure; hold on;
plot(J2,'LineWidth',2); plot(J3,'--','LineWidth',2);
xlim([1 length(J2)]);
title('Cost', 'Interpreter','Latex', 'FontSize',14);
legend('Steady-state LQR',  'Hamiltonian matrix approach');
saveFigureToPdf('fig8',f8);
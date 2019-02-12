% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 02/08/2018
%

clear all, close all, clc

% Plant
A = [0 1; -1 2];
B = [0; 1];

% Simulation Parameters
T = 10;
dt = 0.02;
t = 0:dt:T; 
X0 = [1; 1];

% LQR Parameters
H = [1 0; 0 0];
Q = [1 0; 0 0];
R = 1;
C = 1;
r = [sin(t); t/2];
% r = zeros(2,length(t))

sys=ss(A,B,[],[]);
sysd=c2d(sys,dt,'zoh');
A=sysd.a;
B=sysd.b;

%% Riccati Approach
N = length(t);
K(:,:,N) = C'*H*C;
V(:,:,N) = C'*H*r(end);
[A, B] = c2d(A,B,dt)
for k=N:-1:2
    
    Kd(:,:,k) = -C'*Q*C + K(:,:,k)'*B*inv(R)*B'*K(:,:,k) - K(:,:,k)*A - A'*K(:,:,k);
    Vd(:,:,k) = -(A'-K(:,:,k)'*B*inv(R)*B')*V(:,:,k) - C'*Q*r(:,k);
    
    K(:,:,k-1) = K(:,:,k) - Kd(:,:,k)*dt;
    V(:,:,k-1) = V(:,:,k) - Vd(:,:,k)*dt;
end
ktemp = K(1,1,:); K11(:) = ktemp(:,:,:);
ktemp = K(1,2,:); K12(:) = ktemp(:,:,:);
ktemp = K(2,1,:); K21(:) = ktemp(:,:,:);
ktemp = K(2,2,:); K22(:) = ktemp(:,:,:);

A = [0 1; -1 2];
B = [0; 1];

X = X0;
for k=1:N-1

    u(k)  = -inv(R)*B'*K(:,:,k)*X(:,k);
    SFM(:,N-k+1) = -inv(R)*B'*K(:,:,N-k+1);
    
    XD = A*X(:,k) + B*u(k);
    
    X(:,k+1) = X(:,k) + XD*dt;    
end

% Performance Measure (calc cost)
for k=1:N-2
    J(k) = X(:,end)'*H*X(:,end)/2 + (X(:,k)'*Q*X(:,k) + u(k)*R*u(k))/2;
end

% Performance Measure (calc cost)
for k=1:N-2
    J(k) = (C*X(:,end)-r(:,end))'*H*(C*X(:,end)-r(:,end))/2 + ((C*X(:,k)-r(:,end))'*Q*(C*X(:,k)-r(:,k)) + u(k)*R*u(k))/2;
end

%%
fig1 = figure
ax1 = subplot(211); hold on; plot(t, r(1,:),'LineWidth',2); plot(t, X(1,:),'LineWidth',2);
ax2 = subplot(212); hold on; plot(t, r(2,:),'LineWidth',2); plot(t, X(2,:),'LineWidth',2);
% subplot(211); hold on; plot(t, X(1,:));
% subplot(212); hold on; plot(t, X(2,:));
% figure
% subplot(211); hold on; plot(t, X(1,:),'LineWidth',2);
% subplot(212); hold on; plot(t, X(2,:),'LineWidth',2);
title(ax1, 'State $x_1$','Interpreter','latex','FontSize',14); 
title(ax2, 'State $x_2$','Interpreter','latex','FontSize',14);   
% legend(ax1, 'Method 1', 'Method 2');
% legend(ax2, 'Method 1', 'Method 2');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('lq_lqr',fig1);


fig = figure
plot(t(2:end), SFM(:,2:end),'LineWidth',2);
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
title('State Feedback Matrix', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('stm2',fig);

figure
plot(t(1:end-1), u, 'LineWidth',2);

figure; plot(J)
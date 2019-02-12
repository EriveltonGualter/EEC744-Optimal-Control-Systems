% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 02/08/2018
% LQR

clear all, close all, clc

% Plant
A = [0 1; -1 2];
B = [0; 1];

% LQR Parameters
H = [2 0; 0 0];
Q = [10 0; 0 2];
R = 1;

% Simulation Parameters
T = 10;
dt = 0.02;
t = 0:dt:T; 
X0 = [1; 1];

%% Riccati Approach
N = length(t);
K(:,:,N) = H;
for k=N:-1:2
    Kd(:,:,k) = -Q + K(:,:,k)*B*inv(R)*B'*K(:,:,k) - K(:,:,k)*A - A'*K(:,:,k);
    K(:,:,k-1) = K(:,:,k) - Kd(:,:,k)*dt;
end
ktemp = K(1,1,:); K11(:) = ktemp(:,:,:);
ktemp = K(1,2,:); K12(:) = ktemp(:,:,:);
ktemp = K(2,1,:); K21(:) = ktemp(:,:,:);
ktemp = K(2,2,:); K22(:) = ktemp(:,:,:);

X = X0;
for k=1:N-1

    u(k)  = -inv(R)*B'*K(:,:,k)*X(:,k);
    SFM(:,N-k+1) = -inv(R)*B'*K(:,:,N-k+1);
    
    XD = A*X(:,k) + B*u(k);
    
    X(:,k+1) = X(:,k) + XD*dt;    
end

%%
[Ad, Bd] = c2d(A,B, dt);
% Ad = A; Bd = B;
N = T/dt+2;
P(:,:,1)= eye(2);    
for k=2:N-1
    S = R + Bd'*P(:,:,k-1)*Bd;
    F(:,:,N-k) =  -( inv(S) * Bd' * P(:,:,k-1) * Ad );
    P(:,:,k) = (Ad + Bd*F(:,:,N-k))'*P(:,:,k-1)*(Ad + Bd*F(:,:,N-k)) + F(:,:,N-k)'*R*F(:,:,N-k) + Q;
end

Xric = X0;
for k=1:N-2

    %uric(k) = F(:,:,N-k-1)*Xric(:,k);
    uric(k) = F(:,:,k)*Xric(:,k);
    
    XD = A*Xric(:,k) + B*uric(k);
    
    Xric(:,k+1) = Xric(:,k) + XD*dt;
end

% Performance Measure (calc cost)
for k=1:N-2
    J(k) = X(:,end)'*H*X(:,end)/2 + (X(:,k)'*Q*X(:,k) + u(k)*R*u(k))/2;
end

%%
f1 = figure;
ax1 = subplot(211); hold on; plot(t, X(1,:), t, Xric(1,:),'LineWidth',2);
ax2 = subplot(212); hold on; plot(t, X(2,:), t, Xric(2,:),'LineWidth',2);
title(ax1, 'State $x_1$','Interpreter','latex','FontSize',14); 
title(ax2, 'State $x_2$','Interpreter','latex','FontSize',14);   
legend(ax1, 'Method 1', 'Method 2');
legend(ax2, 'Method 1', 'Method 2');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig1',f1);
    
f2 = figure; hold on;
tu = t(1:end-1);
plot(tu, u, tu, uric,'LineWidth',2);
legend('Method 1', 'Method 2');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
ylabel('Optimal Control Law', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig2',f2);

f3 = figure; hold on;
F_plot(:,:) = F(1,:,:); plot(t(2:end), SFM(:,2:end), t(1:end-1), F_plot,'LineWidth',2);
legend('Method 1', 'Method 1', 'Method 2', 'Method 2');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
title('State Feedback Matrix', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig3',f3);

figure; plot(J)
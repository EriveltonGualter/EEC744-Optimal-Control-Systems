% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 01/22/2018

clear all; close all;

% Plant
A = [0 1; 0 -0.02];
B = [0 ; 0.02];

% Simulation Parametes
X0 = [0.2; 0];  % Initial States
tf = 0.4;       % Final time [s]
dt = 1e-3;      % Sampling time
N = tf/dt+3;
t = 0:dt:tf;
TV = 0; % TV = 1 (time varying) or TV = 0 (steady state)

% LQR Parameters
H = zeros(size(A));
Q = [1 0; 0 0];
R = 1e-8;

%%
[Ad, Bd] = c2d(A,B, dt);
N = tf/dt+2;
P(:,:,1)= eye(2);    
for k=2:N-1 % Backwards loop to find P and K
    S = R + Bd'*P(:,:,k-1)*Bd;
    F(:,:,N-k) =  -( inv(S) * Bd' * P(:,:,k-1) * Ad );
    P(:,:,k) = (Ad + Bd*F(:,:,N-k))'*P(:,:,k-1)*(Ad + Bd*F(:,:,N-k)) + F(:,:,N-k)'*R*F(:,:,N-k) + Q;
end

X = X0;
for k=1:N-2  % Simulate systen using Euler Integration

    if TV 
        u(k) = F(:,:,k)*X(:,k);
    else
        Pss = dare(Ad,Bd,Q,R);
        S = R + Bd'*Pss*Bd;
        Fss = -( inv(S) * Bd' * Pss * Ad );
        u(k) = Fss*X(:,k);
    end   
    
    XD = A*X(:,k) + B*u(k);
    
    X(:,k+1) = X(:,k) + XD*dt;
end

% Performance Measure (calc cost)
jt = X(:,end)'*H*X(:,end)/2;
for k=1:N-2
    J(k) = X(:,end)'*H*X(:,end)/2 + (X(:,k)'*Q*X(:,k) + u(k)*R*u(k))/2;
    jt = jt + (X(:,k)'*Q*X(:,k) + u(k)*R*u(k))/2;
end
jt
%%
N = length(t);
K(:,:,N) = H;
for k=N:-1:2
    Kd(:,:,k) = -K(:,:,k)*A - A'*K(:,:,k) - Q + K(:,:,k)*B*inv(R)*B'*K(:,:,k);
    K(:,:,k-1) = K(:,:,k) - Kd(:,:,k)*dt;
end

ktemp = K(1,1,:); K11(:) = ktemp(:,:,:);
ktemp = K(1,2,:); K12(:) = ktemp(:,:,:);
ktemp = K(2,1,:); K21(:) = ktemp(:,:,:);
ktemp = K(2,2,:); K22(:) = ktemp(:,:,:);

X2 = X0;
for k=1:N-1

    if TV
        u(k)  = -inv(R)*B'*K(:,:,k)*X2(:,k);
        SFM(:,N-k+1) = -inv(R)*B'*K(:,:,N-k+1);
    else
        K = care(A,B,Q,R);
        u(k) = -inv(R)*B'*K*X2(:,k);
    end
    
    XD = A*X2(:,k) + B*u(k);
    
    X2(:,k+1) = X2(:,k) + XD*dt;    
end

%%
close all
f1 = figure;
ax1 = subplot(211); hold on; plot(t, X(1,:),'LineWidth',2); plot(t, X2(1,:),'LineWidth',2);
ax2 = subplot(212); hold on; plot(t, X(2,:),'LineWidth',2); plot(t, X2(2,:),'LineWidth',2);
title(ax1, 'State $x_1$','Interpreter','latex','FontSize',14); 
title(ax2, 'State $x_2$','Interpreter','latex','FontSize',14);   
% legend(ax1, 'Method 1', 'Method 2');
% legend(ax2, 'Method 1', 'Method 2');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig1',f1);
    
f2 = figure; hold on;
tu = t(1:end-1);
plot(tu, u,'LineWidth',2);
% legend('Method 1', 'Method 2');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
ylabel('Optimal Control Law', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig2',f2);

f3 = figure; hold on;
if TV
    F_plot(:,:) = F(1,:,:); 
else
    F_plot = Fss*ones(length(Fss), length(F));
    F_plot(1,:) = Fss(1)*F_plot(1,:);
    F_plot(2,:) = Fss(2)*F_plot(1,:);
end
plot(t(1:end-1), F_plot,'LineWidth',2);
% legend('Method 1', 'Method 1', 'Method 2', 'Method 2');
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
title('State Feedback Matrix', 'Interpreter','Latex', 'FontSize',14);
saveFigureToPdf('fig3',f3);

f4 = figure; plot(J,'LineWidth',2); xlim([1 length(J)]);
title('Cost', 'Interpreter','Latex', 'FontSize',14);

saveFigureToPdf('fig4',f4);
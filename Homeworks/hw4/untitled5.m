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
H = [1 0; 0 1];
Q = [1 0; 0 1];
R = 1;
C = [1 0; 0 1];
r = [sin(t); t/2];
% r = zeros(2,length(t));
% r = ones(2,length(t));

%% Riccati Approach
N = length(t);
S(:,:,N) = C'*H*C;
V(:,:,N) = C'*H*r(end);
K(:,:,N) = inv(R)*B'*S(:,:,N);
for k=N:-1:2
    
    Sd(:,:,k) = -A'*S(:,:,k) - S(:,:,k)*A + S(:,:,k)*B*inv(R)*B'*S(:,:,k) - C'*Q*C;   
    Vd(:,:,k) = - (A-B*K(:,:,k))'*V(:,:,k) - C'*Q*r(:,k);
    
    S(:,:,k-1) = S(:,:,k) - Sd(:,:,k)*dt;
    K(:,:,k) = inv(R)*B'*S(:,:,k-1);
    V(:,:,k-1) = V(:,:,k) - Vd(:,:,k)*dt;
end
% ktemp = K(1,1,:); K11(:) = ktemp(:,:,:);
% ktemp = K(1,2,:); K12(:) = ktemp(:,:,:);
% ktemp = K(2,1,:); K21(:) = ktemp(:,:,:);
% ktemp = K(2,2,:); K22(:) = ktemp(:,:,:);

% A = [0 1; -1 2];
% B = [0; 1];


X = X0;
for k=1:N-1

    u(k)  = -inv(R)*B'*( K(:,:,N-k)*X(:,k)-r(:,k) );        
    XD = A*X(:,k) + B*u(k);
    
    X(:,k+1) = X(:,k) + XD*dt;    
end

% Performance Measure (calc cost)
for k=1:N-2
    J(k) = (C*X(:,end)-r(:,end))'*H*(C*X(:,end)-r(:,end))/2 + ((C*X(:,k)-r(:,end))'*Q*(C*X(:,k)-r(:,k)) + u(k)*R*u(k))/2;
end

%%
subplot(211); hold on; plot(r(1,:)); plot(X(1,:));
subplot(212); hold on; plot(r(2,:)); plot(X(2,:));
% subplot(211); hold on; plot(t, X(1,:));
% subplot(212); hold on; plot(t, X(2,:));
figure
subplot(211); hold on; plot(X(1,:));
subplot(212); hold on; plot(X(2,:));

figure; plot(J)
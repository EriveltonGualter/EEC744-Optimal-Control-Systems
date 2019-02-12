% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 02/08/2018
%

clear all
%close all, clc

% Plant
A = [0 1; -1 -2];
B = [0; 1];

% LQR Parameters
H = [20 0; 0 0];
Q = [10 0; 0 2];
R = 0.01;

% Simulation Parameters
T = 10;
dt = 0.02e-2;
t = 0:dt:T; 
X0 = [1; 1];

%% Riccati Approach
N = length(t);
K(:,:,1) = H;
for h=1:N
    Kd(:,:,N-h+1) = -Q + K(:,:,h)*B*inv(R)*B'*K(:,:,h) - K(:,:,h)*A - A'*K(:,:,h);
    K(:,:,h+1) = K(:,:,h) - Kd(:,:,N-h+1)*dt;
end
ktemp = K(1,1,:); K11(:) = ktemp(:,:,:);
ktemp = K(1,2,:); K12(:) = ktemp(:,:,:);
ktemp = K(2,1,:); K21(:) = ktemp(:,:,:);
ktemp = K(2,2,:); K22(:) = ktemp(:,:,:);

X = X0;
for k=1:N

    u = -inv(R)*B'*K(:,:,N-k+1)*X(:,k);
    SFM(:,k) = -inv(R)*B'*K(:,:,N-k+1);
    
    XD = A*X(:,k) + B*u;
    
    X(:,k+1) = X(:,k) + XD*dt;
end

%% Plot K for parts (a) and (b) and compare the solutions. Plot the states with x(0) = [1, 1]^T for parts (a) and (b) and compare the solutions.
% figure; plot(t,K11(1:end-1), t,K12(1:end-1), t,K21(1:end-1), t,K22(1:end-1), 'LineWidth', 2);
%         title('PLot Ks')
% figure; plot(t,X(1,1:end-1), t,X(2,1:end-1), 'LineWidth', 2);
%         title('states')
% figure; plot(t, SFM(1,:), t, SFM(2,:))

%%
P(:,:,1)= H;
[A,B] = c2d(A,B,dt);
% Backwards loop to find P and K
for k=2:N-1
    S = R + B'*P(:,:,k-1)*B;
    F(:,:,N-k) =  -( inv(S) * B' * P(:,:,k-1) * A );
    P(:,:,k) = (A + B*F(:,:,N-k))'*P(:,:,k-1)*(A + B*F(:,:,N-k)) + F(:,:,N-k)'*R*F(:,:,N-k) + Q;
end

% Simulate systen
x(1,:) = X0;
for k=1:N-2
    if k ~= 1
        x(k,:) = x_new;
    end

    u(k,:) = F(:,:,k)*x(k,:).';

    xdot(k,1) = x(k,2);
    xdot(k,2) = -2*x(k,1) -2*x(k,2) + u(k); 

    x_new = x(k,:) + xdot(k,:)*dt;
end
% figure; 
hold on;
plot(t(1:end-2), x(:,1), t(1:end-2), x(:,2)); title('states 2');
plot(t,X(1,1:end-1), t,X(2,1:end-1), 'LineWidth', 2);
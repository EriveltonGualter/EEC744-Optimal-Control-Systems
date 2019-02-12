function cost = opt( W )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    a = W(1);
    b = 0;
    c = 0;
    d = W(2);
    e = W(3);
    f = 0;
    g = 0;
    h = W(4);
    i = W(5);
    j = W(6);
    k = 0;
    l = 0;
    m = W(7);
    
    
    
    % Plant
    A = [0 1; -1 2];
    B = [0; 1];

    % Simulation Parameters
    T = 10;
    dt = 0.02;
    t = 0:dt:T; 
    X0 = [1; 1];

    % LQR Parameters
    H = [a b; c d];
    Q = [e f; g h];
    R = i;
    C = [j k;l m];
    r = [sin(t); t/2];
    % r = zeros(2,length(t));
    % r = ones(2,length(t));

    %% Riccati Approach
    N = length(t);
    K(:,:,N) = C'*H*C;
    V(:,:,N) = C'*H*r(end);
    for k=N:-1:2

        Kd(:,:,k) = -C'*Q*C + K(:,:,k)'*B*inv(R)*B'*K(:,:,k) - K(:,:,k)*A - A'*K(:,:,k);
        Vd(:,:,k) = -(A'-K(:,:,k)*B*inv(R)*B')*V(:,:,k) - C'*Q*r(:,k);

        K(:,:,k-1) = K(:,:,k) - Kd(:,:,k)*dt;
        V(:,:,k-1) = V(:,:,k) - Vd(:,:,k)*dt;
    end
    ktemp = K(1,1,:); K11(:) = ktemp(:,:,:);
    ktemp = K(1,2,:); K12(:) = ktemp(:,:,:);
    ktemp = K(2,1,:); K21(:) = ktemp(:,:,:);
    ktemp = K(2,2,:); K22(:) = ktemp(:,:,:);

    X = X0;
    for k=1:N-1

        u(k)  = -inv(R)*B'*( K(:,:,k)*X(:,k)-r(:,k) );        
        XD = A*X(:,k) + B*u(k);

        X(:,k+1) = X(:,k) + XD*dt;    
    end

%     % Performance Measure (calc cost)
%     for k=1:N-2
%         J(k) = (C*X(:,end)-r(:,end))'*H*(C*X(:,end)-r(:,end))/2 + ((C*X(:,k)-r(:,end))'*Q*(C*X(:,k)-r(:,k)) + u(k)*R*u(k))/2;
%     end
    cost = norm(rms(r'-X'))
end


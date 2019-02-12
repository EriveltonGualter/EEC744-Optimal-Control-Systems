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
dt = 1e-2;      % Sampling time
N = tf/dt+3;
t = 0:dt:tf;

%% Chapter 5 approach

% F = -inv(R)*B'*K

%% HJB approach

% Initialized parameters
P(:,:,1)= eye(size(A));
Q       = [1 0; 0 0];
R_array = 1e-8;[1e-8 1e-8];

% Plots creation
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;

figure(f1);
ax3 = subplot(211);
ax4 = subplot(212);

% For for differents R. | e.g. R_array = [0.1 1 10];
for i=1 : length(R_array)
    R = R_array(i);     % set R
    
    [Ad, Bd] = c2d(A,B,dt)
    % Backwards loop to find P and K
    for k=2:N-1
        S = R + Bd'*P(:,:,k-1)*Bd;
        F(:,:,N-k) =  -( inv(S) * Bd' * P(:,:,k-1) * Ad );
        P(:,:,k) = (Ad + Bd*F(:,:,N-k))'*P(:,:,k-1)*(Ad + Bd*F(:,:,N-k)) + F(:,:,N-k)'*R*F(:,:,N-k) + Q;
    end

    % Simulate systen
    x(1,:) = X0;
    for k=1:N-2
        if k ~= 1
            x(k,:) = x_new;
        end
        
        if i==length(R_array)
            Pss = dare(A,B,Q,R);
            S = R + B'*Pss*B;
            Fss = -( inv(S) * B' * Pss * A );
            
            u(k,:) = Fss*x(k,:).';
            
            F = Fss.*ones(size(F));
        else
            u(k,:) = F(:,:,k)*x(k,:).';
        end

        xdot(k,1) = x(k,2);
        xdot(k,2) = -0.02*x(k,2) + 0.02*u(k); 

        x_new = x(k,:) + xdot(k,:)*dt;
    end
    
    % Performance Measure (calc cost)
    for k=1:N-2
        J(k) = (x(end,:)*P(:,:,1)*x(end,:)' + x(k,:)*Q*x(k,:)' + u(k,:)*R*u(k,:)')/2;
    end

    % Find Magnitude of the largest closed-loop eigenvalue
    for ii = length(x(:,1))
        cl_sys = A+B*F(:,:,ii);
        max_eigs(ii,:) = max(abs(eig(cl_sys)));
    end
    max_eigen = max(max_eigs)
               
    %% Plots
    figure(f1)
    hold(ax3,'on'); plot(ax3, t, x(:,1),'LineWidth',2); hold(ax3,'off')
    hold(ax4,'on'); plot(ax4, t, x(:,2),'LineWidth',2); hold(ax4,'off')
    title(ax3, 'State $x_1$','Interpreter','latex','FontSize',14); 
    title(ax4, 'State $x_2$','Interpreter','latex','FontSize',14);   
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    ylabel(ax3,'$ \frac{m}{s}$', 'Interpreter','Latex', 'FontSize',14);
    ylabel(ax4,'$ \frac{m}{s^2}$', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig1',f1);

    figure(f2); hold on; 
    F_plot(:,:) = F(1,:,:); plot(t,F_plot,'LineWidth',2);
    title('State Feedback Matrix $F$','Interpreter','latex','FontSize',14); 
    xlabel('N', 'Interpreter','Latex', 'FontSize',14);
    ylabel('F', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig2',f2);

    figure(f3); hold on
    F_plot(:,:) = F(1,:,:); plot(t,F_plot,'LineWidth',2);
    xlim([max(t)-1 max(t)]);
    title('State Feedback Matrix $F$','Interpreter','latex','FontSize',14); 
    xlabel('N', 'Interpreter','Latex', 'FontSize',14);
    ylabel('F', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig3',f3);

    figure(f4); hold on
    plot(t,u,'LineWidth',2); xlim([t(1) t(end)]);
    title('Control Input','Interpreter','latex','FontSize',14); 
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    ylabel('u', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig4',f4);
    
    figure(f5); hold on
    plot(linspace(1,length(t),length(t)), J, 'LineWidth', 2)
    title('Performance Measure LQR (Cost Function)','Interpreter','latex','FontSize',14); 
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    ylabel('Cost Function', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig5',f5);
   
    clear P F
    P(:,:,1)= eye(size(A));
end


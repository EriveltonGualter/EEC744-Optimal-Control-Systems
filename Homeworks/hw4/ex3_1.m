% Book: Optimal Control Theory: An introduxtion by Donald E. Kirk
%
% Erivelton Gualter, 01/22/2018

clear all; close all;

% Mechanical System Description
M = 1; % kg
K = 2; % N/m
B = 2; % N/m/s

A = [ 0 1; -K/M -B/M];
B = [0 ; 1/M];

% Simulation Parametes
X0 = [0.2; 0];  % Initial States
tf = 30;        % Final time [s]
dt = 1e-2;      % Sampling time
N = tf/dt+3;
t = 0:dt:tf;

% Initialized parameters
P(:,:,1)= eye(size(A));
Q       = eye(size(A));
R_array = [1 1];

% Plots creation
%f0 = figure;
f1 = figure;
% f2 = figure;
% f3 = figure;
% f4 = figure;
% f5 = figure;
% f6 = figure;

figure(f1);
ax1 = subplot(211);
ax2 = subplot(212);

% Analytical solutions
syms s;
xt = ilaplace(inv(s*eye(size(A))-A)*X0);
xta =  [(exp(-t).*(cos(t) + sin(t)))/5;
         -(2*exp(-t).*sin(t))/5];

% Simulate systen
x(:,1) = X0;
[A,B] = c2d(A,B,dt);
for k=1:N-2
    u(k,:) = 0;
    XD = A*x(:,k) + B*u(k,:);
    
    x(:,k+1) = x(:,k) + XD*dt;
end

    figure(f1)
    hold(ax1,'on'); plot(ax1, t, x(1,1:end-1), t, xta(1,:), 'LineWidth',2); hold(ax1,'off')
    hold(ax2,'on'); plot(ax2, t, x(2,1:end-1), t, xta(2,:), 'LineWidth',2); hold(ax2,'off')
    title(ax1, 'State Transition Matrix to evaluate discretization','Interpreter','latex','FontSize',14); 
    %legend(ax1, ['Continuos-Time'], ['Discrete-Time']);
    %legend(ax2, ['Continuos-Time'], ['Discrete-Time']);
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    ylabel(ax1,'State $x_1$ $ \frac{m}{s}$', 'Interpreter','Latex', 'FontSize',14);
    ylabel(ax2,'State $x_2$ $ \frac{m}{s^2}$', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig0',f0);
    
%%
figure(f1);
ax3 = subplot(211);
ax4 = subplot(212);

% For for differents R. | e.g. R_array = [0.1 1 10];
for i=1 : length(R_array)
    R = R_array(i);     % set R
    
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
        xdot(k,2) = -2*x(k,1) -2*x(k,2) + u(k); 

        x_new = x(k,:) + xdot(k,:)*dt;
    end
    
    % Performance Measure (calc cost)
    for k=1:N-2
        J(k) = (x(end,:)*P(:,:,1)*x(end,:)' + x(k,:)*Q*x(k,:)' + u(k,:)*R*u(k,:)')/2;
    end

    % Find MMagnitude of the largest closed-loop eigenvalue
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
    %legend(ax3, ['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(2))], ['R = ' num2str(R_array(3))], ['R = ' num2str(R_array(end)), ' SS']);
    %legend(ax4, ['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(2))], ['R = ' num2str(R_array(3))], ['R = ' num2str(R_array(end)), ' SS']);
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    ylabel(ax3,'$ \frac{m}{s}$', 'Interpreter','Latex', 'FontSize',14);
    ylabel(ax4,'$ \frac{m}{s^2}$', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig1',f1);

    figure(f2); hold on; 
    F_plot(:,:) = F(1,:,:); plot(t,F_plot,'LineWidth',2);
    %legend(['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(2))], ['R = ' num2str(R_array(2))],['R = ' num2str(R_array(3))],['R = ' num2str(R_array(3))]);
    title('State Feedback Matrix $F$','Interpreter','latex','FontSize',14); 
    xlabel('N', 'Interpreter','Latex', 'FontSize',14);
    ylabel('F', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig2',f2);

    figure(f3); hold on
    F_plot(:,:) = F(1,:,:); plot(t,F_plot,'LineWidth',2);
    xlim([max(t)-1 max(t)]);
    %legend(['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(2))], ['R = ' num2str(R_array(2))],['R = ' num2str(R_array(3))],['R = ' num2str(R_array(3))]);
    title('State Feedback Matrix $F$','Interpreter','latex','FontSize',14); 
    xlabel('N', 'Interpreter','Latex', 'FontSize',14);
    ylabel('F', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig3',f3);

    figure(f4); hold on
    plot(t,u,'LineWidth',2);
    %legend(['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(2))], ['R = ' num2str(R_array(2))],['R = ' num2str(R_array(3))],['R = ' num2str(R_array(3))]);
    title('Control Input','Interpreter','latex','FontSize',14); 
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    ylabel('u', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig4',f4);
    
    figure(f5); hold on
    plot(t, J, 'LineWidth', 2)
    %legend(['R = ' num2str(R_array(1))], ['R = ' num2str(R_array(2))], ['R = ' num2str(R_array(3))], ['R = ' num2str(R_array(end)), ' SS']);
    title('Performance Measure LQR (Cost Function)','Interpreter','latex','FontSize',14); 
    xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
    ylabel('Cost Function', 'Interpreter','Latex', 'FontSize',14);
    saveFigureToPdf('fig5',f5);
    
    clear P F
    P(:,:,1)= eye(size(A));
end

%%
% figure(f5);
% clear kf_array kp_array F_array P_array
% 
% kp_array = 1:1:N-1;   P_array(:,:) = P(1,:,:); P_array = fliplr(P_array);
% kf_array = 1:1:N-2; F_array(:,:) = F(1,:,:);
% ku_array = 0:1:N-3; 
% 
% subplot(311); plot(kp_array, P_array, 'o'); grid on; xlim([0 N]);
% subplot(312); plot(kf_array, F_array, 'o'); grid on; xlim([0 N]);
% subplot(313); plot(ku_array, u, 'o'); grid on; xlim([0 N]);



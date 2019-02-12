
%%%%%                                                                 %%%%%
%%%%%                         OPTIMAL CONTROL                         %%%%%
%%%%%                                                                 %%%%%
%%%%%                                                                 %%%%%
%%%%%                           Homework  4                           %%%%%


close all
clear
clc


%%%%% Variables:
ts=0.02;                        % Units: s.
t=0:ts:10;                      % Units: s.
N=max(t)/ts;                    % Units: -.


%%%%% State Space:
A=[0 1;-1 2];
B=[0;1];


%%%%% Cost function parameters:
syms x1 x2 u; x=[x1;x2];
H=[20 0;0 0];R=1;Q=[1 0;0 2];
xTHx=0.5*x.'*H*x;
xTQx=x.'*Q*x;
uTRu=u.'*R*u;
K(:,:,N)=H;
K_dot(:,:,N)=-Q+K(:,:,N)*B*inv(R)*B'*K(:,:,N)-K(:,:,N)*A-A'*K(:,:,N);


%%%%% K solving:
for k=1:1:N-1
    K(:,:,N-k)=K(:,:,N-k+1)-K_dot(:,:,N-k+1)*ts;
    K_dot(:,:,N-k)=-Q+K(:,:,N-k)*B*inv(R)*B'*K(:,:,N-k)-...
    K(:,:,N-k)*A-A'*K(:,:,N-k);
end


%%%%% Euler discretization:
x_K(:,1)=[1;1];
for k=1:1:N-1
    u_K(k)=-inv(R)*B'*K(:,:,k)*x_K(:,k);
    x_K(:,k+1)=([0 ts;-1*ts 2*ts]+eye(2))*x_K(:,k)+...
    [0;ts]*u_K(k);
end


%%%%% LQR:
sys=ss(A,B,[],[]);
sysd=c2d(sys,ts,'zoh');
Ad=sysd.a;
Bd=sysd.b;
P0=H;
F(N-1,:)=-inv(R+Bd'*P0*Bd)*Bd'*P0*Ad;
P(:,:,1)=(Ad+Bd*F(N-1,:))'*P0*(Ad+Bd*F(N-1,:))+F(N-1,:)'*R*F(N-1,:)+Q;
for k=2:1:N-1
   F(N-k,:)=-inv(R+Bd'*P(:,:,k-1)*Bd)*Bd'*P(:,:,k-1)*Ad;
   P(:,:,k)=(Ad+Bd*F(N-k,:))'*P(:,:,k-1)*(Ad+Bd*F(N-k,:))+F(N-k,:)'*R*F(N-k,:)+Q;
end


%%%%% Euler discretization:
x_F(:,1)=[1;1];
for k=1:1:N-1
    u_F(k)=F(k,:)*x_F(:,k);
    x_F(:,k+1)=([0 ts;-1*ts 2*ts]+eye(2))*x_F(:,k)+[0;ts]*u_F(k);
end


%%%%% Plots:
t(end)=[];t=t';
figure;subplot(2,1,1);plot(t,x_K(1,:),'b','linewidth',3);hold on;
plot(t,x_F(1,:),'r','linewidth',3);
xlabel('Time (s)');ylabel('x_1');
legend('K - Method','F - Method');
title('First System State.');
subplot(2,1,2);plot(t,x_K(2,:),'b','linewidth',3);hold on;
plot(t,x_F(2,:),'r','linewidth',3);
xlabel('Time (s)');ylabel('x_2');
legend('K - Method','F - Method');
title('Second System State.');

t(end)=[];figure;plot(t,u_K(1,:),'b','linewidth',3);hold on;
plot(t,u_F(1,:),'r','linewidth',3);
xlabel('Time (s)');ylabel('Optimal Control Law Value');
legend('K - Method','F - Method');
title('Optimal Control Law.');


%%%%% RMS error:
RMS_x1=sqrt(mean((x_K(1,:)-x_F(1,:)).^2));
RMS_x2=sqrt(mean((x_K(2,:)-x_F(2,:)).^2));
RMS_u=sqrt(mean((u_K(1,:)-u_F(1,:)).^2));


%%%%% Tracking
%%%%% objective:-----------------------------------------------------------
t=0:ts:10;
C=[1 0; 0 1];H=[1 0;0 0];Q=[10 0;0 2];
r_tf=[sin(t(end));t(end)/2];r=[sin(t);t/2];
term_1=vpa(expand(0.5*(C*x-r_tf).'*H*(C*x-r_tf)),3);
term_2=vpa(expand(0.5*(C*x-r_tf).'*Q*(C*x-r_tf)),3);


%%%%% K solving:
K_track(:,:,N)=C'*H*C;
K_track_dot(:,:,N)=-C'*Q*C-K_track(:,:,N)*Ad-Ad'*K_track(:,:,N)+...
                    K_track(:,:,N)'*Bd*inv(R)*Bd'*K_track(:,:,N);
V_track(:,N)=C'*H*r_tf;
V_track_dot(:,N)=-(Ad'-K_track(:,:,N)'*Bd*inv(R)*Bd')*V_track(:,N)-C'*Q*r(:,N);
for k=1:1:N-1
    K_track(:,:,N-k)=K_track(:,:,N-k+1)-K_track_dot(:,:,N-k+1)*ts;
    K_track_dot(:,:,N-k)=-Q+K_track(:,:,N-k)*Bd*inv(R)*Bd'*K_track(:,:,N-k)-...
                         K_track(:,:,N-k)*Ad-Ad'*K_track(:,:,N-k);
    V_track(:,N-k)=V_track(:,N-k+1)-V_track_dot(:,N-k+1)*ts;
    V_track_dot(:,N-k)=-(Ad'-K_track(:,:,N-k)'*Bd*inv(R)*Bd')*...
                       V_track(:,N-k)-C'*Q*r(:,N-k); 
end


%%%%% Euler discretization:
x_Track(:,1)=[0;0];
for k=1:1:N-1
    u_Track(k)=-inv(R)*B'*(K(:,:,k)*x_Track(:,k)-V_track(:,k));
    x_Track(:,k+1)=([0 ts;-1*ts 2*ts]+eye(2))*x_Track(:,k)+...
    [0;ts]*u_Track(k);
end


%%%%% RMS error:
RMS_x1=sqrt(mean((x_K(1,:)-r(1,1:end-1)).^2));
RMS_x2=sqrt(mean((x_K(2,:)-r(2,1:end-1)).^2));


%%%%% Plots:
t(end)=[];t=t';r(:,end)=[];
figure; hold on;
ax1 = subplot(2,1,1);plot(t,r(1,:),'linewidth',2); hold on; plot(t,x_Track(1,:),'linewidth',2); 
ax2 = subplot(2,1,2);plot(t,r(2,:),'linewidth',2); hold on; plot(t,x_Track(2,:),'linewidth',2);
title(ax1, 'State $x_1$','Interpreter','latex','FontSize',14); 
title(ax2, 'State $x_2$','Interpreter','latex','FontSize',14);   
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);

t(end)=[];figure;plot(t,u_Track(1,:),'linewidth',2);
xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
ylabel('Optimal Control Law', 'Interpreter','Latex', 'FontSize',14);
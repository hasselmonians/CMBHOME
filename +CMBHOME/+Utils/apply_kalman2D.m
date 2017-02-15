function [x,V]=apply_kalman2D(ph,Q,R)
% addpath(genpath([pwd '\KalmanAll']));
% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)
% ph=sin(1:.1:100)+.5*(rand(1,991)-.5);
ss = 6; % state size
os = 2; % observation size
F = [1 1 .5 0 0 0; 
     0 1  1 0 0 0;
     0 0  1 0 0 0;
     0 0  0 1 1 .5;
     0 0  0 0 1 1;
     0 0  0 0 0 1];  %the system matrix
H = [1 0 0 0 0 0;
     0 0 0 1 0 0];   % the observation matrix 
% Q = [.1 0 0;
%       0 .1 0;
%       0 0 .01];%.1*eye(ss);  % Q - the system covariance 
% R = 2;% reliability of measure/ the observation covariance
initx = 1*[ph(1,1); ph(1,3)-ph(1,1); ph(1,3)+ph(1,1)-2*ph(1,2); ph(2,1); ph(2,3)-ph(2,1); ph(2,3)+ph(2,1)-2*ph(2,2)];
initV = 5*eye(ss);

dims=size(Q,1)/2;
F=F([(1:dims) 3+(1:dims)],[ (1:dims) 3+(1:dims)]);
H=H(:,[(1:dims) 3+(1:dims)]);
initx=initx([(1:dims) 3+(1:dims)]);
initV=initV([(1:dims) 3+(1:dims)],[ (1:dims) 3+(1:dims)]);


[x, V] = CMBHOME.Utils.kalman_smoother(ph, F, H, Q, R, initx, initV);
% plot(x(1,:)); hold on;
% plot(5*x(2,:),'g');
% plot(50*x(3,:),'m');
% plot(ph,'.k');
% plot(sin(1:.1:100),'r--'); 
% plot(mean(squeeze(mean(V))),'c');hold off;
% axis([0 991 -2 2]); grid on;



function p=speed_kalman(x,y, fs_video)
defpars=[50 .02 .01 100 2 .1 10];
if ~exist('pars','var')
    pars=defpars;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%% pos estimation
%     pars=[100 1 .03 1000];
%     pars=[100 50 1 1000];
% pars=[108.2592   52.2000    0.9360  762.8400];
    
    R=pars(4)*eye(2);
    Q=zeros(6);
    Q(1,1)=pars(1); Q(4,4)=pars(1);
    Q(2,2)=pars(2); Q(5,5)=pars(2);
    Q(3,3)=pars(3); Q(6,6)=pars(3);
    [out,V]=CMBHOME.Utils.apply_kalman2D([x';y'],Q,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.x=out(1,:)';
p.y=out(4,:)';
% p.vx=50*out(2,:)';
% p.vy=50*out(5,:)';
% p.ax=50^2*out(3,:)';
% p.ay=50^2*out(6,:)';
p.vx=fs_video*out(2,:)';
p.vy=fs_video*out(5,:)';
p.ax=fs_video^2*out(3,:)';
p.ay=fs_video^2*out(6,:)';
p.a=sqrt(p.ax.^2+p.ay.^2);
p.v=sqrt(p.vx.^2+p.vy.^2);
p.vd=atan2(p.vy,p.vx);
p.ad=atan2(p.ay,p.ax);
p.dsp=(p.vx.*p.ax+p.vy.*p.ay)./p.v;     %The derivative of speed
% p.vx_err=50*sqrt(squeeze(V(2,2,:)));
% p.vy_err=50*sqrt(squeeze(V(5,5,:)));
% p.ax_err=50^2*sqrt(squeeze(V(3,3,:)));
% p.ay_err=50^2*sqrt(squeeze(V(6,6,:)));

%%%%%%%%% more complex acceleration quantities
p.a_par_vel=p.a.*(cos(p.vd-p.ad));
p.a_par_vel_pos=p.a_par_vel.*(p.a_par_vel>0);


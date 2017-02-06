function out = Kropff(self)
%Uses Kalman filter from Kropff et al. 2015  Velocity output is stored as
%out.v
import CMBHOME.*

xvec = CMBHOME.Utils.nanInterp(self.b_x,'spline');
yvec = CMBHOME.Utils.nanInterp(self.b_y,'spline');

out = CMBHOME.Session.speed_kalman(xvec, yvec, self.fs_video);

end
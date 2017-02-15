function [xsmooth, Vsmooth, VVsmooth, loglik] = kalman_smoother(y, A, C, Q, R, init_x, init_V, varargin)
% Kalman/RTS smoother.
% [xsmooth, Vsmooth, VVsmooth, loglik] = kalman_smoother(y, A, C, Q, R, init_x, init_V, ...)
%
% The inputs are the same as for kalman_filter.
% The outputs are almost the same, except we condition on y(:, 1:T) (and u(:, 1:T) if specified),
% instead of on y(:, 1:t).

[os T] = size(y);
ss = length(A);

% set default params
model = ones(1,T);
u = [];
B = [];

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i}
   case 'model', model = args{i+1};
   case 'u', u = args{i+1};
   case 'B', B = args{i+1};
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

xsmooth = zeros(ss, T);
Vsmooth = zeros(ss, ss, T);
VVsmooth = zeros(ss, ss, T);

% Forward pass
[xfilt, Vfilt, VVfilt, loglik] = CMBHOME.Utils.kalman_filter(y, A, C, Q, R, init_x, init_V, ...
					       'model', model, 'u', u, 'B', B);

% Backward pass
xsmooth(:,T) = xfilt(:,T);
Vsmooth(:,:,T) = Vfilt(:,:,T);
%VVsmooth(:,:,T) = VVfilt(:,:,T);

for t=T-1:-1:1
  m = model(t+1);
  if isempty(B)
      [xsmooth(:,t), Vsmooth(:,:,t), VVsmooth(:,:,t+1)] = ...
        CMBHOME.Utils.smooth_update(xsmooth(:,t+1), Vsmooth(:,:,t+1), xfilt(:,t), Vfilt(:,:,t), ...
                  Vfilt(:,:,t+1), VVfilt(:,:,t+1), A(:,:,m), Q(:,:,m), [], []);
%       if ~isnan(xsmooth(:,t+1))
%         [xsmooth(:,t), Vsmooth(:,:,t), VVsmooth(:,:,t+1)] = ...
%         smooth_update(xsmooth(:,t+1), Vsmooth(:,:,t+1), xfilt(:,t), Vfilt(:,:,t), ...
%                   Vfilt(:,:,t+1), VVfilt(:,:,t+1), A(:,:,m), Q(:,:,m), [], []);
%       else
%           n = 2;
%           while isnan(xsmooth(:,t+n)) & t+n<=T
%               n = n+1;
%           end
%           [xsmooth(:,t), Vsmooth(:,:,t), VVsmooth(:,:,t+n)] = ...
%         smooth_update(xsmooth(:,t+n), Vsmooth(:,:,t+n), xfilt(:,t), Vfilt(:,:,t), ...
%                   Vfilt(:,:,t+n), VVfilt(:,:,t+n), A(:,:,m), Q(:,:,m), [], []);
%       end
  else
    [xsmooth(:,t), Vsmooth(:,:,t), VVsmooth(:,:,t+1)] = ...
	CMBHOME.Utils.smooth_update(xsmooth(:,t+1), Vsmooth(:,:,t+1), xfilt(:,t), Vfilt(:,:,t), ...
		      Vfilt(:,:,t+1), VVfilt(:,:,t+1), A(:,:,m), Q(:,:,m), B(:,:,m), u(:,t+1));
  end
end

VVsmooth(:,:,1) = zeros(ss,ss);


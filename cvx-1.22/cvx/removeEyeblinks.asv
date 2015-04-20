function[xout]=removeEyeblinks(x,xeye)
%This function computes the forward matrix A that describes the component
%of the EEG that is produced by eyeblinking
%:::INPUT:::
%xeye      a 48xT matrix that is the EEG output of all 48 channels over times t where *eyeblinking was the predominant activity*
%x         a 48xT2 matrix that is the total EEG output that you would like to remove eyeblinking effects from. Here T2>T
%:::OUTPUT:::
%xout      a 48xT2 matrix that is the EEG output cleaned of the effects of eyeblinking

%Compute the blinking covariance
xavg=mean(xeye,1); 
R=(xeye-xavg)*(xeye-xavg)'; %(48x48)
sqrtR=sqrtm(R);

%Solve for forward model 
cvx_begin
variable w(48,1);
minimize( norm(sqrtR*w) );
subject to norm(w)<=1;          % Shouldn't norm(w) be = 1?
cvx_end
a=w/norm(w)^2; %(48x1)
ainv=inv(a'*a)*a';

%clean the EEG data in x by removing the effect of eye blinking
xout=(eye(48)-a*ainv)*x;









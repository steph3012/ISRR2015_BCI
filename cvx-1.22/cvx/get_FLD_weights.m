function[w]=get_FLD_weights(X,d)
%This function computes the forward matrix A that describes the component
%of the EEG that is produced by eyeblinking
%:::INPUT:::
%X       a 48xNxT matrix that is the EEG output of all 48 channels over N
%        total trials (total=correct+incorrect) where each trial has T time samples.
%d       a Nx1 vector of labels {0,1} for each trial where 0=correct,
%        1=incorrect
%:::OUTPUT:::
%w       a 48x1 vector of weights for each channel 

N1=sum(d==0);   %number of correct trials
N2=sum(d==1);   %number of incorrect trials

%Compute the cumulative difference between incorrect and correct trials
Xcorr=X(:,d==0,:);
Xincorr=X(:,d==1,:);
delx=(1/N1)*sum(Xcorr,2)-(1/N2)*sum(Xincorr,2);    %(48xT)
delx_bar=sum(delx,2);                              %should return a 48x1 vector

%Compute within-class covariances 
Xcorr_mean=mean(mean(Xcorr,2),2);       %(48x1)
Xincorr_mean=mean(mean(Xincorr,2),2);   %(48x1)
Rcorr=0; Rincorr=0;
for(n=1:N)
    for(t=1:T)
        Rcorr=Rcorr+(X(:,n,t)-Xcorr_mean)*(X(:,n,t)-Xcorr_mean)'; 
        Rincorr=Rincorr+(X(:,n,t)-Xincorr_mean)*(X(:,n,t)-Xincorr_mean)';
    end
end
w=inv(R1+R2)*delx_bar;










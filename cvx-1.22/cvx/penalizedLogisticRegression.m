function [w]=penalizedLogisticRegression(X,d)
%This function performs the iterations described in Table 1 of 
%Parra et. al. 2005 in order to return a vector of learned weights w for
%online classification of EEG signals.  
%:::Inputs:::
%X     a (48xN) matrix where each column is a 48-channel EEG trace of a
%      "correct" or "incorrect" baxter iteration
%d     a (Nx1) vector of labels {0,1} where a 1 in element i means that the
%      ith column of X is the EEG trace of a correct trial and vice versa
%:::Outputs:::
%w     a (48x1) vector of weights that is learned from the data
%NOTE: The returned w is for one sample in time! This approach assumes that
%all samples are i.i.d. over time.  Therefore each column are the
%48-channel outputs of the EEG signal at a specific instant in time t for
%different labeled trials (either correct or incorrect).

X(49,:)=ones(1,N);
w=zeros(49,1);
tol=1e-4;
diff=inf;
Lambda=[eye(48) 1; ones(1,49)];
while(diff>tol)
    p=1/(1+exp(-w'*X));
    g=X*(d-p)-Lambda*w;
    H=X*diag(p.*(1-p))*X'+Lambda;
    wlast=w;
    w=w+inv(H)*g;
    diff=sqrt((wlast-w)'*(wlast-w));
end

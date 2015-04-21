function[w]=l1LogisticRegression_BCI_ISRR(X,d)
%This function performs l1 penalized logistic regression using the toolbox
%"L1 General" available online in order to return a vector of learned weights w for
%online classification of EEG signals.  
%:::Inputs:::
%X     a (49xN) matrix where each column is a 48-channel EEG trace of a
%      "correct" or "incorrect" baxter iteration
%d     a (Nx1) vector of labels {-1,1} where a 1 in element i means that the
%      ith column of X is the EEG trace of a correct trial and vice versa
%:::Outputs:::
%w     a (49x1) vector of weights that is learned from the data where the
%      last element is b (y-intercept)
%NOTE: The returned w is for one sample in time! This approach assumes that
%all samples are i.i.d. over time.  Therefore each column are the
%48-channel outputs of the EEG signal at a specific instant in time t for
%different labeled trials (either correct or incorrect).
addpath(genpath(pwd))
X(49,:)=ones(1,N);  %add row for y-intercept
tol=1e-4;
diff=inf;
Lambda=[eye(48) 1; ones(1,49)];
funObj = @(w)LogisticLoss(w,X,d);
w_init = zeros(49,1);
fprintf('\nComputing L1-Regularized Logistic Regression Coefficients...\n');
w = L1General2_PSSgb(funObj,w_init,lambda);
stem(wLogL1,'g');
xlim([1 49]);
title('L1-Regularized Logistic Regression');

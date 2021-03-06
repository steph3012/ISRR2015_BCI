addpath(genpath(pwd))

%Generate some data
nInstances = 250;
nVars = 50;
X = randn(nInstances,nVars);
y = X*((rand(nVars,1) > .5).*randn(nVars,1)) + randn(nInstances,1);
%Least Squares Solution
wLS = X\y;

%Ridge Regression
lambda = 0.5*ones(nVars,1); % Penalize each element by the same amount
R = chol(X'*X + diag(lambda));
wRR = R\(R'\(X'*y));

X = [ones(nInstances,1) X]; % Add Bias element to features
lambda = 10*ones(nVars+1,1); % Penalize each element by the same amount
y = sign(y); % Convert y to binary {-1,1} representation

funObj = @(w)LogisticLoss(w,X,y);
w_init = zeros(nVars+1,1);

% Maximum Likelihood
fprintf('\nComputing Maximum Likelihood Logistic Regression Coefficients\n');
mfOptions.Method = 'newton';
wLogML = minFunc(funObj,w_init,mfOptions);

%L1-Regularized Logistic Regression
fprintf('\nComputing L1-Regularized Logistic Regression Coefficients...\n');
wLogL1 = L1General2_PSSgb(funObj,w_init,lambda);

figure;
clf;hold on;
subplot(2,2,1);
stem(wLogML,'r');
xlim([1 nVars+1]);
title('Maximum Likelihood Logistic Regression');
% subplot(2,2,2);
% stem(wLogL2,'b');
% xlim([1 nVars+1]);
% title('L2-Regularized Logistic Regression');
subplot(2,2,3);
stem(wLogL1,'g');
xlim([1 nVars+1]);
title('L1-Regularized Logistic Regression');
% subplot(2,2,4);
% stem(wLogL1L2,'c');
% xlim([1 nVars+1]);
% title('Elastic-Net Logistic Regression');

fprintf('Number of Features Selected by Maximum Likelihood Logistic Regression classifier: %d (out of %d)\n',nnz(wLogML(2:end)),nVars);
%fprintf('Number of Features Selected by L2-regualrized Logistic Regression classifier: %d (out of %d)\n',nnz(wLogL2(2:end)),nVars);
fprintf('Number of Features Selected by L1-regualrized Logistic Regression classifier: %d (out of %d)\n',nnz(wLogL1(2:end)),nVars);
%fprintf('Number of Features Selected by Elastic-Net Logistic Regression classifier: %d (out of %d)\n',nnz(wLogL1L2(2:end)),nVars);
fprintf('Classification error rate on training data for L1-regularied Logistic Regression: %.2f\n',sum(y ~= sign(X*wLogL1))/length(y));
%pause;
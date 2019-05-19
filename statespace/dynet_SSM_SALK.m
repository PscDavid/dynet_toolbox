function SALK = dynet_SSM_SALK(Y,p,lambda)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The Sparse Adaptive Least-squares Kalman filter with self-tuning memory
%                                       D. Pascucci, University of Fribourg
%                                       M. Rubega,   University of Geneve
% Last update: 16.05.2019
% updated 06.05.2019 by Maria Rubega    -SVD
% updated 09.05.2019 by David Pascucci  -SVD+norm max lambda+self-tuning c
%--------------------------------------------------------------------------
% INPUT:
% -Y:      matrix
%           [Trials X Channels X Time-samples]
%           Sample data
% -p:      scalar positive
%           p-order
% -lambda: scalar positive,        if lambda 0 -> no regularization
%           Regularization parameter, default 0.001
%--------------------------------------------------------------------------
% OUTPUT:   SALK, structure with fields:
% -AR:     matrix
%           [Channels X Channels X p X Time-samples]
%           MVAR model coefficients
% -R:     Measurement Noise Covariance Matrix (innovation based)
%           [Channels X Channels X Time]
% -PY:     matrix
%           [Trials X Channels X Time-samples]
%           One-step PredictioR on Y
% -PYe:    matrix
%           [Trials X Channels X Time-samples]
%           One-step Residuals
% -c:      self-tuning c for each k
%==========================================================================
% References:
% [1] Nilsson, M. (2006). Kalman filtering with unknown noise covariances.
%     In Reglermöte 2006.
% [2] Bach, F., Jenatton, R., Mairal, J., & OboziRki, G. (2011). 
%     Convex optimization with sparsity-inducing norms.
%     Optimization for Machine Learning, 5, 19-53.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -check input
if numel(size(Y))<3;error('Check the dimeRioR of Y');end
if nargin<3;lambda   = 0.001;end

% -Preallocating main variable
[trl,dim,tm] = size(Y);
AR           = zeros(dim,dim*p,tm);          % AR matrix
R           = zeros(dim,dim,tm);             % Additive observation noise
PY           = zeros(size(Y));               % One-step predictioR

% -Initializing variables
xm           = zeros(dim*p,dim)+1e-4;        % Prior state estimates
trEk         = zeros(tm,1);                  % Innovation monitoring
allc         = zeros(tm,1);                  % Self-tuning c, for all k

% -Loop from p(lag)+1 through time
for k = (p+1):tm
    
    % TraRition matrix H and observation z at time k
    data_back  = squeeze(Y(:,:,k-1:-1:k-p));
    H          = reshape(data_back,trl,dim*p);
    
    % Measurement and Innovation
    Z          = Y(:,:,k);
    PY(:,:,k)  = H*xm;
    vk         = Z - PY(:,:,k); % Innovation at time k

    % Measurements innovation monitoring
    tmp        = (vk'*vk);
    R(:,:,k)   = tmp./max(trl-1,1);
    trEk(k,:)  = trace(tmp);

    % SVD Tikhonov (Spectral decomposition)
    % - determine lambda
    lambda_max = norm(H'*Z,'inf'); % convex rule [2]
    lambda_k   = lambda*lambda_max;
    % economy-size decomposition of trl-by-dim H
    [U,S,V]    = svd(H,'econ');
    % only the first dim column of U are computed, and S is dim-by-dim
    d          = diag(S);
    D          = diag(d./(d.^2+lambda_k));
    Hinv       = V*D*U';
    betas      = Hinv*Z;
    
    % self-tuning adaptation constant
    if k>(p+1)*2 % 0.05<=c<=0.95
        c      = min(0.05+abs((mean(trEk(k:-1:k-p,:))-       ...
                            mean(trEk(k-p-1:-1:k-p*2,:)))    ...
                         ./(mean(trEk(k-p-1:-1:k-p*2,:)))),0.95);
    else
        c      = 0.05;
    end
    allc(k,:)  = c;

    % Kalman update
    xp         = (xm + c*betas) ./ (1+c);  % [1]
    xm         = xp;
    AR(:,:,k)  = xm';
       
end

% -Saving output variables
AR            = reshape(AR,[dim dim p tm]);
SALK.AR       = AR;
SALK.R        = R;
SALK.PY       = PY;
SALK.c        = allc; 


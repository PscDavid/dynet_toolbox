function STOK = dynet_SSM_STOK(Y,p,ff)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The Self-Tuning optimized Kalman filter
%                                       D. Pascucci, University of Fribourg
%                                       M. Rubega,   University of Geneve
% Last update: 23.08.2019
%--------------------------------------------------------------------------
% INPUTs:
% -Y:      matrix
%           [Trials X Channels X Time]
%           Sample data
% -p:      scalar positive
%           p-order
% -ff:     percentage of variance explained for setting the filter factor
%           Regularization parameter, default 0.99
%--------------------------------------------------------------------------
% OUTPUT:   STOK, structure with fields:
% -AR:     matrix
%           [Channels X Channels X p X Time]
%           MVAR model coefficients
% -R:      Measurement Noise Covariance Matrix (innovation based)
%           [Channels X Channels X Time]
% -PY:     matrix
%           [Trials X Channels X Time]
%           One-step PredictioR on Y
% -c:      self-tuning c for each k
% -FFthr   filtering factor threshold for each k [2]
%==========================================================================
% References:
% [1] Nilsson, M. (2006). Kalman filtering with unknown noise covariances.
%     In Reglermöte 2006.
% [2] Hansen, P. C. (1987). The truncated svd as a method for 
%     regularization. BIT Numerical Mathematics, 27(4), 534-553.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -check input
if numel(size(Y))<3;error('Check the dimensions of Y');end
if nargin<3;ff   = 0.99;end

% -Preallocating main variable
[trl,dim,tm] = size(Y);
AR           = zeros(dim,dim*p,tm);          % AR matrix
R            = zeros(dim,dim,tm);            % Additive observation noise
PY           = zeros(size(Y));               % One-step predictioR

% -Initializing variables
xm           = zeros(dim*p,dim)+1e-4;        % Prior state estimates
trEk         = zeros(tm,1);                  % Innovation monitoring
allc         = zeros(tm,1);                  % Self-tuning c, for all k

% -Loop from p(lag)+1 through time
for k = (p+1):tm
    
    % Transition matrix H and observation z at time k
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
    % economy-size decomposition of trl-by-dim H
    [U,S,V]    = svd(H,'econ');
    % only the first dim column of U are computed, and S is dim-by-dim
    d          = diag(S);
    if isnumeric(ff)
    % determine filtering factor threshold [2]
    relv       = d.^2./sum(d.^2);
    filtfact   = find(double(cumsum(relv)<ff),1,'last');
    if isempty(filtfact)
        filtfact = 1;
    end
    lambda_k   = d(filtfact).^2;
    STOK.FFthr(k,:) = filtfact;
    %diag(1./d.*((d.^2)./(d.^2+lambda_k)));
    D          = diag(d./(d.^2+lambda_k)); 
    else
    D          = diag(1./d); 
    STOK.FFthr(k,:) = NaN;
    end
    
    Hinv       = V*D*U';
    betas      = Hinv*Z;
        
    % self-tuning adaptation constant
    if k>(p+1)*2 % 0.05<=c<=0.95
        ntrEk  = trEk(k-1:-1:(k-p*2));
        e_k    = mean(ntrEk(1:p));
        e_p    = mean(ntrEk(p+1:end));
        x      = abs(e_k-e_p)./e_p;
        c      = min(0.05 + x ,0.95); % or  x./max percentage allowed
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
STOK.AR       = AR;
STOK.R        = R;
STOK.PY       = PY;
STOK.c        = allc; 


function PDC = dynet_ar2pdc(KF,srate,freqs,measure,univ,flow,PSD)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Obtain PDC, sPDC, info-PDC from tvAR coeffients
%                                          M.Rubega, D.Pascucci, 17.10.2018
% Last update: 22.10.2019
%--------------------------------------------------------------------------
% INPUTs
% - KF:        Output from dynet_SSM_KF or dynet_SSM_STOK
%              containing the matrix of AR coefficients
%              [n x n x order x time]
%              and the estimated measurement noise covariance matrix R
%              [n x n x time]
% - srate:     Sampling rate
% - freqs:     Frequency vector (column)
% - metric:    see OUTPUT
% - univ:      Remove (0, default) or Keep (1) the diagonal elements
% - flow:      normalization per columns (1) or rows (2)
% - PSD:       (1) Add the normalized parametric PSD on diagonals
%              (0) none (only for graphical purpose)
%--------------------------------------------------------------------------
% OUTPUTs
% - PDC:       [Nodes X Nodes X Freq X Time]
%              one of
%             'PDC'     % non-squared       Eq. (18) in [2]
%             'sPDC'    % squared           Eq. (7) in [3]
%             'PDCnn'   % non-normalized
%             'sPDCnn'  % squared non-normalized
%             'iPDC'    % info-PDC          Eq. (5.25) in [4]
%             'iPDCs'   % info-PDC squared
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m; dynet_parspsd.m 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% References:
% [1] Milde, T., Leistritz, L., ..., & Witte, H. (2010), Neuroimage, 50(3),
%     960-969. A new Kalman filter approach for the estimation of
%     high-dimensional time-variant multivariate AR models and its
%     application in analysis of laser-evoked brain potentials.
% [2] Baccalá, L. & Sameshima, K. (2001) Biol Cybern, 84 (6), 463–474
%     Partial directed coherence: a new concept in neural structure
%     determination.
% [3] Astolfi, L., ..., & Babiloni, F. (2006), IEEE Transactions on
%     Biomedical Engineering, 53(9), 1802-1812
%     Assessing cortical functional connectivity by partial directed
%     coherence: simulations and application to real data.
% [4] K. Sameshima and L. A. Baccala, (2014), CRC Press
%     Methods in brain connectivity inference through multivariate time
%     series analysis.
% [5] Toppi, J., ..., & Astolfi, L. (2013), 35th Annual International
%     Conference of the IEEE EMBS, 4346-4349
%     The Effect of Normalization of Partial Directed Coherence on the
%     Statistical Assessment of Connectivity Patterns: A Simulation Study
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% check input
default('measure','sPDC');
default('univ',0);          
default('flow',1);          % normalization per columns
default('PSD',0);

[nodes,~,order,time] = size(KF.AR);

% check R
if isfield(KF,'R') && numel(size(KF.R))<3
    KF.R       = repmat(KF.R,[1 1 time]);
elseif ~isfield(KF,'R')
    KF.R       = repmat(eye(nodes),[1 1 time]);
end

% -Transform AR (time-domain) into A(f) (frequency-domain)
%  freqs should be a column vector
if size(freqs,2)>size(freqs,1); freqs=freqs'; end

Z        = exp(-2*pi*1i*freqs/srate).^(1:order);
A        = repmat(eye(nodes), [1 1 length(freqs) time]);
for k = 1:order
    tmp  = repmat(-KF.AR(:, :, k,:), [1 1 length(freqs) 1]);
    A    = A + bsxfun(@times,tmp,reshape(Z(:,k),1,1,[]));
end

% -Metrics
switch measure
    
    case 'PDC' % Eq. (18) in [2]
        
        PDC = bsxfun(@rdivide,abs(A), sqrt(sum(abs(A).^2,flow)));
        
    case 'sPDC' % Eq. (7) in [3]
        
        PDC = bsxfun(@rdivide,abs(A).^2, (sum(abs(A).^2,flow)));
        
    case 'PDCnn'
        PDC = abs(A);
        
    case 'sPDCnn'
        PDC = abs(A).^2;  
        
    case 'iPDC' % Eq. (5.25) in [4] -> if R is diagonal, 'generalized PDC'
        
        % For construction, R should be costant over time:
        R         = mean(KF.R(:,:,max(round(end/2),order):end),3);
        SIGMA     = pinv(R);
        w_ii      = repmat(diag(R),[1 nodes length(freqs) time]);
        
        den_2     = zeros(1,nodes,length(freqs),time);
        for j = 1:nodes
            a_j   = A(:,j,:,:);
            den_1 = permute(sum(bsxfun(@times,repmat(a_j,[1 nodes 1 1]), ...
                SIGMA),1),[2 1 3 4]);
            den_2(1,j,:,:) = sum(bsxfun(@times,den_1,conj(a_j)));
        end
        den       = repmat(den_2,[nodes 1 1 1]);
        PDC       = abs(bsxfun(@rdivide,A,sqrt( bsxfun(@times,w_ii,den) )));
        
    case 'iPDCs' % Eq. (5.25) in [4] -> if R is diagonal, 'generalized PDC'
        
        % For construction, R should be costant over time:
        R         = mean(KF.R(:,:,max(round(end/2),order):end),3);
        SIGMA     = pinv(R);
        w_ii      = repmat(diag(R),[1 nodes length(freqs) time]);
        
        den_2     = zeros(1,nodes,length(freqs),time);
        for j = 1:nodes
            a_j   = A(:,j,:,:);
            den_1 = permute(sum(bsxfun(@times,repmat(a_j,[1 nodes 1 1]), ...
                SIGMA),1),[2 1 3 4]);
            den_2(1,j,:,:) = sum(bsxfun(@times,den_1,conj(a_j)));
        end
        den       = repmat(den_2,[nodes 1 1 1]);
        PDC       = abs(bsxfun(@rdivide,abs(A).^2,( bsxfun(@times,w_ii,den) )));
end

% -Remove or keep the autoregressive component
if univ==0
    for dg = 1:size(PDC,1)
        PDC(dg,dg,:,:) = NaN;
    end
end

% -If required, store the parametric PSD on the diagonal
%  ( only for plotting purposes )
if PSD
    KF              = dynet_parpsd(KF,srate,freqs,2);
    pPSD            = abs(KF.SS).^2;
    [~,d,f,t]       = size(PDC);
    dg              = repmat(eye(d),[1 1 f t]);
    pPSD            = (pPSD-min(pPSD(:)))./range(pPSD(:))...
        *range(PDC(dg==0))+min(PDC(dg==0));
    for k = 1:d
        PDC(k,k,:,:)= pPSD(k,k,:,:);
    end
end


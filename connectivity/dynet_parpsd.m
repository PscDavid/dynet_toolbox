function KF = dynet_parpsd(KF,srate,f_range,SSout,t_win)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Get A(f) and SS(optional) from AR coefficients
% adds frequency-related info on the Kalman filter estimates
%                                       D. Pascucci, University of Fribourg
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INPUTs
% - KF:         Output from SSM_KF/SSM_STOK
% - srate:      Sampling rate
% - f_range:    Frequency range for A(f) and SS computation, column vector
%               - linear spacing if size(f_range)==2
%               - user-specified spacing if size(f_range)>2
% - SSout:      (0) A(f),  (1) A (f ) and SS, (2) SS
% - t_win:      Temporal window (in samples) to restrict A(f) and SS
%--------------------------------------------------------------------------
% OUTPUTs (additional fields in KF)
% - AF:         A(f) matrix
%               [Channels X Channels X Freq X Time]
% - SS:         Spectral-Autospectral matrix (optional)
%               [Channels X Channels X Freq X Time]
% - AuS:        Autospectral matrix (optional)
%               [Channels X Freq X Time]
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -check input
PAR          = KF.AR;
if isfield(KF,'R')
    R        = KF.R;
    if numel(size(KF.R))==2
        R    = repmat(KF.R,[ 1 1 size(KF.AR,4)]);
    end
elseif isfield(KF,'R')
    R        = KF.R;
end
if numel(f_range)>2; freqs = f_range;
else; freqs  = linspace(f_range(1), f_range(end), diff(f_range')+1)';
end
if nargin==5
    PAR      = PAR(:,:,:,t_win);
    R        = R(:,:,t_win);
    KF.t_win_f   = t_win;
end
[nodes,~,order,time] = size(PAR);

% -append freq info to KF main struct
KF.srate    = srate;
KF.freqs     = freqs;

% -Z-transform -> A(f)
Z            = exp(-2*pi*1i*freqs/srate).^(1:order);
A            = repmat(eye(nodes), [1 1 numel(freqs) time]);
for k = 1:order
    tmp      = repmat(-PAR(:, :, k,:), [1 1 numel(freqs) 1]);
    A        = A + bsxfun(@times,tmp,reshape(Z(:,k),1,1,[]));
end
KF.Af        = A;

% -SS
S_ft         = median(R(:,:,round(time/2):end),3);
if SSout > 0
    SS       = NaN(nodes,nodes,numel(freqs),time);
    for t = 1:time
        for f = 1:numel(freqs)
            A_ft        = squeeze(A(:,:,f,t));
            SS(:,:,f,t) = A_ft\S_ft/A_ft';
        end
    end
    KF.SS        = SS;
    if SSout < 2
        % -AuS
        AuS      = NaN(nodes,numel(freqs),time);
        for k = 1:nodes
            AuS(k,:,:) = SS(k,k,:,:);
        end
        KF.AuS       = AuS;
    end
end

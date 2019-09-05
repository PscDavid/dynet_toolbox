function [Y,noise] = addnoise(Y,SNRdb,type)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Add white or 1/f noise
%
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INPUTs
% - Y:      matrix
%           [Trials X Channels X Time]
%           Noiseless data
% - SNRdb   positive scalar
%           signal-to-noise ratio in db
% - type    'w' (white noise) or '1/f' pink noise
%--------------------------------------------------------------------------
% OUTPUTs
% - Y:      matrix
%           [Trials X Channels X Time]
%           Noisy data
% - noise   additive noise
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

As          = std(Y(:));
Aw          = As/10^((SNRdb/20));
switch type
    case 'w'    
        noise        = Aw.*randn(size(Y));
        Y            = Y + noise;
    case '1/f'
        noise        = Aw.*randnd(-1,size(Y));
        Y            = Y + noise;
end            


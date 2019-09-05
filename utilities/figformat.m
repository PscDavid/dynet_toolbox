function figformat(h,v,scale,fsize)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define default figure format 
%
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INPUTs
% - h:     default 0 otherwise 1, draw horizontal black line in inf,0
% - v:     default 0 otherwise 1, draw vertical black line in 0,inf
% - scale: scale
% - fsize: font size
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
default('h',0);
default('v',0);
default('fsize',15);
set(gca,'fontsize',fsize,'fontname','Helvetica',...
    'linewidth',1.5,'TickDir','out')
set(gcf,'color','w')
if h
    hline(0,'k--');
end
if v
    vline(0,'k--');
end
% scale  = 0.1;
if nargin==3
    pos    = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
end

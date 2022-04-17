function figformat(h,v,scale,fsize)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define default figure format 
%
% Last update: 31.03.2021
%--------------------------------------------------------------------------
% INPUTs
% - h:     default 0 otherwise x, draw horizontal black line in inf,x
% - v:     default 0 otherwise x, draw vertical black line in x,inf
% - scale: scale
% - fsize: font size
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
default('h',0);
default('v',0);
default('scale',NaN);
default('fsize',15);
% set(gca,'fontsize',fsize,'fontname','Helvetica',...
%     'linewidth',1.5,'TickDir','out')
set(gca,'fontsize',fsize,'fontname','Helvetica',...
    'linewidth',1,'TickLength',[0 0])
set(gcf,'color','w')
if ~isnan(h)
    hline(h,'k--');
end
if ~isnan(v)
    vline(v,'k--');
end
% scale  = 0.1;
if nargin==3 && ~isnan(scale)
    pos    = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
end

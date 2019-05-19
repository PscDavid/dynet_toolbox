function figformat(h,v,scale,fsize)

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

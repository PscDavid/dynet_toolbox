function formatlegend(pl,labels,lw,fsz)
default('lw',5);
default('fsz',20);
[~, hobj, ~, ~] = legend(pl,labels,'box','off');
hl = findobj(hobj,'type','line');set(hl,'LineWidth',lw);
ht = findobj(hobj,'type','text');set(ht,'FontSize',fsz);

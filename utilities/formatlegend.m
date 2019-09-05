function formatlegend(pl,labels,lw,fsz)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define default legend figure format 
%
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INPUTs
% - pl,labels: legend labels
% - lw:        line width
% - fsz:       font size
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% default.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
default('lw',5);
default('fsz',20);
[~, hobj, ~, ~] = legend(pl,labels,'box','off');
hl = findobj(hobj,'type','line');set(hl,'LineWidth',lw);
ht = findobj(hobj,'type','text');set(ht,'FontSize',fsz);

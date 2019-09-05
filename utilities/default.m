function default(argname,def_val)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Confirm existence of 'argname' in the ws and assign the value of
% 'def_val'
%
% Last update: 22.08.2019
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
varexist    = evalin('caller',['exist(''' argname ''',''var'')']);
makedefault = 0;
if ~varexist
    makedefault = 1;
else
    varisempty = evalin('caller',['isempty(' argname ')']);
    if varisempty
        makedefault = 1;
    end
end
if makedefault
    assignin('caller',argname,def_val);
end
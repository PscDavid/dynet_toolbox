function shuffled_x = shuffling(x)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
% Random permutation of the vector x
%
% Last update: 22.08.2019
%--------------------------------------------------------------------------
% INVOKED FUNCTIONs
% randperms.m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

shuffled_x =  x(randperms(length(x))); 
% MR: shuffled_x =  x(randperm(length(x)));
end
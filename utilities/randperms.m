function p = randperms(varargin)
%==========================================================================
% Random permutation by Po Hu 11.08.2005
%     
% Based on RANDPERM                                         %
%   RANDPERM(n) is a random permutation of the integers from 1 to n.  %
%   For example, RANDPERM(6) might be [2 4 5 6 1 3].                  %    
%   Copyright 1984-2004 The MathWorks, Inc.                           %
%   $Revision: 5.10.4.1 $  $Date: 2004/03/02 21:48:27 $               %
%--------------------------------------------------------------------------
% Description:
%     randperms(n) is a random permutation of the integers from 1 to n.
%     randperms(m,n) when m > n, returns the nth element of randperm(m);
%                    when m < n, returns the permutation of the intergers
%                    from m to n.
%     randperms(m,n,idx) returns the idx th element of permutated intergers
%                        from m to n. 
% For example:
%     randperms(6)  might be [2 4 5 6 1 3].
%     randperms(6,1) returns 1st element of randperm(6) might be 2.
%     randperms(6,1:3) returns 1st-3rd elements of randperm(6) might be 
%[2 4 5].
%     randperms(2,6) returns might be [2 4 5 6 3].
%     randperms(2,6,2) returns 2nd element of randperms(2,6) might be 4.
%     randperms(2,6,2:4) returns 2nd-4th elements of randperms(2,6) might 
% be [4 5 6 3].
%--------------------------------------------------------------------------

switch nargin
    case 0
        error('Not enough input arguments.')
    case 1
        n = varargin{1};
        [~,p] = sort(rand(1,n)); % MR: ~ instead of ignore
    case 2
        if varargin{1} >= varargin{2}
        n = varargin{1};
        idx = varargin{2};
        [~,p] = sort(rand(1,n)); % MR: ~ instead of ignore
        p = p(idx);
        else
            n = varargin{2}-varargin{1}+1;
            [~,p] = sort(rand(1,n)); % MR: ~ instead of ignore
            p = p + (varargin{1}-1)*ones(size(p));
        end
    case 3
        if varargin{1} >= varargin{2}
            error('randperms(m,n,idx) requires positive range (n-m) > 0.')
        else
            n = varargin{2}-varargin{1}+1;
            [~,p] = sort(rand(1,n)); % MR: ~ instead of ignore
            p = p + (varargin{1}-1)*ones(size(p));
            p = p(varargin{3});
        end
    otherwise
        error('Too many input arguments')
end

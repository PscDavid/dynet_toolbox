function dosc = damposc(freqs,t,A,B,tau)
% tau             = (obj.order*obj.dt)*.25;
% t                 = (1:obj.order)*obj.dt; 
f                 = 2*pi*freqs; 
if numel(f)==1 || isrow(f)    
    A             = A*B; % in percentage scaling
    if numel(f)==1
        osc  = cos(f.*t');
    else
        osc  = sum(cos(f.*t'),2);
        osc  = osc./max(osc);
    end
    dosc = A+(((B).*exp(-t./tau)') .* osc);
else
    dosc = zeros(numel(f(:,1)),numel(t));
    for k = 1:numel(f(:,1))
        if numel(B)==1
            A_i         = A*B; % in percentage scaling
            B_i         = B;
        else
            A_i         = A*B(k);
            B_i         = B(k);
        end
        if size(f,2)==1
            osc  = cos(f(k).*t');
        else
            fline= 2*pi*(freqs(k,1):freqs(k,2));
            osc  = sum(cos(fline.*t'),2);
            osc  = osc./max(osc);
        end
        dosc(k,:) = A_i+(((B_i).*exp(-t./tau)') .* osc);
    end
end
% dosc = dosc.*sign(dosc(1));



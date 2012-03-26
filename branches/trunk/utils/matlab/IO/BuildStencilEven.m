function fx=BuildStencilEven(alpha,der)
%
% A function to compute finite-difference coefficients for the der^th 
% derivative of a function on an evenly spaced grid.  2 alpha+1 sets of 
% coefficients are returned where set 1,2,...,alpha are one-sided schemes 
% at the left end, alpha+1 is the centered scheme, and alpha+2,...,rank 
% are the one-sided schemes at the right end.  
% 
rank=2*alpha+1;

% One-sided schemes for the left end-points.  
for ip=1:alpha
    for m=-ip+1:rank-ip
        for n=1:rank
            mat(m+ip,n)=(m)^(n-1)/factorial(n-1);
        end
    end 
    mat;
    minv=inv(mat);
    fx(1:rank,ip)=minv(der+1,1:rank)';
end
% The centered scheme
for m=-alpha:alpha
    for n=1:rank
        mat(m+alpha+1,n)=(m)^(n-1)/factorial(n-1);
    end
end 
minv=inv(mat);
fx(1:rank,alpha+1)=minv(der+1,1:rank)';

% Reflect the one-sided schemes from the left end.
if mod(der,2)==0
    for ip=1:alpha
        fx(1:rank,rank-ip+1)=flipud(fx(1:rank,ip));
    end
else
    for ip=1:alpha
        fx(1:rank,rank-ip+1)=-flipud(fx(1:rank,ip));
    end
end
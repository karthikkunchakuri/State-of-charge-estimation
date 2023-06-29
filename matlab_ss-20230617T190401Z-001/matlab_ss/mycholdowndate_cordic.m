function [L] = mycholdowndate_cordic(L, x)
n = length(x);
L1 = zeros(n);
x1 = zeros(n,1);
for k = 1:n
    L1(k,k) = (L(k,k));
    x1(k) = (x(k));
%     if(abs(L1(k,k))<abs(x1(k)))
%         break;
%     end
    r = abs(sqrt(L1(k,k)^2 - x1(k)^2));
    theta = atanh(x1(k)/L1(k,k));
    if(imag(theta)~=0)
        theta=abs(theta);
    end
    
    L1(k, k) = r;
    
    if k < n
       L1((k+1):n, k) = (L((k+1):n, k) *cosh(theta) - sinh(theta) * x((k+1):n));
       x1((k+1):n) = cosh(theta) * x((k+1):n) - sinh(theta) * L((k+1):n, k);
    end
    
    L(:,k) = (L1(:,k));
    x(:,1) = x1(:,1);
end
L=L';
end

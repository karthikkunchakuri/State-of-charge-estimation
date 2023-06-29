function [L] = mycholupdate_cordic(L, x) % x=[]

n = length(x);  
L1 = zeros(n);
x1 = zeros(n,1);

for k = 1:n
    L1(k,k) = L(k,k);
    x1(k) = x(k);
    r = sqrt(L1(k,k)^2 + x1(k)^2);
    theta = atan(x1(k)/L1(k,k));
    L1(k, k) = r;
    if k < n
       L1((k+1):n, k) = (L((k+1):n, k) *cos(theta) + sin(theta) * x((k+1):n));
       x1((k+1):n) = cos(theta) * x((k+1):n) - sin(theta) * L((k+1):n, k);
    end
    L(:,k) = (L1(:,k));
    x(:,1) = x1(:,1);      
end
L=L';
end

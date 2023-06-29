function [L,p] = mycholupdowndate_cordic(L, x, sign)
p=0;
if(sign=='+')
    L=mycholupdate_cordic(L,x);
elseif(sign=='-')
    L=mycholdowndate_cordic(L,x);
else
   error("valid sign is either '+' or '-'");
end
end
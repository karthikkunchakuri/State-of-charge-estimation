function yhat = cellOutputEqn(celldz,I,Z,dR,offset,T,model)
yhat = OCVfromSOCtemp(Z+celldz,T,model); % OCV part
yhat = yhat - I * dR; % delta resistance part
yhat = yhat + offset; % polarization, hysteresis, bar resistance
% note: sensor noise handled separately
end
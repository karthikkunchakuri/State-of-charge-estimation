function [zkbar,zkbarbnd,dzk,dzkbnd,ekfData] = iterEKF_ESC(vk,ik,Tk,deltat,ekfData)
model = ekfData.model;
% Load the cell model parameters

G = getParamESC('GParam',Tk,model);
M = getParamESC('MParam',Tk,model);
M0 = getParamESC('M0Param',Tk,model);
Q  = getParamESC('QParam',Tk,model);
RC = exp(-deltat./abs(getParamESC('RCParam',Tk,model)))';
RC1 = RC(1);
RC2 = RC(2);
R  = getParamESC('RParam',Tk,model)';
R1 = R(1);
R2 = R(2);
R0 = getParamESC('R0Param',Tk,model);
eta = getParamESC('etaParam',Tk,model);
if ik<0
    ik=ik*eta;
end

% Get data stored in ekfData structure
I      = ekfData.priorI;
SigmaX = ekfData.SigmaX;
SigmaV = ekfData.SigmaV;
SigmaW = ekfData.SigmaW;
xhat   = ekfData.xhat;
irInd1 = ekfData.irInd1;
irInd2 = ekfData.irInd2;
hkInd  = ekfData.hkInd;
zkInd  = ekfData.zkInd;
signIk = ekfData.signIk;


% EKF Step 0: Compute Ahat[k-1], Bhat[k-1]
nx = length(xhat); Ahat = zeros(nx,nx); Bhat = zeros(nx,1);
Ahat(zkInd,zkInd) = 1; Bhat(zkInd) = -deltat/(3600*Q);
Ahat(irInd1,irInd1) = diag(RC1); Bhat(irInd1) = 1-RC1(:);
Ahat(irInd2,irInd2) = diag(RC2); Bhat(irInd2) = 1-RC2(:);
Ah  = exp(-abs(I*G*deltat/(3600*Q)));                                      % hysteresis factor
Ahat(hkInd,hkInd) = Ah;
B = [Bhat, 0*Bhat];
Bhat(hkInd) = -abs(G*deltat/(3600*Q))*Ah*(1+sign(I)*xhat(hkInd));
B(hkInd,2) = Ah-1;

% Step 1a: State estimate time update
xhat = Ahat*xhat + B*I;

% Step 1b: Error covariance time update
%          sigmaminus(k) = Ahat(k-1)*sigmaplus(k-1)*Ahat(k-1)' + ...
%                          Bhat(k-1)*sigmawtilde*Bhat(k-1)'
SigmaX = Ahat*SigmaX*Ahat' + Bhat*SigmaW*Bhat';

% Step 1c: Output estimate
yhat = OCVfromSOCtemp(xhat(zkInd),Tk,model) + M0*signIk + M*xhat(hkInd) ...
                               - R1*xhat(irInd1) - R2*xhat(irInd2) - R0*ik;

% Step 2a: Estimator gain matrix
Chat = zeros(1,nx);
Chat(zkInd)  = dOCVfromSOCtemp(xhat(zkInd),Tk,model);
Chat(hkInd)  = M;
Chat(irInd1) = -R1;
Chat(irInd2) = -R2;
Dhat = 1;
SigmaY = Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat';
L = SigmaX*Chat'/SigmaY;

% Step 2b: State estimate measurement update
r = mean(vk) - yhat; % residual.  Use to check for sensor errors...
if r^2 > 100*SigmaY, L(:)=0.0; end
xhat = xhat + L*r;
xhat(zkInd) = min(1.05,max(-0.05,xhat(zkInd)));

% Step 2c: Error covariance measurement update
SigmaX = SigmaX - L*SigmaY*L';
%   % Q-bump code
if r^2 > 4*SigmaY % bad voltage estimate by 2 std. devs, bump Q
    fprintf('Bumping SigmaX\n');
    SigmaX(zkInd,zkInd) = SigmaX(zkInd,zkInd)*ekfData.Qbump;
end
[~,S,V] = svd(SigmaX);
HH = V*S*V';
SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness

% Save data in ekfData structure for next time...
ekfData.priorI = ik;
ekfData.SigmaX = SigmaX;
ekfData.xhat = xhat;
zkbar = xhat(zkInd);
zkbarbnd = 3*sqrt(SigmaX(zkInd,zkInd));

% The "bar" filter update is complete. Now, work on "delta" filter updates
offset       = M*xhat(hkInd) - R1*xhat(irInd1) - R2*xhat(irInd2) - R0*(ik);
celldz       = ekfData.celldz;
cellSdz      = ekfData.cellSdz;
for thecell  = 1:length(vk)
    % Implement SPKF for delta-soc
    I = ekfData.priorI;
    % First, update the delta-SOC SPKFs
    % Step 1a - State prediction time update
    celldz(thecell) = Ahat(zkInd,zkInd)*celldz(thecell) + B(zkInd,1)*I;
    % Step 1b - Error covariance time update
    cellSdz(thecell) = Ahat(zkInd,zkInd)*cellSdz(thecell)*Ahat(zkInd,zkInd)' ...
                       + Bhat(zkInd,1)*SigmaW*Bhat(zkInd,1)';
    % Step 1c - output estimate
    I = ik;
    cellyhat = cellOutputEqn(celldz,I,zkbar,ekfData.dR0(thecell),offset,Tk,model);
    % Step 2a - Estimator gain matrix
    chatY     = dOCVfromSOCtemp(celldz(thecell)+zkbar,Tk,model);dhat = 1;
    cellY     = chatY*cellSdz(thecell)*chatY' + dhat*SigmaV*dhat';
    cellL     = cellSdz(thecell)*chatY'/cellY;
    % Step 2b - State estimate measurement update
    celldz(thecell) = celldz(thecell) + cellL*(vk(thecell) - cellyhat(thecell));
    % Step 2c - Error covariance measurement update
    cellSdz(thecell) = cellSdz(thecell) - cellL*cellY*cellL;
end

% Save data in spkfData structure for next time...
ekfData.celldz = celldz;
ekfData.cellSdz = cellSdz;
dzk = celldz;
dzkbnd = 3*sqrt(cellSdz);

end
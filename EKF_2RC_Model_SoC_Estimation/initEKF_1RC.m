function ekfData = initEKF_1RC(v0,T0,SigmaX0,SigmaV,SigmaW,model)

  % Initial state description
  ir1   = 0;                           ekfData.irInd1 = 1;
  SOC0  = mean(SOCfromOCVtemp(v0,T0,model)); ekfData.zkInd = 2;
  ekfData.xhat  = [ir1 SOC0]'; % initial state

  % Covariance values
  ekfData.SigmaX = SigmaX0;
  ekfData.SigmaV = SigmaV;
  ekfData.SigmaW = SigmaW;
  ekfData.Qbump = 5;
  
  % previous value of current
  ekfData.priorI = 0;
  ekfData.signIk = 0;
  
  % store model data structure too
  ekfData.model = model;

  % Now, initialize variables for the "delta" filters
  ekfData.celldz = 0*v0(:);
  ekfData.cellSdz = SigmaX0(ekfData.zkInd,ekfData.zkInd)*ones(size(v0(:)));  
end
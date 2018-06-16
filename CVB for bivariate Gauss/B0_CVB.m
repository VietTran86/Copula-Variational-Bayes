%//////////////////////////////////////////////////////////////// GAUSS (Ground truth)
GAUSS.sigmaX = 2;
GAUSS.sigmaY = 1;
GAUSS.rho = 0.8;

SIGMA.GAUSS = sigma2SIGMA(GAUSS.sigmaX,GAUSS.sigmaY,GAUSS.rho);
SIGMA.GAUSS0= sigma2SIGMA(GAUSS.sigmaX,GAUSS.sigmaY,0);

%//////////////////////////////////////////////////////////////// CVB approx
nVB  = 100; %maximum number of iVB

Rho  = -0.99:0.01:0.99;

nRho = length(Rho);
%------------------------------------- declaration
CVB.sigmaX = zeros(nRho,nVB);
CVB.sigmaY = zeros(nRho,nVB);
CVB.rho    = zeros(nRho,nVB);
CVB.KLD   = zeros(nRho,nVB);

MR_KLD0 = zeros(nRho,1);
CopKLD0 = zeros(nRho,1);
MR_KLDend = zeros(nRho,1);
CopKLDend = zeros(nRho,1);

nCVB = zeros(1,nRho);
KLDend = zeros(1,nRho);
%------------------------------------- initialization
CVB.sigmaX(:,1) = 1;
CVB.sigmaY(:,1) = 1;
CVB.rho(:,1) = Rho';
%//////////////////////////////////////////

SIGMA.MR = sigma2SIGMA(CVB.sigmaX(1,1),CVB.sigmaY(1,1),0);
MR_KLD0(:)  = Func_KLDMultiGauss(SIGMA.MR  ,SIGMA.GAUSS0);

for iRho = 1:nRho
    
    SIGMA.CVB = sigma2SIGMA(CVB.sigmaX(1,1),CVB.sigmaY(1,1),CVB.rho(iRho,1));
    
    CVB.KLD(iRho,1)= Func_KLDMultiGauss(SIGMA.CVB ,SIGMA.GAUSS);
    CopKLD0(iRho) = CVB.KLD(iRho,1) - MR_KLD0(1);
        
for iVB = 2:2:(nVB-2) 
    
    %------------------------------------- X to Y
    
    [CVB.sigmaX(iRho,iVB)  ,CVB.sigmaY(iRho,iVB)  ,CVB.rho(iRho,iVB),CVB.KLD(iRho,iVB)] = Func_CVB_biGauss(CVB.sigmaX(iRho,iVB-1),CVB.sigmaY(iRho,iVB-1),CVB.rho(iRho,iVB-1),GAUSS.sigmaX,GAUSS.sigmaY,GAUSS.rho);

    %------------------------------------- Y to X
    
    [CVB.sigmaY(iRho,iVB+1),CVB.sigmaX(iRho,iVB+1),CVB.rho(iRho,iVB+1),CVB.KLD(iRho,iVB+1)] = Func_CVB_biGauss(CVB.sigmaY(iRho,iVB),CVB.sigmaX(iRho,iVB),CVB.rho(iRho,iVB),GAUSS.sigmaY,GAUSS.sigmaX,GAUSS.rho);
    
    if (CVB.KLD(iRho,iVB)-CVB.KLD(iRho,iVB+1)) <= (0.01 * CVB.KLD(iRho,iVB))
        nCVB(iRho) = iVB;
        KLDend(iRho) = CVB.KLD(iRho,iVB+1);
        
        SIGMA.CVB = sigma2SIGMA(CVB.sigmaX(iRho,iVB+1),CVB.sigmaY(iRho,iVB+1),CVB.rho(iRho,iVB+1));
        SIGMA.MR  = sigma2SIGMA(CVB.sigmaX(iRho,iVB+1),CVB.sigmaY(iRho,iVB+1),0);
        
        MR_KLDend(iRho) = Func_KLDMultiGauss(SIGMA.MR  ,SIGMA.GAUSS0);
        CopKLDend(iRho) = KLDend(iRho) - MR_KLDend(iRho);
        
        break;
    end
end
end
%////////////////////////////////////////////////////////////////


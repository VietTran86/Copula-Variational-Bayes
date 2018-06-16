%///////////////////////////////////////////////////////// GAUSS

SIGMA.GAUSS = sigma2SIGMA(GAUSS.sigmaX,GAUSS.sigmaY,GAUSS.rho);

%///////////////////////////////////////////////////////// VB

rhoVB = find(Rho == 0);
numVB = 1;

SIGMA.VB0 = sigma2SIGMA(CVB.sigmaX(rhoVB,numVB),CVB.sigmaY(rhoVB,numVB),CVB.rho(rhoVB,numVB));

%-------------------------------------
numVB = nCVB(rhoVB);

SIGMA.VB = sigma2SIGMA(CVB.sigmaX(rhoVB,numVB),CVB.sigmaY(rhoVB,numVB),CVB.rho(rhoVB,numVB));

nuVB = numVB;
KLDVB = CVB.KLD(rhoVB,1:nuVB+1);

%///////////////////////////////////////////////////////// CVB
rhos = [0.5,-0.5]; % "Fig. 11" in https://arxiv.org/abs/1803.10998 (v1)

for irho = 1:2
    
    plot_rho =  rhos(irho);
    
    rhoCVB = find(abs(Rho-plot_rho) < 1e-3);
    numCVB = 1;
    
    SIGMA.CVB0 = sigma2SIGMA(CVB.sigmaX(rhoCVB,numCVB),CVB.sigmaY(rhoCVB,numCVB),CVB.rho(rhoCVB,numCVB));
    
    %-------------------------------------
    
    rhoCVB = find(abs(Rho-plot_rho) < 1e-3);
    numCVB = nCVB(rhoCVB);
    
    SIGMA.CVB = sigma2SIGMA(CVB.sigmaX(rhoCVB,numCVB),CVB.sigmaY(rhoCVB,numCVB),CVB.rho(rhoCVB,numCVB));
    
    nuCVB = numCVB;
    KLDCVB = CVB.KLD(rhoCVB,1:nuCVB+1);
    
    %-------------------------------------
    %////////////////////////////////////////////////////////////////
    
    MU = [0,0];
    
    %-------------------------------------
    [x1,y1,F.GAUSS] = Func_contourGauss(MU,SIGMA.GAUSS);
    [x2 ,y2 ,F.VB]   = Func_contourGauss(MU,SIGMA.VB);
    [x3a,y3a,F.CVB0] = Func_contourGauss(MU,SIGMA.CVB0);
    [x3 ,y3 ,F.CVB]  = Func_contourGauss(MU,SIGMA.CVB);
    
    v1 = [.003 .01 .1]; %levels
    v2 = [.003 .01 .1]; %levels
    v3 = [.003 .01 .1]; %levels
    
    h1 = figure(irho);
    hold on
    
    xoffset = -8;
    
    plot(x1,10*normpdf(x1,0,sqrt(SIGMA.GAUSS(1,1)))+ xoffset,'k','LineWidth',0.7);
    plot(x2,10*normpdf(x2,0,sqrt(SIGMA.VB(1,1)))  + xoffset,'b-.');
    plot(x3,10*normpdf(x3,0,sqrt(SIGMA.CVB0(1,1)))+ xoffset,'r:');
    plot(x3,10*normpdf(x3,0,sqrt(SIGMA.CVB(1,1))) + xoffset,'r--','LineWidth',0.7);
    
    yoffset = -15;
    
    plot(10*normpdf(y1,0,sqrt(SIGMA.GAUSS(2,2))) + yoffset,y1,'k','LineWidth',0.7);
    plot(10*normpdf(y2,0,sqrt(SIGMA.VB(2,2)))   + yoffset,y2,'b-.');
    plot(10*normpdf(y3,0,sqrt(SIGMA.CVB0(2,2))) + yoffset,y3,'r:');
    plot(10*normpdf(y3,0,sqrt(SIGMA.CVB(2,2)))  + yoffset,y3,'r--','LineWidth',0.7);
    
    contour(x1 ,y1 ,F.GAUSS,v1,'k','LineWidth',0.7);
    contour(x2 ,y2 ,F.VB ,v2,'b-.');
    contour(x3a,y3a,F.CVB0,v3,'r:');
    contour(x3 ,y3 ,F.CVB ,v3,'r--','LineWidth',0.7);
    
    set(h1, 'defaultAxesTickLabelInterpreter','latex');
    set(h1, 'defaultLegendInterpreter','latex');
    legend('Gauss','VB (converged)',['CVB (inital), $\tilde{\rho_{[0]}}$ = ',num2str(plot_rho)],'CVB (converged)')
    
    ylabel('\theta_2')
    xlabel('\theta_1')
    
    if irho == 1
        movegui(h1,'northwest')
        
    else
        movegui(h1,'north')
    end
end
%////////////////////////////////////////////////////////////////



h2 = figure(3);
hold on
plot(0,CVB.KLD(rhoVB,2),'kx');
plot(0,KLDend(rhoVB),'ko');
plot(Rho,CVB.KLD(:,1),'k--');
plot(Rho,KLDend,'k');

legend('VB (inital)','VB (converged)','CVB (inital)','CVB (converged)')
ylabel('KL divergence')
xlabel('$\tilde{\rho}_{[0]}$','Interpreter','latex')

ax = gca;
ax.YAxisLocation = 'origin'; 

movegui(h2,'northeast')
%-------------------------------------
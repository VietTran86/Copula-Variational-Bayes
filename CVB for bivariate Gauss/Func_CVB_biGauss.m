function [X1,Y1,rho1,KLD] = Func_CVB_biGauss(X0,Y0,rho0,sigmaX,sigmaY,rho)
%//////////////////////////////////////////////////////////////////////
%--------------------------------- TRUE values
betaYX  = rho * sigmaY/sigmaX;
sigmaYX = sigmaY * sqrt(1-rho^2);

%--------------------------------- CVB
bYX = rho0 * Y0/X0;
YX = Y0 * sqrt(1-rho0^2);
%---------------------------------
X1 = 1/sqrt(1/sigmaX^2 + (bYX - betaYX)^2/sigmaY^2/(1-rho^2));

KLD = -log(X1/sigmaX * YX/sigmaYX * exp((sigmaYX^2-YX^2)/2/sigmaYX^2));
%---------------------------------
rho1 = rho0/sqrt(rho0^2 +(X0/X1)^2 * (1-rho0^2));

Y1 = Y0 * sqrt(rho0^2 * (X1/X0)^2 + (1-rho0^2));
%---------------------------------
%//////////////////////////////////////////////////////////////////////
end
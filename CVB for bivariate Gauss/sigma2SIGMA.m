function SIGMA = sigma2SIGMA(sigmaX,sigmaY,rho)

rX = sigmaX^2;
rY = sigmaY^2;
rXY= rho * sigmaY * sigmaX;

SIGMA = [rX rXY; rXY rY];

end
function [x1,x2,F] = Func_contourGauss(MU,SIGMA)
% openExample('stats/ComputeCumulativeProbabilitiesOverRegionsExample')
% https://www.mathworks.com/help/stats/multivariate-normal-distribution.html

sigma1 = sqrt(SIGMA(1,1));
sigma2 = sqrt(SIGMA(2,2));

x1 = -5*sigma1:0.05:5*sigma1; 
x2 = -5*sigma2:0.05:5*sigma2; 
[X1,X2] = meshgrid(x1,x2);

F = mvnpdf([X1(:) X2(:)],MU,SIGMA);
F = reshape(F,length(x2),length(x1));

end
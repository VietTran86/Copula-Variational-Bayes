function KLD = Func_KLDMultiGauss(SIGMA0,SIGMA1)
N = size(SIGMA1,1);

KLD = 1/2*(trace(SIGMA1\SIGMA0)+log(det(SIGMA1)/det(SIGMA0))-N);
%http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Kullback.E2.%80.93Leibler_divergence

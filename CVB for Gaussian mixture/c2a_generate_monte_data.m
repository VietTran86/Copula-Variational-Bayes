%//////////////////////////////////////////////////////////////////////////

X_gauss = zeros(setting.M,setting.N,setting.K);

for k=1:setting.K
    X_gauss(:,:,k) = mvnrnd(para.MU(:,k),[1 1],setting.N)';
end

L_multinomial = mnrnd(1,ones(1,setting.K)/setting.K,setting.N)';

%-------------------------------------

data.L = (1:setting.K) * L_multinomial;

data.X = zeros(setting.M,setting.N);

for i=1:setting.N
    data.X(:,i) = X_gauss(:,i,data.L(i));
end

%//////////////////////////////////////////////////////////////////////////
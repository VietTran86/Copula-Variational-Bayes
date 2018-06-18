function [iPurity,iMSE,iEntropy] = Func_ClusterEvaluation(estL,estMU,trueL,trueMU)


in.K = size(estMU,2); % number of initial clusters
in.N = size(estL,2);  % number of data (time points)

%------------------------------
jCorrect = zeros(1,in.K);
iCorrect = zeros(1,in.K);

for i=1:in.K
    for j=1:in.K
        jCorrect(j) = sum(estL(trueL == j) == i);
    end
    iCorrect(i) = max(jCorrect);
end
iPurity = sum(iCorrect)/in.N * 100; % percent

%------------------------------

Eij = zeros(1,in.K);
Pij = zeros(1,in.K);
Ei  = zeros(1,in.K);
Pi  = zeros(1,in.K);

for ik=1:in.K
    Ei(ik) = sum(estL == ik);
    if Ei(ik) ~= 0
        for j=1:in.K
            Eij(j) = sum(estL(trueL == j) == ik);
            Pij(j) = Eij(j)/Ei(ik);
        end
        Pi(ik) = -sum(Pij(Pij ~=0) .* log(Pij(Pij ~=0)));
    else
        Pi(ik) = 0;
    end
end
iEntropy = sum(Pi.*Ei)/in.N;

%------------------------------

numPerms = factorial(in.K);
kMSE = zeros(1,numPerms);

inxPerms = perms(1:in.K);

for jp=1:numPerms
    colPerm = inxPerms(jp,:);
    kMSE(jp) = mean(Euclidean(estMU(:,colPerm),trueMU));
end

iMSE = min(kMSE);


%//////////////////////////////////////////////////////////////////////////
return

    function EU = Euclidean(MU1,MU2)
        EU = sum((MU1 - MU2).^2);
    end

%//////////////////////////////////////////////////////////////////////
end
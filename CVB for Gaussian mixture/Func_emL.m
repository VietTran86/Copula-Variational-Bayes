function [emL_L,emL_MU,emL_ELBO,emL_numLoop] = Func_emL(init_pos,dataX,maxLoop)
%//////////////////////////////////////////////////////////////////////////

in.M = size(init_pos,1); % number of data's dimension
in.K = size(init_pos,2); % number of initial clusters
in.N = size(dataX,2);     % number of data (time points)

%------------------------------ declaration

emL.MU   = zeros(in.M,in.K,maxLoop);
emL.sigma =zeros(  in.K,maxLoop);
emL.L   = zeros(maxLoop,in.N);
emL.P   = zeros(in.K,in.N);
logP    = zeros(in.K,in.N);
ELBO    = zeros(1,maxLoop);
%------------------------------ initialization

emL.MU(:,:,1) = init_pos;

emL.sigma(:,1) = 1;
%//////////////////////////////////////////////////////////////////////////

Flag.stopLoop = false;
iLoop = 1;

while (~Flag.stopLoop)
    %------------------------------
    iLoop = iLoop + 1;
    %------------------------------
    
    for k = 1:in.K
        
        logP(k,:) = logNormalPDF(dataX,emL.MU(:,k,iLoop-1)) - emL.sigma(k,iLoop-1)^2;
    end
    
    [~,emL.L(iLoop,:)] = max(logP);
    
    for k = 1:in.K
        Lk = emL.L(iLoop,:) == k;
        
        if any(Lk)
            emL.sigma(  k,iLoop) = 1/sqrt(sum(Lk));
            emL.MU   (:,k,iLoop) = sum(dataX(:,Lk),2)/sum(Lk);
        else
            emL.sigma(  k,iLoop) = emL.sigma(  k,iLoop-1);
            emL.MU   (:,k,iLoop) = emL.MU   (:,k,iLoop-1);
        end
    end
    
    %------------------------------
    
    ELBO(iLoop) = ClusterELBO(emL.L(iLoop,:),emL.MU(:,:,iLoop),emL.sigma(:,iLoop));
    
    %------------------------------
    if (iLoop == 2)
        
    elseif (iLoop == maxLoop) || all(emL.L(iLoop,:) == emL.L(iLoop-1,:))
        
        Flag.stopLoop = true;
        
        results()
    end
    
    %------------------------------
    
end

%//////////////////////////////////////////////////////////////////////////
return

%//////////////////////////////////////////////////////////////////////

    function results()
        
        emL_L  = emL.L(iLoop,:);
        emL_MU = emL.MU(:,:,iLoop);
        emL_numLoop = iLoop;
        emL_ELBO = ELBO(iLoop);
    end

%//////////////////////////////////////////////////////////////////////

    function ELBO = ClusterELBO(L,MU,sigma)
        
        iELBO = sum(logNormalPDF(dataX,MU(:,L)));
        
        logKLDGauss = log(prod(sigma)^2) + in.K * (log(2*pi) + 1);
        
        Denom = zeros(1,in.K);
        for jk=1:in.K
            Denom(jk) = sigma(jk)^2 * sum(L==jk);
        end
        
        ELBO = iELBO + logKLDGauss - sum(Denom);
    end

    %//////////////////////////////////////////////////////////////////////

    function LOGmvnpdf = logNormalPDF(X,MU)
        LOGmvnpdf = -log(2*pi) - Euclidean(X,MU)/2;
    end

    function EU = Euclidean(MU1,MU2)
        EU = sum((MU1 - MU2).^2);
    end

%//////////////////////////////////////////////////////////////////////
end

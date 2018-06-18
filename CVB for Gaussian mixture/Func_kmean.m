function [kmean_L,kmean_MU,kmean_ELBO,kmean_numLoop] = Func_kmean(init_pos,dataX,maxLoop)
%//////////////////////////////////////////////////////////////////////////        

in.M = size(init_pos,1); % number of data's dimension
in.K = size(init_pos,2); % number of initial clusters
in.N = size(dataX,2);    % number of data (time points)

%------------------------------ declaration

distance = zeros(in.K,in.N);
kmean.L  = zeros(maxLoop,in.N);   % L  is "labels of clusters"
kmean.MU = zeros(in.M,in.K,maxLoop); % MU is "mean vectors" 
kmean.ELBO = zeros(1,maxLoop);

%------------------------------ initialization

kmean.MU(:,:,1) = init_pos; 

%//////////////////////////////////////////////////////////////////////////

Flag.stopLoop = false;
iLoop  = 1;

while (~Flag.stopLoop)
    %------------------------------
    iLoop = iLoop + 1;
    %------------------------------
    
    for k = 1:in.K
        distance(k,:) = Euclidean(dataX,kmean.MU(:,k,iLoop-1));
    end
    
    [~,kmean.L(iLoop,:)] = min(distance);
    
    for k = 1:in.K
        Lk = (kmean.L(iLoop,:) == k);
        
        if any(Lk)
            kmean.MU(:,k,iLoop) = sum(dataX(:,Lk),2)/sum(Lk);
        else
            kmean.MU(:,k,iLoop) = kmean.MU(:,k,iLoop-1);
        end
    end
    
    %------------------------------ 
    
    kmean.ELBO(iLoop) = ClusterELBO(kmean.L(iLoop,:),kmean.MU(:,:,iLoop));
    
    %------------------------------
    if (iLoop == 2)
        
    elseif (iLoop == maxLoop) || all(kmean.L(iLoop,:) == kmean.L(iLoop-1,:))
        
        results();
        
        Flag.stopLoop = true;
    end

    %------------------------------
end

%//////////////////////////////////////////////////////////////////////////
return

    %//////////////////////////////////////////////////////////////////////
    
    function results() 

        kmean_numLoop = iLoop;
        kmean_L    = kmean.L(iLoop,:);
        kmean_MU   = kmean.MU(:,:,iLoop);
        kmean_ELBO = kmean.ELBO(iLoop);
        
    end

    %//////////////////////////////////////////////////////////////////////

    function iELBO = ClusterELBO(L,MU)
        
        iMU = MU(:,L);
        
        iELBO = sum(logNormalPDF(dataX,iMU));
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

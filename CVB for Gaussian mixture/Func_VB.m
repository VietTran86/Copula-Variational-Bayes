function [VB_L,VB_MU,VB_ELBO,VB_numLoop] = Func_VB(init_pos,dataX,maxLoop,ELBOthresh)
%//////////////////////////////////////////////////////////////////////////        

in.M = size(init_pos,1); % number of data's dimension
in.K = size(init_pos,2); % number of initial clusters
in.N = size(dataX,2);     % number of data (time points)

%------------------------------ declaration

VB.MU    = zeros(in.M,in.K,maxLoop);
VB.sigma = zeros(  in.K,maxLoop);
VB.L     = zeros(   1,in.N);
VB.P     = zeros(in.K,in.N);
logP     = zeros(in.K,in.N);
ELBO     = zeros(1,maxLoop);
%------------------------------ initialization

VB.MU(:,:,1) = init_pos;

VB.sigma(:,1) = 1;
%//////////////////////////////////////////////////////////////////////////

Flag.stopLoop = false;
iLoop  = 1;


while (~Flag.stopLoop)
    %------------------------------
    iLoop = iLoop + 1;
    %------------------------------
    
    for k = 1:in.K
        
        logP(k,:) = logNormalPDF(dataX,VB.MU(:,k,iLoop-1)) - VB.sigma(k,iLoop-1)^2;
        
    end
    
    VB.P = exp(logP - repmat(logSumExp(logP,1),in.K,1));
    
    for k=1:in.K
        
        sumPk = sum(VB.P(k,:));
        
        if sumPk == 0
            VB.MU(:,k,iLoop) = VB.MU(:,k,iLoop-1);
            
            VB.sigma(k,iLoop) = VB.sigma(k,iLoop-1);
        else
            VB.MU(:,k,iLoop) = dataX * (VB.P(k,:)' / sumPk);
            
            VB.sigma(k,iLoop) = 1/sqrt(sumPk);
        end
    end
    
    [~,VB.L] = max(VB.P);  
    %------------------------------ 
        
     ELBO(iLoop) = ClusterELBO(VB.P,VB.MU(:,:,iLoop),VB.sigma(:,iLoop));
    
    %------------------------------
    
    ELBOincrease = ELBO(iLoop) - ELBO(iLoop-1);
    
    if (iLoop == 2)
        
    elseif (iLoop == maxLoop) || (ELBOincrease >= 0 && ELBOincrease <= ELBOthresh)
        
        Flag.stopLoop = true;
        
        results()
    end

end


%//////////////////////////////////////////////////////////////////////////
return

    %//////////////////////////////////////////////////////////////////////
    
    function results()    
        
        VB_L  = VB.L;
        VB_MU = VB.MU(:,:,iLoop);
        VB_ELBO    = ELBO(iLoop);
        VB_numLoop = iLoop;
    end

    %//////////////////////////////////////////////////////////////////////

    function ELBO = ClusterELBO(P,MU,sigma)
                
        kELBO = zeros(1,in.K);
        for j=1:in.K
            kELBO(j) = log(2*pi) + 2*log(sigma(j)) + sum(P(j,:).*logNormalPDF(dataX,MU(:,j)));            
        end
        ELBO = sum(kELBO) - sum(sum(log(P.^P)));
    end
    
    %//////////////////////////////////////////////////////////////////////

    function LOGmvnpdf = logNormalPDF(X,MU)
        LOGmvnpdf = -log(2*pi) - Euclidean(X,MU)/2;
    end

    function EU = Euclidean(MU1,MU2)
        EU = sum((MU1 - MU2).^2); 
    end

    %//////////////////////////////////////////////////////////////////////

    function [LOGX,logMaxX] = logSumExp(logX,dim)
        
        if dim > 0
        
            logMaxX = max(logX,[],dim);
        
            LOGX = log(sum(exp(logX - logMaxX),dim)) + logMaxX;
            
        else
            
            logx = logX(:);
            
            logMaxX = max(logx);
           
            LOGX = log(sum(exp(logx - logMaxX))) + logMaxX;
        end
    end

    %//////////////////////////////////////////////////////////////////////
   
    
end

function [emMU_L,emMU_MU,emMU_ELBO,emMU_numLoop] = Func_emMU(init_pos,dataX,maxLoop,ELBOthresh)
%//////////////////////////////////////////////////////////////////////////        

in.M = size(init_pos,1); % number of data's dimension
in.K = size(init_pos,2); % number of initial clusters
in.N = size(dataX,2);     % number of data (time points)

%------------------------------ declaration

emMU.MU = zeros(in.M,in.K,maxLoop);
emMU.L  = zeros(   1,in.N);
emMU.P  = zeros(in.K,in.N);
  logP  = zeros(in.K,in.N);
ELBO    = zeros(1,maxLoop);
%------------------------------ initialization

emMU.MU(:,:,1) = init_pos;

%//////////////////////////////////////////////////////////////////////////

Flag.stopLoop = false;
iLoop  = 1;

while (~Flag.stopLoop)
    %------------------------------
    iLoop = iLoop + 1;
    %------------------------------
    
    for k = 1:in.K
        
        logP(k,:) = logNormalPDF(dataX,emMU.MU(:,k,iLoop-1));
    end
    
    emMU.P = exp(logP - repmat(logSumExp(logP,1),in.K,1));
        
    for k=1:in.K
        
        sumPk = sum(emMU.P(k,:));
        
        if sumPk == 0
            
            emMU.MU(:,k,iLoop) = emMU.MU(:,k,iLoop-1);
            
        else
            emMU.MU(:,k,iLoop) = dataX * (emMU.P(k,:)' / sumPk);

        end
    end
    
    %------------------------------ 
    
    ELBO(iLoop)  = ClusterELBO(emMU.P,emMU.MU(:,:,iLoop));
        
    [~,emMU.L] = max(emMU.P);
       
    %------------------------------
    
    ELBOincrease = ELBO(iLoop) - ELBO(iLoop-1);
    
    if (iLoop == 2)
        
    elseif (iLoop == maxLoop) || (ELBOincrease >= 0 && ELBOincrease <= ELBOthresh)
        
        Flag.stopLoop = true;
        
        results()
    end
    %------------------------------
end

%//////////////////////////////////////////////////////////////////////////
return

    function results()    
        
        emMU_L  = emMU.L;
        emMU_MU = emMU.MU(:,:,iLoop);
        emMU_ELBO    = ELBO(iLoop);
        emMU_numLoop = iLoop;
    end

    %//////////////////////////////////////////////////////////////////////

    function ELBO = ClusterELBO(P,MU)
        
        kELBO = zeros(1,in.K);
        for j=1:in.K
            kELBO(j) = sum(P(j,:) .* logNormalPDF(dataX,MU(:,j)));
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

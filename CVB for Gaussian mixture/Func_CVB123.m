function [CVB1_L,CVB1_MU,CVB1_ELBO,...
          CVB2_L,CVB2_MU,CVB2_ELBO,...
          CVB3_L,CVB3_MU,CVB3_ELBO,CVB_numLoop] = Func_CVB123(init_pos,dataX,maxLoop,ELBOthresh)
%//////////////////////////////////////////////////////////////////////////        

in.M = size(init_pos,1); % number of data's dimension
in.K = size(init_pos,2); % number of initial clusters
in.N = size(dataX,2);     % number of data (time points)

numPick = in.N;
%------------------------------ declaration

jMU    = zeros(in.M,in.K,numPick);
jP     = zeros(in.K,in.N,numPick);
jELBO  = zeros(1,numPick);

CVB.MU    = zeros(in.M,in.K,maxLoop);
CVB.P     = zeros(in.K,in.N);
CVB.L     = zeros(   1,in.N);
CVB.jL    = zeros(   1,numPick);
CVB.Q     = zeros(   1,numPick);
%------------------------------ 
for iPick = 1:numPick
    
    CVBj.MU    = zeros(in.M,in.K,in.K,maxLoop);
    CVBj.sigma = zeros(  in.K,in.K,maxLoop);
    CVBj.L     = zeros(in.K,in.N);
    CVBj.W     = zeros(in.K,in.K,in.N);
      logW     = zeros(in.K,in.K,in.N);
    CVBj.P     = zeros(in.K,in.N);
    CVBj.ELBO  = zeros(1,maxLoop);

CVBj.sigma(:,:,1) = 1;

for m = 1:in.K
    
    CVBj.MU(:,:,m,1) = init_pos;
end
%//////////////////////////////////////////////////////////////////////////
    
Flag.stopLoop = false;

iLoop = 1;
    
while (~Flag.stopLoop)
    %------------------------------
    iLoop   = iLoop + 1;
    %------------------------------
    
    for k = 1:in.K
    for m = 1:in.K
        
        logW(k,m,:) = logNormalPDF(dataX,CVBj.MU(:,k,m,iLoop-1)) - CVBj.sigma(k,m,iLoop-1)^2;
    end
    end
    
    %------------------------------
    CVBj.W = exp(logW - repmat(logSumExp(logW,1),in.K,1,1));

    CVBj.W(:,:,iPick) = eye(in.K);
    %------------------------------ 
    logGammak = zeros(in.K,in.K);
    
    for k=1:in.K    
    for m=1:in.K
        
        Wk = CVBj.W(k,m,:);
        Wk = Wk(:);
        
        sumWk = sum(Wk);
        
        if sumWk == 0
            CVBj.MU(:,k,m,iLoop) = CVBj.MU(:,k,m,iLoop-1);
            
            CVBj.sigma(k,m,iLoop) = CVBj.sigma(k,m,iLoop-1);
            
            logGammak(k,m) = 0;
            
        else
            
            CVBj.MU(:,k,m,iLoop) = dataX * (Wk / sumWk);
            
            CVBj.sigma(k,m,iLoop) = 1/sqrt(sumWk);
            
            %------------------------------ 
            MU    = CVBj.MU(:,k,m,iLoop);
            
            W_inotj        = Wk;
            W_inotj(iPick) = [];
            
            term_inotj = sum(log(W_inotj.^W_inotj));
            
            logGammak(k,m) = log(2*pi) - log(sumWk) + sum(Wk'.*(logNormalPDF(dataX,MU))) - term_inotj;
            
            %------------------------------ 
        end
    end
    end
    
    logGamma = sum(logGammak);
    
    logGammaSum = logSumExp(logGamma,2);
    
    CVBj.P(:,iPick) = exp(logGamma - logGammaSum)';
    
    CVBj.ELBO(iLoop) = logGammaSum;
    
    %------------------------------
        
    ELBOincrease = CVBj.ELBO(iLoop) - CVBj.ELBO(iLoop-1);
    
    if (iLoop == 2)
    
    elseif (iLoop == maxLoop) || (ELBOincrease >= 0 && ELBOincrease <= ELBOthresh)
        
        CVBj.nLoop(iPick) = iLoop;
        
        Flag.stopLoop = true;
        
        %------------------------------
        mMU = zeros(2,in.K,in.K);
        for m=1:in.K
            mMU(:,:,m) = CVBj.P(m,iPick) * CVBj.MU(:,:,m,iLoop-1);
        end
        kMU = sum(mMU,3);
        
        for i=1:in.N
            if i ~= iPick
                CVBj.P(:,i) = CVBj.W(:,:,i) * CVBj.P(:,iPick);
            end
        end
        %------------------------------
        
        jELBO(iPick)   = CVBj.ELBO(iLoop);
        jMU(:,:,iPick) = kMU;
        jP (:,:,iPick) = CVBj.P;
        
       
    end

    %------------------------------
    
end

end

nLoopmean = mean(CVBj.nLoop);

CVB_numLoop = nLoopmean;

%/////////////////////////////////////

Pick.MU   = mean(jMU,3); %METHOD1: average uniform
Pick.ELBO = mean(jELBO);

Pick.P = zeros(in.K,in.N);
for iPick=1:in.N
    Pick.P(:,iPick) = jP (:,iPick,iPick);
end
[~,Pick.L] = max(Pick.P);

%------------------------------------
CVB1_L    = Pick.L;
CVB1_MU   = Pick.MU;
CVB1_ELBO = Pick.ELBO;
%------------------------------------

%/////////////////////////////////////

CVB2.Q = zeros(1,numPick);  %METHOD2: pick maximum ELBO
[~,inxMaxELBO] = max(jELBO);
CVB2.Q(inxMaxELBO) = 1;

[CVB2_L,CVB2_MU,CVB2_ELBO] = augmentCVB(CVB2.Q,jMU,jP);

%/////////////////////////////////////

CVB3.Q = exp(jELBO - logSumExp(jELBO,0)); %METHOD3: average ELBO

[CVB3_L,CVB3_MU,CVB3_ELBO] = augmentCVB(CVB3.Q,jMU,jP);

%//////////////////////////////////////////////////////////////////////////
return
    
    function [L,qMU,ELBO] = augmentCVB(Q,jMU,jP)
        
        tamMU    = zeros(   2,in.K,numPick);
        tamP     = zeros(in.K,in.N,numPick);
        for jk = 1:numPick
            tamMU(:,:,jk) = jMU(:,:,jk) * Q(jk);
            tamP (:,:,jk) = jP (:,:,jk) * Q(jk);
        end
        
        qMU = sum(tamMU,3);
        qP  = sum(tamP,3);
        
        [~,L] = max(qP);
        
        ELBO = jELBO * Q';
        
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

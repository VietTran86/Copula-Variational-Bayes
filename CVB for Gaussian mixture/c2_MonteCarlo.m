
%//////////////////////////////////////////// delaration

kmean_ELBO    = zeros(1,setting.MonteCarlo);
kmean_numLoop = zeros(1,setting.MonteCarlo);
kmean_Purity  = zeros(1,setting.MonteCarlo);
kmean_MSE     = zeros(1,setting.MonteCarlo);
kmean_Entropy = zeros(1,setting.MonteCarlo);

emL_ELBO    = zeros(1,setting.MonteCarlo);
emL_numLoop = zeros(1,setting.MonteCarlo);
emL_Purity  = zeros(1,setting.MonteCarlo);
emL_MSE     = zeros(1,setting.MonteCarlo);
emL_Entropy = zeros(1,setting.MonteCarlo);

emMU_ELBO    = zeros(1,setting.MonteCarlo);
emMU_numLoop = zeros(1,setting.MonteCarlo);
emMU_Purity  = zeros(1,setting.MonteCarlo);
emMU_MSE     = zeros(1,setting.MonteCarlo);
emMU_Entropy = zeros(1,setting.MonteCarlo);

VB_ELBO    = zeros(1,setting.MonteCarlo);
VB_numLoop = zeros(1,setting.MonteCarlo);
VB_Purity  = zeros(1,setting.MonteCarlo);
VB_MSE     = zeros(1,setting.MonteCarlo);
VB_Entropy = zeros(1,setting.MonteCarlo);

CVB1_ELBO    = zeros(1,setting.MonteCarlo);
CVB1_Purity  = zeros(1,setting.MonteCarlo);
CVB1_MSE     = zeros(1,setting.MonteCarlo);
CVB1_Entropy = zeros(1,setting.MonteCarlo);

CVB2_ELBO    = zeros(1,setting.MonteCarlo);
CVB2_Purity  = zeros(1,setting.MonteCarlo);
CVB2_MSE     = zeros(1,setting.MonteCarlo);
CVB2_Entropy = zeros(1,setting.MonteCarlo);

CVB3_ELBO    = zeros(1,setting.MonteCarlo);
CVB3_Purity  = zeros(1,setting.MonteCarlo);
CVB3_MSE     = zeros(1,setting.MonteCarlo);
CVB3_Entropy = zeros(1,setting.MonteCarlo);

CVB_numLoop = zeros(1,setting.MonteCarlo);

%////////////////////////////////////////////

for monte = 1:(setting.MonteCarlo)
    
    Flag.monte.init = (monte == 1);
        
    c2a_generate_monte_data;
 
    %////////////////////////////////////////////
    
    [kmean_L,kmean_MU,kmean_ELBO(monte),kmean_numLoop(monte)] = Func_kmean(setting.init_pos,data.X,setting.maxLoop);
    [emL_L  ,emL_MU  ,emL_ELBO(monte)  ,emL_numLoop(monte)  ] = Func_emL  (setting.init_pos,data.X,setting.maxLoop);
    [emMU_L ,emMU_MU ,emMU_ELBO(monte) ,emMU_numLoop(monte) ] = Func_emMU (setting.init_pos,data.X,setting.maxLoop,setting.ELBOthresh);
    [VB_L   ,VB_MU   ,VB_ELBO(monte)   ,VB_numLoop(monte)   ] = Func_VB   (setting.init_pos,data.X,setting.maxLoop,setting.ELBOthresh);
    
    [CVB1_L,CVB1_MU,CVB1_ELBO,...
     CVB2_L,CVB2_MU,CVB2_ELBO,...
     CVB3_L,CVB3_MU,CVB3_ELBO,CVB_numLoop] = Func_CVB123(setting.init_pos,data.X,setting.maxLoop,setting.ELBOthresh);

    %////////////////////////////////////////////
    
    trueL  = data.L;
    trueMU = para.MU;
    
    [kmean_Purity(monte),kmean_MSE(monte),kmean_Entropy(monte)] = Func_ClusterEvaluation(kmean_L,kmean_MU,trueL,trueMU);
    [emL_Purity(monte)  ,emL_MSE(monte)  ,emL_Entropy(monte)  ] = Func_ClusterEvaluation(emL_L  ,emL_MU  ,trueL,trueMU);
    [emMU_Purity(monte) ,emMU_MSE(monte) ,emMU_Entropy(monte) ] = Func_ClusterEvaluation(emMU_L ,emMU_MU ,trueL,trueMU);
    [VB_Purity(monte)   ,VB_MSE(monte)   ,VB_Entropy(monte)   ] = Func_ClusterEvaluation(VB_L   ,VB_MU   ,trueL,trueMU);
    [CVB1_Purity(monte) ,CVB1_MSE(monte) ,CVB1_Entropy(monte) ] = Func_ClusterEvaluation(CVB1_L ,CVB1_MU ,trueL,trueMU);
    [CVB2_Purity(monte) ,CVB2_MSE(monte) ,CVB2_Entropy(monte) ] = Func_ClusterEvaluation(CVB2_L ,CVB2_MU ,trueL,trueMU);
    [CVB3_Purity(monte) ,CVB3_MSE(monte) ,CVB3_Entropy(monte) ] = Func_ClusterEvaluation(CVB3_L ,CVB3_MU ,trueL,trueMU);
    
    %//////////////////////////////////////////// plot the case of "Radius = 4"
    
    if (monte == 1) && (Radius == setting.plotRadius)
        
        c2b_clusterplot;
        
    end
    
    
end

%//////////////////////////////////////////// result

mean_kmean_ELBO(inxRun)   = mean(kmean_ELBO);
mean_kmean_numLoop(inxRun) = mean(kmean_numLoop);
mean_kmean_Purity(inxRun)  = mean(kmean_Purity);
mean_kmean_MSE(inxRun)     = mean(kmean_MSE);
mean_kmean_Entropy(inxRun) = mean(kmean_Entropy);

mean_emL_ELBO(inxRun)    = mean(emL_ELBO);
mean_emL_numLoop(inxRun) = mean(emL_numLoop);
mean_emL_Purity(inxRun)  = mean(emL_Purity);
mean_emL_MSE(inxRun)     = mean(emL_MSE);
mean_emL_Entropy(inxRun) = mean(emL_Entropy);

mean_emMU_ELBO(inxRun)    = mean(emMU_ELBO);
mean_emMU_numLoop(inxRun) = mean(emMU_numLoop);
mean_emMU_Purity(inxRun)  = mean(emMU_Purity);
mean_emMU_MSE(inxRun)     = mean(emMU_MSE);
mean_emMU_Entropy(inxRun) = mean(emMU_Entropy);

mean_VB_ELBO(inxRun)    = mean(VB_ELBO);
mean_VB_numLoop(inxRun) = mean(VB_numLoop);
mean_VB_Purity(inxRun)  = mean(VB_Purity);
mean_VB_MSE(inxRun)     = mean(VB_MSE);
mean_VB_Entropy(inxRun) = mean(VB_Entropy);

mean_CVB1_ELBO(inxRun)    = mean(CVB1_ELBO);
mean_CVB1_Purity(inxRun)  = mean(CVB1_Purity);
mean_CVB1_MSE(inxRun)     = mean(CVB1_MSE);
mean_CVB1_Entropy(inxRun) = mean(CVB1_Entropy);

mean_CVB2_ELBO(inxRun)    = mean(CVB2_ELBO);
mean_CVB2_Purity(inxRun)  = mean(CVB2_Purity);
mean_CVB2_MSE(inxRun)     = mean(CVB2_MSE);
mean_CVB2_Entropy(inxRun) = mean(CVB2_Entropy);

mean_CVB3_ELBO(inxRun)    = mean(CVB3_ELBO);
mean_CVB3_Purity(inxRun)  = mean(CVB3_Purity);
mean_CVB3_MSE(inxRun)     = mean(CVB3_MSE);
mean_CVB3_Entropy(inxRun) = mean(CVB3_Entropy);

mean_CVB_numLoop(inxRun) = mean(CVB_numLoop);

%////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////

h.MSE = figure;
hold on;
%-----------------------------------------
plot(setting.Radius,mean_kmean_MSE,spec.kmean{:});
plot(setting.Radius,mean_emL_MSE,spec.emL{:});
plot(setting.Radius,mean_emMU_MSE,spec.emMU{:});
plot(setting.Radius,mean_VB_MSE,spec.VB{:});
plot(setting.Radius,mean_CVB1_MSE,spec.CVB1{:});
plot(setting.Radius,mean_CVB2_MSE,spec.CVB2{:});
plot(setting.Radius,mean_CVB3_MSE,spec.CVB3{:});
%-----------------------------------------
title(['Number of Monte Carlo = ',num2str(setting.MonteCarlo)]);
xlabel('Radius');
ylabel('Mean squared error (MSE) of cluster means');

legend
movegui(h.MSE,'northwest')

%//////////////////////////////////////////////////////////////////////////

h.Purity = figure;
hold on;
%-----------------------------------------
plot(setting.Radius,mean_kmean_Purity,spec.kmean{:});
plot(setting.Radius,mean_emL_Purity,spec.emL{:});
plot(setting.Radius,mean_emMU_Purity,spec.emMU{:});
plot(setting.Radius,mean_VB_Purity,spec.VB{:});
plot(setting.Radius,mean_CVB1_Purity,spec.CVB1{:});
plot(setting.Radius,mean_CVB2_Purity,spec.CVB2{:});
plot(setting.Radius,mean_CVB3_Purity,spec.CVB3{:});
%-----------------------------------------
title(['Number of Monte Carlo = ',num2str(setting.MonteCarlo)]);
xlabel('Radius');
ylabel('Percentage of successful classifications (Purity) (%)');

legend
movegui(h.Purity,'north')

%//////////////////////////////////////////////////////////////////////////

h.ELBO = figure;
hold on;
%-----------------------------------------
plot(setting.Radius,mean_kmean_ELBO,spec.kmean{:});
plot(setting.Radius,mean_emL_ELBO,spec.emL{:});
plot(setting.Radius,mean_emMU_ELBO,spec.emMU{:});
plot(setting.Radius,mean_VB_ELBO,spec.VB{:});
plot(setting.Radius,mean_CVB1_ELBO,spec.CVB1{:});
plot(setting.Radius,mean_CVB2_ELBO,spec.CVB2{:});
plot(setting.Radius,mean_CVB3_ELBO,spec.CVB3{:});
%-----------------------------------------
title(['Number of Monte Carlo = ',num2str(setting.MonteCarlo)]);
xlabel('Radius');
ylabel('Evidence lower bound (ELBO)');

legend
movegui(h.ELBO,'northeast')

%//////////////////////////////////////////////////////////////////////////

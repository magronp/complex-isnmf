% Plot the influence of the phase parameters on the separation quality for
% complex ISNMF, complex NMF and consistent anisotropic Wiener filters

clc; clear all; close all;
test_or_dev = 'Dev';
set_settings_cisnmf;

%%% CISNMF
load(strcat(metrics_path,'learning_cisnmf.mat'));
SDR = squeeze(mean(score(:,:,:,end,1),3)) ;
SIR = squeeze(mean(score(:,:,:,end,2),3)) ;
SAR = squeeze(mean(score(:,:,:,end,3),3)) ;

Nt = length(Tau);
SDR(1,:) = repmat(SDR(1,1),[1 Nt]);
SIR(1,:) = repmat(SIR(1,1),[1 Nt]);
SAR(1,:) = repmat(SAR(1,1),[1 Nt]);

figure;
gcmap = colormap(autumn); gcmap = gcmap(end:-1:1,:); 
gcmap(:,2,:) = 1 - ((0:63)/64).^3;
colormap(gcmap);

subplot(1,3,1); imagesc(SDR); axis xy;
title('SDR (dB)','fontsize',14); set(gca,'xticklabel',Tau,'yticklabel',Kappa); xlabel('\tau','fontsize',16); ylabel('\kappa','fontsize',16);
subplot(1,3,2); imagesc(SIR); axis xy;
title('SIR (dB)','fontsize',14); set(gca,'xticklabel',Tau,'yticklabel',[]); xlabel('\tau','fontsize',16); 
subplot(1,3,3); imagesc(SAR); axis xy;
title('SAR (dB)','fontsize',14); set(gca,'xticklabel',Tau,'yticklabel',[]); xlabel('\tau','fontsize',16); 


%Without SAR
figure;
gcmap = colormap(autumn); gcmap = gcmap(end:-1:1,:); 
gcmap(:,2,:) = 1 - ((0:63)/64).^3;
colormap(gcmap);

subplot(1,2,1); imagesc(SDR); axis xy;
title('SDR (dB)','fontsize',14); set(gca,'xticklabel',Tau,'yticklabel',Kappa); xlabel('\tau','fontsize',16); ylabel('\kappa','fontsize',16);
subplot(1,2,2); imagesc(SIR); axis xy;
title('SIR (dB)','fontsize',14); set(gca,'xticklabel',Tau,'yticklabel',[]); xlabel('\tau','fontsize',16); 


%%% Wieners
load(strcat(metrics_path,'learning_wieners.mat'));
score_av = mean(score,4);

figure;
subplot(1,3,1); imagesc(score_av(:,:,1)); title('SDR (dB)','FontSize',16);
subplot(1,3,2); imagesc(score_av(:,:,2)); title('SIR (dB)','FontSize',16);
subplot(1,3,3); imagesc(score_av(:,:,3)); title('SAR (dB)','FontSize',16);


%%% CNMF
load(strcat(metrics_path,'learning_cnmf.mat'));
score_av = squeeze(mean(score,3));

figure;
subplot(1,3,1); plot(score_av(:,1)); title('SDR (dB)','FontSize',16); set(gca,'xticklabel',{'10^{-3}','10^{-1}','10^{1}','10^{3}'}); xlabel('\sigma_u');
subplot(1,3,2); plot(score_av(:,2)); title('SIR (dB)','FontSize',16); set(gca,'xticklabel',{'10^{-3}','10^{-1}','10^{1}','10^{3}'}); xlabel('\sigma_u');
subplot(1,3,3); plot(score_av(:,3)); title('SAR (dB)','FontSize',16); set(gca,'xticklabel',{'10^{-3}','10^{-1}','10^{1}','10^{3}'}); xlabel('\sigma_u');

clear all; clc; close all;

Fs = 44100; Nfft = 2048; Nw = Nfft; hop = Nfft/4;
dataset_path = 'datasets/DSD100/';
testordev = 'Dev';

% Get data
Nsongs = 50;
kappa = zeros(Nsongs,4);

% Estimate the VM concentration parameters
for it=1:Nsongs
    clc; fprintf('Song %d / %d \n',it,Nsongs);
    [~,~,Sm] = get_data_DSD(dataset_path,testordev,it,Fs,Nfft,Nw,hop);
    
    for j=1:4
        kaux = estim_kappa_vm(Sm(:,:,j),Nfft,hop);
        kappa(it,j) = kaux;
    end
    
% For Harmo and Percu
%      Saux1 = Sm(:,:,2);
%      Saux2 = sum(Sm(:,:,[1 3 4]),3);
%      kappa(it,1) = estim_kappa_vm(Saux1,Nfft,hop);
%      kappa(it,2) = estim_kappa_vm(Saux2,Nfft,hop);
    
end

% Plot
kappa(sum(isnan(kappa),2)==1,:) = [];
figure;
boxplot(kappa,'symbol',''); ylabel('\kappa','fontsize',16); set(gca,'xticklabel',{'bass';'drums';'other';'vocals'},'fontsize',16);


clear all; clc; close all;

Fs = 44100;
winleng = 2.^(9:13);
Nwinleg = length(winleng);

Nfft = 4096; Nw = Nfft; hop = Nfft/4;
dataset_path = 'datasets/DSD100/';
testordev = 'Dev';
Nsongs = 50;
kappa = zeros(4,Nwinleg);

for n=1:Nwinleg
    
    clc; fprintf('win length %d / %d \n',n,Nwinleg);
    Nw = winleng(n);
    Nfft = Nw;
    hop = Nfft/4;
    
    % Get data
    it=1;
    [~,~,Sm] = get_data_DSD(dataset_path,testordev,it,Fs,Nfft,Nw,hop);
    S = Sm;
    for it=1:Nsongs
        [~,~,Sm] = get_data_DSD(dataset_path,testordev,it,Fs,Nfft,Nw,hop);
        S = cat(2,S,Sm);
    end
    
    % Estimate kappa
    for j=1:4
        kaux = estim_kappa_vm(S(:,:,j),Nfft,hop);
        kappa(j,n) = kaux;
    end
    
end

save('vonmises/results/kappa_vs_nw_dsd.mat','kappa','winleng');

% Plot kappa vs. window length
figure;
for j=1:4
    subplot(1,4,j);
    plot(1:Nwinleg,kappa(j,:),'b-*');
    set(gca,'xticklabel',round(winleng/Fs*1000));
    xlabel('Window length (ms)','fontsize',16); ylabel('\kappa','fontsize',16);
end;

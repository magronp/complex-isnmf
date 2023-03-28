clear all; clc; close all;
%Fs = 44100; winleng = 2.^(9:14); Nwinleg = length(winleng);
Fs = 8000; winleng = 96*2.^(0:5); Nwinleg = length(winleng);


kappa_piano = zeros(1,Nwinleg);
kappa_guitar = zeros(1,Nwinleg);

for n=1:Nwinleg
    
    clc; fprintf('win length %d / %d \n',n,Nwinleg);
    
    Nw = winleng(n);
    Nfft = Nw;
    hop = Nfft/4;

    % Get piano data
    X = [];
    for it=1:30
        [~,Y] = get_data_MAPS_piece(Fs,Nfft,Nw,hop,it,[0 30]);
        X = [X Y];
    end

    % Estimate the VM concentration parameter
    kaux = estim_kappa_vm(X,Nfft,hop);
    kappa_piano(n) = kaux;
    
    % Get guitar data
    X = [];
    for it=1:6
        [~,Y] = get_data_guitar_piece(Fs,Nfft,Nw,hop,it);
        X = [X Y];
    end
    
    % Estimate the VM concentration parameter
    kaux = estim_kappa_vm(X,Nfft,hop);
    kappa_guitar(n) = kaux;
end

save('vonmises/results/kappa_vs_nw.mat','kappa_guitar','kappa_piano','winleng');

% Plot kappa vs. window length
figure;
subplot(1,2,1); plot(1:Nwinleg,kappa_piano,'b-*');
set(gca,'xticklabel',round(winleng/Fs*1000)); title('Piano');
xlabel('Window length (ms)','fontsize',16); ylabel('\kappa','fontsize',16);
subplot(1,2,2); plot(1:Nwinleg,kappa_guitar,'b-*');
set(gca,'xticklabel',round(winleng/Fs*1000)); title('Guitar');
xlabel('Window length (ms)','fontsize',16); ylabel('\kappa','fontsize',16);


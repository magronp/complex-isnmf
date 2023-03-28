clear all; clc; close all;
Fs = 8000; Nfft = 1000; Nw=Nfft; hop = Nfft/4;

% Gen data
[x,X,F,T,ts,freq] = get_data_MAPS_piece(Fs,Nfft,Nw,hop,1,[60 70]);
V = abs(X); phX = angle(X);

% Estim kappa and get peaks
[kappa,ph_centered,phaux,f_centr] = estim_kappa_vm(X,Nfft,hop);

Tagueule = cell(1,F); Kappa_freq = zeros(1,F);
for f=1:F
    phc = phaux(f,f_centr(f,:)==1);
    Tagueule{f} = phc;
    % ML estimation of kappa
    if ~isempty(phc)
        kappa = fzero(@(x) besseli(1,x)/besseli(0,x)-mean(cos(phc)),1);
    else
        kappa =nan;
    end
    Kappa_freq(f) = kappa;
end

Kappa_freq(isnan(Kappa_freq))=[];
boxplot(Kappa_freq);
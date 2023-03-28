clear all; clc; close all;
Fs = 8000; Nfft = 744; Nw=Nfft; hop = Nfft/4;
Nbins = 100;

% Gen data
%[sm,x,Sm,X,F,T,ts,freq,f1_Hz,f2_Hz] = gen_synthetic('no_ol',Fs,Nfft,Nw,hop);

[x,X,F,T,ts,freq] = get_data_MAPS_piece(Fs,Nfft,Nw,hop,1,[60 70]);
V = abs(X); phX = angle(X);

% Estim kappa and get peaks
[kappa,ph_centered,phaux,f_centr] = estim_kappa_vm(X,Nfft,hop);

% Spectrogram
h = figure;
imagesc(ts,freq,log10(V(:,1:end)+0.01)); axis xy; ylabel('Frequency (Hz)','fontsize',16); xlabel('Time (s)','fontsize',16);
h.Position = [980 680 387 285];

% Phasogram
phi = phX(f_centr==1);
h = figure; colormap(hsv);
imagesc(ts,freq,phX(:,1:end)); axis xy; ylabel('Frequency (Hz)','fontsize',16); xlabel('Time (s)','fontsize',16);
colorbar;
h.Position = [980 680 387 285];

% Phase histogram
h = figure;
histogram(phi,Nbins,'Normalization','pdf'); xlabel('\phi','fontsize',16); ylabel('Relative frequency','fontsize',16);
h.Position = [980 680 387 285];

% Centered phasogram
h = figure; colormap(hsv);
imagesc(ts,freq,phaux(:,1:end)); axis xy; ylabel('Frequency (Hz)','fontsize',16); xlabel('Time (s)','fontsize',16);
phiplot = linspace(-pi,pi,Nbins);
p_theo = exp(kappa * cos(phiplot))/(2*pi*besseli(0,kappa));
colorbar;
h.Position = [980 680 387 285];

% Centered phase histogram
h = figure;
histogram(ph_centered,Nbins,'Normalization','pdf'); hold on;
plot(phiplot,p_theo,'r'); xlabel('\psi','fontsize',16); ylabel('Relative frequency','fontsize',16);
h.Position = [980 680 387 285];

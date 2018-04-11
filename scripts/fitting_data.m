% This script performs the fitting test to find kappa such that the AG
% distribution best represent the data

clc; clear all; close all;
test_or_dev = 'Test';
set_settings_cisnmf;

%%% Load data
num_piece = 2;
[sm,x,Sm,X] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
[F,T,J] = size(Sm);

%%% Dictionnary learning on isolated sources
fprintf('Dico learning... \n');
Ktot = K*J;
W = cell(1,J);
W_matrix = zeros(F,Ktot);
for j=1:J
    waux = NMF(abs(Sm(:,:,j)).^2,rand(F,K),rand(K,T),iter_dico,0,0); 
    W{j}=waux;
    W_matrix(:,(j-1)*K+1:j*K) = waux;
end

%%% Initial ISNMF - H
fprintf('Initial ISNMF... \n');
Hinit = rand(Ktot,T);
[~,H_nmf_matrix] = NMF(abs(X).^2,W_matrix,Hinit,iter_init,0,0,1,ones(F,T),0);

% Get H onto cell form
H_nmf = cell(1,J);
for j=1:J
    H_nmf{j} = H_nmf_matrix((j-1)*K+1:j*K,:);
end

%%% Estimate normalized variables and get the histograms
Kappa = [0 0.01 0.05 0.1:0.2:1];
Nk = length(Kappa);
Nbins = 300;
histcis = zeros(Nk,Nbins); binscis = zeros(Nk,Nbins);

for kap=2:Nk
    kappa = Kappa(kap);

    % CISNMF
    clc; fprintf('Complex ISNMF... \n Kappa %d / %d \n',kap,Nk);
    [Xcisnmf,Wcis,Hcis,mucis] = complex_isnmf(X,W,H_nmf,angle(Sm),kappa,0,hop,iter_sep,zeros(F,T,J),0);

    % Variances
    vcis = zeros(F,T,J);
    for j=1:J
        vcis(:,:,j) = Wcis{j}*Hcis{j};
    end

    % AG moments
    lambda = besseli(1,kappa) / besseli(0,kappa) *sqrt(pi)/2;
    rho=besseli(2,kappa)./besseli(0,kappa) - lambda.^2;
    m = lambda * sqrt(vcis) .* exp(1i*mucis);
    gamma = (1-lambda^2)*vcis;
    c = rho*vcis .* exp(2i*mucis);
    m_X = sum(m,3);
    gamma_X = sum(gamma,3);      
    c_X = sum(c,3);
    detGamma = gamma_X.^2 - abs(c_X).^2;

    % normalized variables
    Xnormalized_cis = 2*(gamma_X .* abs(X-m_X).^2 - real(conj(c_X).*(X-m_X).^2) ) ./ (detGamma+eps);
    
    % histograms
    [auxh,auxb] = hist(Xnormalized_cis(:),Nbins);
    histcis(kap,:) =auxh; binscis(kap,:) =auxb;

end

 % ISNMF - Wiener
fprintf('ISNMF-Wiener... \n');
[~,H_isnmf] = NMF(abs(X).^2,W_matrix,Hinit,iter_init+iter_sep,0,0,1,ones(F,T),0);
vis = zeros(F,T,J);
for j=1:J
    vis(:,:,j) = W_matrix(:,(j-1)*K+1:j*K)*H_isnmf((j-1)*K+1:j*K,:);
end

Xnormis = 2*abs(X).^2 ./(sum(vis,3)+eps);
[auxh,auxb] = hist(Xnormis(:),Nbins);
histcis(1,:) =auxh; binscis(1,:) =auxb;

% Plot
xmax = 15; rnge = 0:0.01:xmax;
chi2density = chi2pdf(rnge,2) /0.5;

figure;
semilogy(binscis',(histcis./repmat(histcis(:,1),[1 Nbins]))'); V =axis; axis([0 xmax 10^-4 1]); hold on;
semilogy(rnge,chi2density,'k--');
hl=legend(string(Kappa)); set(hl,'fontsize',14);
xlabel('y','fontsize',16); ylabel('p(y)','fontsize',16);
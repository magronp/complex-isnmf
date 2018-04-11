% Influence of kappa and delta (anisotropy and consistency weight)
% for the consistent anisotropic Wiener filter

clc; clear all; close all;
test_or_dev = 'Dev';
set_settings_cisnmf;

% Learning parameters
Delta = [0 10.^(-3:3)]; Nd = length(Delta);
Kappa = [0 10.^(-2:1:2)]; Nk = length(Kappa);

score = zeros(Nk,Nd,3,Nsongs);

for ind=1:Nsongs
    
    clc; fprintf('Data %d / %d \n',ind,Nsongs);
    % Load data
    num_piece = datavec(ind);
    [sm,x,Sm,X] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);

    %%% Variances estimates - semi-informed setting
    
    % dictionaries from the isolated power spectra
    W = zeros(F,Ktot);
    for j=1:J
        waux = NMF(abs(Sm(:,:,j)).^2,rand(F,K),rand(K,T),iter_dico,0,0);   % ISNMF in the semi-informed setting
        W(:,(j-1)*K+1:j*K)=waux;
    end

    % separation: NMF on the mixture
    [~,H] = NMF(abs(X).^2,W,rand(Ktot,T),iter_init+iter_sep,0,0,1,ones(F,T),0);
    variances = zeros(F,T,J);
    for j=1:J
        variances(:,:,j) = W(:,(j-1)*K+1:j*K)*H((j-1)*K+1:j*K,:);
    end
    
    magnitude = sqrt(variances);
    Sm_approx = magnitude .* exp(1i*angle(X));
    
    % Consistent Anisotropic Wiener Filtering
    for kap=1:Nk
        for d=1:Nd
            clc; fprintf('Data %d / %d \n Kappa %d / %d \n Delta %d / %d \n',ind,Nsongs,kap,Nk,d,Nd);
            delta = Delta(d); kappa = Kappa(kap);
            
            Xaw = anisotropic_wiener(X,Sm_approx,kappa*ones(F,T,J),hop);
            Xcaw = caw(Xaw,variances,kappa,delta,Nw,hop,wtype);

            se = real(iSTFT(Xcaw,Nfft,hop,Nw,wtype));
            [sd,si,sa] = GetSDR(se,sm);
            score(kap,d,:,ind) = [mean(sd) mean(si) mean(sa)];        
        end
    end
end

% Record scores
save(strcat(metrics_path,'learning_wieners.mat'),'score','Delta','Kappa');

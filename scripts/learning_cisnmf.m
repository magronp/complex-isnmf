% Influence of kappa and tau (anisotropy and phase constraint parameter)
% for complex ISNMF

clc; clear all; close all;
test_or_dev = 'Dev';
set_settings_cisnmf;

% Learning parameters
Kappa = [0 0.1 0.5 1]; Nkappa = length(Kappa);
Tau = [0 0.1 0.5 1 5]; Ntau = length(Tau);

score = zeros(Nkappa,Ntau,Nsongs,iter_sep+1,3);

for ind=1:Nsongs

    % Load data
    num_piece = datavec(ind);
    [sm,x,Sm,X] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);

    % Dictionnary learning on isolated sources
    clc; fprintf('Data %d / %d \n Dico learning \n',ind,Nsongs);
    W = cell(1,J); Weuc = cell(1,J);
    W_matrix = zeros(F,Ktot);
    for j=1:J
        waux = NMF(abs(Sm(:,:,j)).^2,rand(F,K),rand(K,T),iter_dico,0,0); 
        W{j}=waux; W_matrix(:,(j-1)*K+1:j*K) = waux;
    end

    % Initial random activation matrices and ISNMF on the mixture  
    [~,Hnmf_matrix] = NMF(abs(X).^2,W_matrix,rand(Ktot,T),iter_init,0,0,1,ones(F,T),0);
    
    % Store the H matrices into cells for CISNMF
    Hnmf = cell(1,J);
    for j=1:J
        Hnmf{j} = Hnmf_matrix((j-1)*K+1:j*K,:);
    end
    
    % Initial phases and normalized frequencies (and H)
    muini = repmat(angle(X),[1 1 J]);
    nu = zeros(F,T,J);
    for j=1:J
        nu(:,:,j) = get_frequencies_qifft(W_matrix(:,(j-1)*K+1:j*K)*Hnmf_matrix((j-1)*K+1:j*K,:))/Nfft;
    end

    % CISNMF with kappa = 0
    clc; fprintf('Data %d / %d \n Complex ISNMF \n kappa= tau =0 \n',ind,Nsongs);
    [~,~,~,~,~,bss,se] = complex_isnmf(X,W,Hnmf,muini,0,0,hop,iter_sep,nu,0,1,sm,Nw,wtype);
    score(1,1,ind,:,:) = bss;

    % CISNMF - Influence of kappa and tau   
    for pk=2:Nkappa
        for pt=1:Ntau
            clc; fprintf('Data %d / %d \n Complex ISNMF \n kappa %d / %d \n tau %d / %d \n',ind,Nsongs,pk,Nkappa,pt,Ntau);
            kappa = Kappa(pk);
            tau = Tau(pt);
            
            [~,~,~,~,~,bss,se] = complex_isnmf(X,W,Hnmf,muini,kappa,tau,hop,iter_sep,nu,0,1,sm,Nw,wtype);
            score(pk,pt,ind,:,:) = bss;
        end
    end

end

save(strcat(metrics_path,'learning_cisnmf.mat'),'score','Kappa','Tau');
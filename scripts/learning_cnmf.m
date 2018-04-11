% Influence of sigma_u (phase constraint parameter) for complex NMF with
% the Euclidean distance

clc; clear all; close all;
test_or_dev = 'Dev';
set_settings_cisnmf;

% Learning parameters
Sigma_u = [0 10.^(-3:3)]; Nsigu = length(Sigma_u);
score = zeros(Nsigu,3,Nsongs);


for ind=1:Nsongs

    % Load data
    num_piece = datavec(ind);
    [sm,x,Sm,X] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);

    % Dictionnary learning on isolated sources
    clc; fprintf('Data %d / %d \n Dico learning \n',ind,Nsongs);
    W = cell(1,J);
    W_matrix = zeros(F,Ktot);
    for j=1:J
        waux = NMF(abs(Sm(:,:,j)),rand(F,K),rand(K,T),iter_dico,2,0);
        W{j}=waux; W_matrix(:,(j-1)*K+1:j*K) = waux;
    end

    
    % Initial random activation matrices and Euclidean NMF on the mixture  
    [~,Hnmf_matrix] = NMF(abs(X),W_matrix,rand(Ktot,T),iter_init,2,0,1,ones(F,T),0);
    
    % Store the H matrices into cells for complex NMF
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
    
    % Complex NMF
    for sig=1:Nsigu
        clc; fprintf('Data %d / %d \n Complex NMF \n sigma_u %d / %d \n',ind,Nsongs,sig,Nsigu);
        sigma_u = Sigma_u(sig);
        
        Xe=cnmf_wisa_ph(X,W,Hnmf,iter_sep,sigma_u,0,nu,hop);
        se = real(iSTFT(Xe,Nfft,hop,Nw,wtype));
        [sd,si,sa] = GetSDR(se,sm);
        score(sig,:,ind) = [mean(sd) mean(si) mean(sa)];
    end

end

% Record scores
save(strcat(metrics_path,'learning_cnmf.mat'),'score','Sigma_u');

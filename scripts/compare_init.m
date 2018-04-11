clc; clear all; close all;
test_or_dev = 'Dev';
set_settings_cisnmf;

score_i = zeros(iter_sep+1,3,Nsongs);
score_r = zeros(iter_init+iter_sep+1,3,Nsongs);

for ind=1:Nsongs

    clc; fprintf('Song %d / %d... \n',ind,Nsongs);
    
    % Source generation
    num_piece = datavec(ind);
    [sm,x,Sm,X] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);

    % Dictionnary learning on isolated sources
    fprintf('Dico learning... \n');
    Ktot = K*J;
    W = cell(1,J);
    W_matrix = zeros(F,Ktot);
    for j=1:J
        waux = NMF(abs(Sm(:,:,j)).^2,rand(F,K),rand(K,T),iter_dico,0,0); 
        W{j} = waux;
        W_matrix(:,(j-1)*K+1:j*K) = waux;
    end

    % ISNMF initialization on the mix
    fprintf('ISNMF init... \n');
    Hrand_matrix = rand(Ktot,T);
    [~,Hnmf_matrix] = NMF(abs(X).^2,W_matrix,Hrand_matrix,iter_init,0,0,1,ones(F,T),0);

    % Store the activations matrix into cells
    Hnmf = cell(1,J); Hrand = cell(1,J);
    for j=1:J
        Hnmf{j} = Hnmf_matrix((j-1)*K+1:j*K,:);
        Hrand{j} = Hrand_matrix((j-1)*K+1:j*K,:);
    end   
    
    % Initial phases and normalized frequencies
    muini = repmat(angle(X),[1 1 J]);
    nu = zeros(F,T,J);
    for j=1:J
        nu(:,:,j) = get_frequencies_qifft(W_matrix(:,(j-1)*K+1:j*K)*Hnmf_matrix((j-1)*K+1:j*K,:))/Nfft;
    end

    % Complex ISNMF - random initialization
    fprintf('Complex ISNMF 1 /2 ... \n');
    [~,~,~,~,~,aux] = complex_isnmf(X,W,Hnmf,muini,0.5,0.5,hop,iter_sep,nu,0,1,sm,Nw,wtype);
    score_i(:,:,ind) = aux;
    
    % Complex ISNMF - initialization with ISNMF
    fprintf('Complex ISNMF 2 /2 ... \n');
    [~,~,~,~,~,aux] = complex_isnmf(X,W,Hrand,muini,0.5,0.5,hop,iter_init+iter_sep,nu,0,1,sm,Nw,wtype);
    score_r(:,:,ind) = aux;

end

% Record scores
save(strcat(metrics_path,'bss_init.mat'),'score_i','score_r');


% Plot SDR over iterations
score_i_av = mean(score_i,3);
score_r_av = mean(score_r,3);

figure;
subplot(1,2,1); plot(0:iter_sep,score_i_av(:,1));
title('ISNMF initialization'); ylabel('SDR (dB)','fontsize',16); xlabel('Iterations','fontsize',16);
subplot(1,2,2); plot(0:iter_sep+iter_init,score_r_av(:,1));
title('Random initialization'); xlabel('Iterations','fontsize',16);

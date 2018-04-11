% Perform the separation on the test set for various algorithms

clc; clear all; close all;
test_or_dev = 'Test';
set_settings_cisnmf;

for ind=1:Nsongs

    clc; fprintf('Data %d / %d \n',ind,Nsongs);
    
    % Get data
    num_piece = datavec(ind);
    [sm,x,Sm,X] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);
    
    % Dictionnary learning on isolated sources
    fprintf('Dico learning... \n');
    Ktot = K*J;
    Wis = cell(1,J); Weu = cell(1,J);
    Wis_matrix = zeros(F,Ktot); Weu_matrix = zeros(F,Ktot);
    for j=1:J
        Wini = rand(F,K); Hini=rand(K,T);
        waux = NMF(abs(Sm(:,:,j)).^2,Wini,Hini,iter_dico,0,0);  Wis{j}=waux; Wis_matrix(:,(j-1)*K+1:j*K) = waux;
        waux = NMF(abs(Sm(:,:,j)),Wini,Hini,iter_dico,2,0);  Weu{j}=waux; Weu_matrix(:,(j-1)*K+1:j*K) = waux;
    end


    % Initial ISNMF and EuNMF on the mixture
    fprintf('Initial NMFs... \n');
    Hini_matrix = rand(Ktot,T);
    [~,His_nmf_matrix] = NMF(abs(X).^2,Wis_matrix,Hini_matrix,iter_init,0,0,1,ones(F,T),0);
    [~,Heu_nmf_matrix] = NMF(abs(X),Weu_matrix,Hini_matrix,iter_init,2,0,1,ones(F,T),0);

    
    % Store activation matrices in cells for complex NMFs
    His_nmf = cell(1,J); Heu_nmf = cell(1,J);
    for j=1:J
        His_nmf{j} = His_nmf_matrix((j-1)*K+1:j*K,:);
        Heu_nmf{j} = Heu_nmf_matrix((j-1)*K+1:j*K,:);
    end
    
    %%% Initial phases and normalized frequencies
    muini = repmat(angle(X),[1 1 J]);
    nu_is = zeros(F,T,J); nu_eu = zeros(F,T,J);
    for j=1:J
        nu_is(:,:,j) = get_frequencies_qifft(Wis_matrix(:,(j-1)*K+1:j*K)*His_nmf_matrix((j-1)*K+1:j*K,:))/Nfft;
        nu_eu(:,:,j) = get_frequencies_qifft(Weu_matrix(:,(j-1)*K+1:j*K)*Heu_nmf_matrix((j-1)*K+1:j*K,:))/Nfft;
    end

    %%% Separation

    Se = zeros(F,T,J,Nalgo);
    
    fprintf('ISNMF-Wiener... \n');
    [~,H_isnmf] = NMF(abs(X).^2,Wis_matrix,Hini_matrix,iter_init+iter_sep,0,0,1,ones(F,T),0);
    variances = zeros(F,T,J);
    for j=1:J
        variances(:,:,j) = Wis_matrix(:,(j-1)*K+1:j*K)*H_isnmf((j-1)*K+1:j*K,:);
    end
    Se(:,:,:,1) = variances ./ (sum(variances,3)+eps) .* X;
    Sm_approx = sqrt(variances) .* exp(1i*repmat(angle(X),[1 1 J]));
    
    fprintf('ISNMF-AW... \n');
    Se(:,:,:,2) = anisotropic_wiener(X,Sm_approx,kappa_aw*ones(F,T,J),hop);
    
    fprintf('ISNMF-CAW... \n');
    aux = anisotropic_wiener(X,Sm_approx,kappa_caw*ones(F,T,J),hop);
    Se(:,:,:,3) = caw(aux,variances,kappa_caw,delta_caw,Nw,hop,wtype);
    
    fprintf('CNMF... \n');
    Se(:,:,:,4)=cnmf_wisa_ph(X,Weu,Heu_nmf,iter_sep,sigma_u,0,nu_eu,hop);
        
     % CISNMF
    fprintf('Complex ISNMF... \n');
    Se(:,:,:,5) = complex_isnmf(X,Wis,His_nmf,muini,kappa_cisnmf,tau_cisnmf,hop,iter_sep,nu_is,0);
    
    %%% Synthesis and record
    
    % Time-domain synthesis
    s_estim = zeros(J,length(sm),Nalgo);
    for al=1:Nalgo
        s_estim(:,:,al) = real(iSTFT(Se(:,:,:,al),Nfft,hop,Nw));
    end
    
    % Record
    for j=1:J
        audiowrite(strcat(audio_path,'song',int2str(ind),'_source',int2str(j),'_orig.wav'),sm(j,:),Fs)
        for al = 1:Nalgo
            audiowrite(strcat(audio_path,'song',int2str(ind),'_source',int2str(j),'_',algos{al},'.wav'),s_estim(j,:,al),Fs);
        end
    end

end
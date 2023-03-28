clc; clearvars; close all;

% Parameters
Fs = 44100;
kappa_aw=1.6;

Nw = 4096; kappa_aw_var=[2.28 1.26 1.51 1.30];    % with Nw = 4096
%Nw = 2048; kappa_aw_var=[2.16 1.33 1.7 1.39];    % with Nw = 2048

Nfft = Nw; hop = Nw/4;
direc = 'vonmises/audio_files/';

% Initialize data vec (remove useless mixtures)
Nsongs = 50;
datavec = 1:Nsongs;
dataNaN=[6 34 36 40]; datavec(dataNaN)=[];
datavec = datavec(1:min(Nsongs,length(datavec))); Nsongs = length(datavec);
score = zeros(4,3,3,Nsongs); SIR = zeros(4,3,3,Nsongs); SAR = zeros(4,3,3,Nsongs);


for it =1:Nsongs
   
    % Source generation
    clc; fprintf('Data %d / %d \n',it,Nsongs);
    num_piece = datavec(it);
    [sm,x,Sm,X,ts,freq] = get_data_DSD('datasets/DSD100/','Test',num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);
    v = abs(Sm).^2;
    
    % Wiener filtering
    Xw = v .* X./(sum(v,3)+eps);
    
    %  Anisotropic Wiener
    Xaw = anisotropic_wiener(X,v,kappa_aw*ones(F,T,J),hop);
    
    %  Anisotropic Wiener - variable kappa
    kap = zeros(F,T,J);
    for j=1:J
        kap(:,:,j) = kappa_aw_var(j);
    end
    Xaw_var = anisotropic_wiener(X,v,kap,hop);
    
    % Synthesis
    sw = real(iSTFT(Xw,Nfft,hop,Nw));
    saw = real(iSTFT(Xaw,Nfft,hop,Nw));
    saw_var = real(iSTFT(Xaw_var,Nfft,hop,Nw));
    
    % Score
    [sdr,sir,sar] = GetSDR(sw,sm); score(:,1,:,it) = [(sdr) (sir) (sar)];
    [sdr,sir,sar] = GetSDR(saw,sm); score(:,2,:,it) = [(sdr) (sir) (sar)];
    [sdr,sir,sar] = GetSDR(saw_var,sm); score(:,3,:,it) = [(sdr) (sir) (sar)];

end

save('vonmises/results/ssep_bss.mat','score');
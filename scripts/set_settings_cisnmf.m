% Set the settings used in the experiments for complex ISNMF

Fs = 44100;

% Data
Nsongs = 50;
datavec = 1:Nsongs;
switch test_or_dev
    case 'Dev'
        dataNaN=[1 18 29 38 49];
    case 'Test'
        dataNaN=[6 34 36 40];
end
datavec(dataNaN)=[];
datavec = datavec(1:min(Nsongs,length(datavec)));
Nsongs = length(datavec);

J = 4;

% STFT parameters
Nfft = 4096;
Nw = 4096;
hop = Nw/4;
wtype = 'hann';

% Paths
dataset_path = 'datasets/DSD100/';
audio_path = 'complex-isnmf/audio_files/';
metrics_path = 'complex-isnmf/metrics/';

% Algorithms
algos = {'isnmf-W','isnmf-AW','isnmf-CAW','cnmf','cisnmf'};
Nalgo = length(algos);

% NMF
K = 50; Ktot = K*J;

% Phase parameters
kappa_caw = 0.1; delta_caw = 10^-3; 
kappa_aw = 1;
sigma_u = 10^-2;
kappa_cisnmf = 0.5; tau_cisnmf=1;

% Iterations
iter_dico = 200;
iter_init = 50;
iter_sep = 100;
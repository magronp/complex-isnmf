% Compute the score for the benchmarked algorithms

clc; clear all; close all;
test_or_dev = 'Test';
set_settings_cisnmf;

% Init signals and scores
L = 440320;
sm = zeros(J,L);
se = zeros(J,L,Nalgo);
score = zeros(Nalgo,Nsongs,J,3);

for ind=1:Nsongs
    
    clc; fprintf('Data %d / %d \n',ind,Nsongs);
     
    % Original Files path
    for j=1:J
       sm(j,:) = audioread(strcat(audio_path,'song',int2str(ind),'_source',int2str(j),'_orig.wav'));
       for al=1:Nalgo
           se(j,:,al) = audioread(strcat(audio_path,'song',int2str(ind),'_source',int2str(j),'_',algos{al},'.wav'));
       end
    end
    
    % BSS Eval
    for al=1:Nalgo
       [auxd,auxi,auxa] = GetSDR(squeeze(se(:,:,al)),sm);
       score(al,ind,:,:) = [auxd auxi auxa];
    end   
    
end

% Record scores
save(strcat(metrics_path,'separation_benchmark.mat'),'score','algos');

% Consistent anisotropic Wiener filtering.
%
% Ref:
% "Consistent anisotropic Wiener filtering for audio source separation",
% Paul Magron, Jonathan Le Roux and Tuomas Virtanen
% Proc. of the IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA)
% October 2017
% 
% Inputs:
%     Xaw : F*T*J initial STFTs obtained e.g. with anisotropic Wiener filtering
%     variances : F*T*J estimates of the sources' variances
%     kappa >= 0 : anisotropy parameter
%     delta >= 0 : consistency parameter
%     Nw : STFT window length
%     hop : hop size of the STFT
%     wtype : analysis window type
%     tol_pcg : stopping criterion for the pcg algorithm
%     max_iter : maximum number of iterations for pcg
%     sm : ground truth sources for computing the score over iterations
%
% Outputs:
%     Xe : estimated STFT sources
%     iter : number of iterations of PCG for each source
%     score : BSS eval over iterations for each source

function [Xe,iter,score] = caw(Xaw,variances,kappa,delta,Nw,hop,wtype,tol_pcg,max_iter,sm)

[F,T,J] = size(variances);

bss = 1;
if nargin<10
    sm=zeros(J,1); bss=0;
end

if nargin<9
    max_iter = 100;
end

if nargin<8
    tol_pcg = 1e-6;
end

% Prepare outputs
Xe = zeros(F,T,J);
iter = cell(1,J-1);
score = cell(1,J-1);

var_single = zeros(F,T,2);
mu_single = zeros(F,T,2);
sm_single = zeros(2,length(sm));

% Loop over the sources (one vs. the other in turns)
for j=1:J-1

    % Get variance, phase and time domain sources for the current 2 spirces
    var_single(:,:,1) = variances(:,:,j)+eps; 
    var_single(:,:,2) = sum(variances(:,:,j+1:end),3)+eps; 
    mu_single(:,:,1) = angle(Xaw(:,:,j));
    mu_single(:,:,2) = angle(sum(Xaw(:,:,j+1:end),3));
    
    sm_single(1,:) = sm(j,:); 
    sm_single(2,:) = sum(sm(j+1:end,:),1); 
    
    % CAW for 2 sources
    [SE,iter_single,score_single] = caw2sources(Xaw(:,:,j),kappa,mu_single,var_single,delta,Nw,hop,wtype,tol_pcg,max_iter,bss,sm_single);
    
    % Collect outputs
    iter{j} = iter_single;
    score{j} = score_single;
    Xe(:,:,j) = SE;
end

% Last source
Xe(:,:,J) = sum(Xaw,3)-sum(Xe(:,:,1:J-1),3);

end

% CAW for 2 sources only
function [SE,iter,score] = caw2sources(Sini,kappa,mu,variances,delta,Nw,hop,wtype,tol_pcg,max_iter,bss,sm)

[F,T]=size(Sini);
Nfft=2*(F-1);

% Anisotropy weights 
lambda = besseli(1,kappa) ./ besseli(0,kappa);
rho = (besseli(2,kappa).*besseli(0,kappa) - besseli(1,kappa).^2 )./ besseli(0,kappa).^2;

% Covariance in the AG mdoel
gamma = (1-lambda.^2).* variances;
c = rho.*variances .*exp(2*1i*mu) ;

% Mixture and posterior covariance
gamma_X = sum(gamma,3);      
c_X = sum(c,3);
detGX = gamma_X.^2 - abs(c_X).^2+eps;
        
gamma_post = abs(gamma(:,:,1) - ( gamma_X .* (gamma(:,:,1).^2+abs(c(:,:,1)).^2) - 2 * gamma(:,:,1) .* real(c(:,:,1).*conj(c_X))  )  ./ detGX);
c_post = c(:,:,1) - (2*gamma(:,:,1).*gamma_X.*c(:,:,1) - gamma(:,:,1).^2 .* c_X - c(:,:,1).^2 .* conj(c_X)  )  ./ detGX;
detG = gamma_post.^2-abs(c_post).^2;


%%% Wiener filter initialization %%%
SE=Sini;

%%% Preconditioned conjugate gradient %%%
wei=repmat([1; 2*ones(F-2,1); 1],[1 T]);
se=real(iSTFT(SE,Nfft,hop,Nw,wtype));
FSE=SE-STFT(se,Nfft,hop,Nw,wtype);
r=-delta*FSE;

z =  detG ./ (1+2*delta*(Nw-hop)/Nw*gamma_post + (delta*(Nw-hop)/Nw)^2 * detG) .* ( (gamma_post./(detG+eps) + delta*(Nw-hop)/Nw ).* r + c_post./(detG+eps) .* conj(r) );
P=z;
rsold=real(sum(sum(wei.*conj(r).*z)));

converged=false;

score = zeros(2,3,max_iter+1);
if bss
    se2 = sum(sm,1)-se;
    [sd,si,sa] = GetSDR([se;se2],sm);
    score(:,:,1) = [sd si sa];
end

% Loop
iter=0;
while and(~converged,iter<max_iter)
    iter=iter+1;
    p=real(iSTFT(P,Nfft,hop,Nw,wtype));
    
        
    FP=P-STFT(p,Nfft,hop,Nw,wtype);
    AP=( gamma_post .* P - c_post .* conj(P)  )./(detG+eps) +delta*FP;
    alpha=rsold/real(sum(sum(wei.*conj(P).*AP))+realmin);
    SE=SE+alpha*P;
    
    if bss
        se1 = real(iSTFT(SE,Nfft,hop,Nw,wtype)); se2 = sum(sm,1)-se1;
        [sd,si,sa] = GetSDR([se1;se2],sm);
        score(:,:,iter+1) = [sd si sa];
    end
        
    
    converged=(sum(sum(alpha^2*real(P.*conj(P)))) < tol_pcg*sum(sum(real(SE.*conj(SE)))));
    r=r-alpha*AP;
    z =  detG ./ (1+2*delta*(Nw-hop)/Nw*gamma_post + (delta*(Nw-hop)/Nw)^2 * detG) .* ( (gamma_post./(detG+eps) + delta*(Nw-hop)/Nw ).* r + c_post./(detG+eps) .* conj(r) );
    
    rsnew=real(sum(sum(wei.*conj(r).*z)));
    beta=rsnew/(rsold+realmin);
    P=z+beta*P;
    rsold=rsnew;
    
end

% Remove non-relevant scores
score = score(:,:,1:max_iter+1);

end
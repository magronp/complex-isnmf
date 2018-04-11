% Complex ISNMF - EM algorithm based on the anisotropic Gaussian model
%
% Ref:
% "Complex ISNMF: a phase-aware model for monaural audio source separation",
% Paul Magron and Tuomas Virtanen,
% IEEE Transactions on Audio, Speech, and Language Processing, 2018
%
% Inputs:
%     X : F*T mixture
%     Wini : a cell of J F*K initial dictionaries
%     Hini : a cell of J K*T initial activations matrices
%     muini : F*T*J initial phase fields
%     kappa : anisotropy parameter
%     tau : prior phase weight
%     hop : hop size (samples) of the STFT
%     Niter : number of iterations
%     nu : F*T*J normalized frequencies
%     compute_cost : if =1, the cost and BSSeval scores are computed
%     sm : J*L original time-domain sources
%     Nw : analysis window length
%     upW : W is updated only if =1
%
% Outputs:
%     m_post : estimated components (posterior means)
%     W : dictionaries
%     H : activation matrices
%     mu : estimated phases
%     C_MAP : MAP functional values over iterations
%     scorebss : SDR, SIR and SAR over iterations
%     se : estimated time-domain sources

function [m_post,W,H,mu,C_MAP,scorebss,se] = complex_isnmf(X,Wini,Hini,muini,kappa,tau,hop,Niter,nu,upW,compute_cost,sm,Nw,wtype)

[F,T,J] = size(muini);
Nfft = 2*(F-1);

if nargin<11
    compute_cost = 0;
end

if nargin<10
    upW=1;
end

% Anisotropy parameters
lambda = besseli(1,kappa) / besseli(0,kappa) *sqrt(pi)/2;
rho=besseli(2,kappa)./besseli(0,kappa) - lambda.^2;

% Initial parameters
W = Wini; H = Hini;
v = zeros(size(muini));
for j=1:J
   v(:,:,j) = W{j}*H{j}; 
end
mu = muini;

% Initial moments
m = lambda * sqrt(v) .* exp(1i*mu);
gamma = (1-lambda^2)*v;
c = rho*v .* exp(2i*mu);

% Cost
if compute_cost
    C_MAP = zeros(1,Niter+1);
    scorebss = zeros(Niter+1,3);
end


% EM algo
for iter=1:Niter
    
    % --------------------------------- E step

    % ----- update estimated sum moments
    m_X = sum(m,3);
    gamma_X = sum(gamma,3);      
    c_X = sum(c,3);
    detGamma = gamma_X.^2 - abs(c_X).^2;

    % ----- posterior moments
    m_post = m + ( (gamma.*gamma_X - c.*conj(c_X)).*(X-m_X) + (c.*gamma_X - gamma.*c_X).* conj(X-m_X)  ) ./ (detGamma+eps) ;
    gamma_post = gamma - ( gamma_X .* (gamma.^2+abs(c).^2) - 2 * gamma .* real(c.*conj(c_X))  )  ./ (detGamma+eps);
    c_post = c - (2*gamma.*gamma_X.*c - gamma.^2 .* c_X - c.^2 .* conj(c_X)  )  ./ (detGamma+eps);


    % Compute the cost
    if compute_cost
        se = real(iSTFT(m_post,Nfft,hop,Nw,wtype));
        C_MAP(iter) = compute_cmap(X,m,gamma,c,mu,tau,hop,nu) ;
        [sd,si,sa] = GetSDR(se,sm);
        scorebss(iter,:) = [mean(sd) mean(si) mean(sa) ];
    end
    
    
    % --------------------------------- M step

    % ----- auxiliary parameters
    p = ((1-lambda^2)*(gamma_post+abs(m_post).^2)-rho*real(exp(-2*1i*mu).*(c_post+m_post.^2))  )/((1-lambda^2)^2-rho^2)  ;
    q = 2*lambda*(rho-(1-lambda^2))* real( m_post .* exp(-1i*mu) ) /((1-lambda^2)^2-rho^2);
    q(q>0)=0;
    
    % NMF updates
    for j=1:J
        
        % Get current source parameters
        vj = v(:,:,j)+eps;
        pj = p(:,:,j);
        qj = q(:,:,j);

        % Update H
        H{j} =  H{j} .* sqrt( ( W{j}' *(pj.*vj.^(-2)))  ./ (  W{j}' * ( vj.^(-1) - qj.*vj.^(-1.5) ) ) );
        vj= W{j}*H{j}+eps;
        
        % Update W
        if upW
            W{j} =  W{j} .* sqrt( ( (pj.*vj.^(-2))* H{j}' )  ./ (  (vj.^(-1) - qj.*vj.^(-1.5) )* H{j}'  )  );
            sumWj = sqrt(sum(W{j}).^2);
            W{j} = W{j} ./repmat(sumWj,[F 1]) ;
            H{j} = repmat(sumWj',[1 T]) .* H{j};           
        end
        
        v(:,:,j)= W{j}*H{j};
    end
    
        
    % ----- phase param
    beta_part = 2*lambda*(1-lambda^2-rho) ./  (  ( (1-lambda^2)^2-rho^2 ) * sqrt(v) +eps ) .* m_post   ;
    for t=2:T-1
       mu(:,t,:)=angle( beta_part(:,t,:) + tau*(  exp(1i* (mu(:,t-1,:)+2*pi*hop*nu(:,t,:) ) ) +  exp(1i* (mu(:,t+1,:)-2*pi*hop*nu(:,t+1,:)))     ) );
    end

    % ----- update m, gamma and c
    m = lambda * sqrt(v) .* exp(1i*mu);
    gamma = (1-lambda^2)*v;
    c = rho*v .* exp(2i*mu);

end

% Finally, one E step for having the most up to date posterior mean
m_X = sum(m,3);
gamma_X = sum(gamma,3);      
c_X = sum(c,3);
detGamma = gamma_X.^2 - abs(c_X).^2;
m_post = m + ( (gamma.*gamma_X - c.*conj(c_X)).*(X-m_X) + (c.*gamma_X - gamma.*c_X).* conj(X-m_X)  ) ./ (detGamma+eps);
    
% Last cost
if compute_cost
    se = real(iSTFT(m_post,Nfft,hop,Nw,wtype));
    C_MAP(iter+1) = compute_cmap(X,m,gamma,c,mu,tau,hop,nu) ;
    [sd,si,sa] = GetSDR(se,sm);
    scorebss(iter+1,:) = [mean(sd) mean(si) mean(sa) ];
end
    
end



function [cmap] = compute_cmap(X,m,gamma,c,mu,tau,hop,nu) 

[F,~,J] =size(m);

m_X = sum(m,3);
gamma_X = sum(gamma,3);      
c_X = sum(c,3);
detGamma = gamma_X.^2 - abs(c_X).^2;

logc = log(detGamma+eps)/2 + ( gamma_X .* abs(X-m_X).^2  + real( (X-m_X).^2 .* conj(c_X) ) ) ./ (detGamma+eps);

logphi = 0;
for j=1:J
    mu_m = [zeros(F,1) mu(:,1:end-1,j)  ];
    mu_p = [mu(:,2:end,j) zeros(F,1)]; 
    logphi =logphi +  tau* real( exp(-1i*mu(:,:,j)) .* ( exp( 1i* (mu_m+2*pi*hop*nu(:,:,j) ) ) +exp(1i* (mu_p-2*pi*hop*nu(:,:,j))) ) );
end

cmap = sum(logphi(:)-logc(:)) ;
   
end
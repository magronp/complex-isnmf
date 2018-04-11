% Source separation using anisotropic gaussian model with phase unwrapping
% prior
% Ref: "Phase-dependent anisotropic Gaussian model for audio source
% separation", P. Magron, R. Badeau, B. David, in Proc. ICASSP 2017
%
% Inputs:
%     X : F*T mixture
%     v : F*T*J variances
%     kappa : F*T*J concentration parameters
%     hop : STFT overlap (in samples)
%     UN : K*T onset indicatrix function
%     muini : F*T*J initial phases (if some are "known")
%
% Outputs:
%     Shat : estimated sources
%     mu : unwrapped phase (sinusoidal model)


function [Se,mu] = anisotropic_wiener(X,v,kappa,hop,UN,muini)

[F,T,J] = size(v);

if nargin<6
    muini = repmat(angle(X),[1 1 J]);
end

if nargin<5
    UN = zeros(J,T);
end


% Unwrapping Parameters
Nfft = (F-1)*2;
Ct = 2*pi*hop/Nfft;

% Weight parameters
lambda = besseli(1,kappa) ./ besseli(0,kappa);
rho = (besseli(2,kappa).*besseli(0,kappa) - besseli(1,kappa).^2 )./ besseli(0,kappa).^2;

% Initialization
mu = muini;
Se=sqrt(v) .* exp(1i * mu);

% Loop over time frames
for t=2:T
    if (sum(UN(:,t))<J) % if the current frame is onset for all sources, do nothing
        
        % Initialisation : if onset frame, do nothing, if not, unwrapping
        for j=1:J
             if (UN(j,t)==0) % non-onset frame for source k
                 f_inf = get_frequencies_qifft_frame(abs(Se(:,t,j))+eps);
                 mu(:,t,j) = angle(Se(:,t-1,j))+Ct*f_inf;
             end
            Se(:,t,j) = abs(Se(:,t,j)) .* exp(1i * mu(:,t,j));
        end
        
        % Compute prior, means, variances and covariances
        Xtilde = squeeze(Se(:,t,:));
        
        lambda_aux = squeeze(lambda(:,t,:));
        rho_aux = squeeze(rho(:,t,:));
        
        m = lambda_aux.*Xtilde;
        gamma = (1-lambda_aux.^2).* abs(Xtilde).^2;
        c = rho_aux.*Xtilde.^2;
        
        m_X = repmat(sum(m,2),[1 J]);
        gamma_X = repmat(sum(gamma,2),[1 J]);        
        c_X = repmat(sum(c,2),[1 J]);
        
        Xaux = repmat(X(:,t),[1 J]);
        
        % Get the MMSE estimator
        m_post = m + ( (gamma.*gamma_X - c.*conj(c_X)).*(Xaux-m_X) + (c.*gamma_X - gamma.*c_X).* conj(Xaux-m_X)    ) ./ (gamma_X.^2 - abs(c_X).^2+eps);
        
        Se(:,t,:) = m_post;
        
    end
end

end


% Inst. frequency and regions of influence
function [f_inf,f_centr,f_harm] = get_frequencies_qifft_frame(v)

v = v(:)';

%Central peaks
[~,f_centr] = findpeaks(max(v,10^-6),'MINPEAKHEIGHT',max(v)*10^-2);

Nfreq = length(f_centr);
f_harm = zeros(1,Nfreq);

if (Nfreq >0)
    % Quadratic interpolation of frequencies
    for ind = 1:Nfreq
        f = f_centr(ind);
        f_harm(ind) = qint(log10(v(f-1)),log10(v(f)),log10(v(f+1)))+f;
    end

    % Frequencies in Regions of influence
    f_inf = zeros(length(v),1);
    deb = 1;
    index_lim = zeros(1,Nfreq-1);
    
    for ind = 1:(Nfreq-1)
        f = f_centr(ind);
        fp = f_centr(ind+1);
        fin = floor((v(fp)*f+v(f)*fp)/(v(fp)+v(f)));
        f_inf(deb:fin) = f_harm(ind);
        deb = fin+1;
        index_lim(ind) = fin;
    end

    f_inf(deb:end) = f_harm(end)-1;
    
else
    f_inf = (1:length(v))'-1;
end

end

% Quadratic Interpolated FFT
function [p,b,a] = qint(ym1,y0,yp1)

p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
b = y0 - 0.25*(ym1-yp1)*p;
a = 0.5*(ym1 - 2*y0 + yp1);

end
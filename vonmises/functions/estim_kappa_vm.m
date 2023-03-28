function [kappa,ph_centered,phaux,f_centr] = estim_kappa_vm(X,Nfft,hop)

[F,T] = size(X);
phX = angle(X);

% Phase in the sinusoidal model
[f_inf,f_centr] = freq_unwrap(abs(X));
nu = f_inf / Nfft;
deltphi = 2*pi*hop*nu;
mu = [zeros(F,1) phX(:,1:end-1)] + deltphi;


% Get TF bins corresponding to peaks + first neighbors
f_centr_up = [zeros(1,T); f_centr(1:end-1,:)];
f_centr_down = [f_centr(2:end,:); zeros(1,T);];
f_centr = min(f_centr+f_centr_up+f_centr_down,1);


% Centered phases
phaux =  angle(exp(1i*(phX - mu)));
ph_centered = phaux(f_centr==1);
ph_centered(isnan(ph_centered)) = [];

% ML estimation of kappa
if ~isempty(ph_centered)
    kappa = fzero(@(x) besseli(1,x)/besseli(0,x)-mean(cos(ph_centered)),1);
else
    kappa =nan;
end


end




function [f_inf,f_centr,f_harm] = freq_unwrap(V)

[F,T] = size(V);

% Instantaneous frequencies
f_inf = zeros(F,T);
f_centr = zeros(F,T);
f_harm = cell(1,T);

for t=1:T
    [inf,centr,harm] = freq_influence(abs(V(:,t)));
    f_inf(:,t) = inf-1;
    if isempty(centr)
        f_centr(:,t) = 0;
    else
        c = centr;
        f_centr(c(:),t) = 1;
    end
    f_harm{t} = harm;
end


end



% PU: frequencies and regions of influence
function [f_inf,f_centr,f_harm] = freq_influence(v)

v = v(:)';

%Central peaks
%[~,f_centr] = findpeaks(v,'MINPEAKHEIGHT',0.01*max(v));
[~,f_centr] = findpeaks(max(v,10^-6),'MINPEAKHEIGHT',max(v)*10^-2);

Nfreq = length(f_centr);
f_harm = zeros(1,Nfreq);

if (Nfreq >0)
    % QIFFT
    for ind = 1:Nfreq
        f = f_centr(ind);
        f_harm(ind) = qint(log10(v(f-1)),log10(v(f)),log10(v(f+1)))+f;
    end

    % Regions of influence
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

    f_inf(deb:end) = f_harm(end);
    
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

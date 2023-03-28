%%% Plot von Mises densities

clear all; close all; clc;

Np=100;
phi = linspace(0,2*pi,Np);

mui=40;
mu = phi(mui);
kappa1 = 0;
kappa2 = 2;
kappa3 = 5;
kappa4 = 100;

p1 = exp(kappa1 * cos(phi - mu))/(2*pi*besseli(0,kappa1));
p2 = exp(kappa2 * cos(phi - mu))/(2*pi*besseli(0,kappa2));
p3 = exp(kappa3 * cos(phi - mu))/(2*pi*besseli(0,kappa3));
p4 = exp(kappa4 * cos(phi - mu))/(2*pi*besseli(0,kappa4));

plot(p1,'k+-'); hold on;

xpi = {'0';'';'';'';'';'2 \pi'};
ax=gca;ax.XTickLabel = xpi; 

plot(p2,'b.-'); hold on;
plot(p3,'ro-'); hold on;
plot(p4,'g*-'); hold on;

ha=legend('\kappa=0','\kappa=2','\kappa=5','\kappa=100'); set(ha,'FontSize',14);

plot([mui mui], [0 p4(40)],'k--'); hold on;
ha = text(mui,-0.1,'\mu'); set(ha,'FontSize',16);
ha=xlabel('\phi'); set(ha,'FontSize',16);
ha=ylabel('p(\phi|\mu,\kappa)'); set(ha,'FontSize',16);



%%% Plot densities of cos(phi) and sin(phi) where phi is von Mises
%clear all; close all; clc;

V = 20;
kappa = 2;
x = -V:0.0001:V;

% densities
densityx = @(x) exp(kappa*x/V) ./ (pi*besseli(0,kappa)*sqrt(V.^2-x.^2)  );
densityy = @(y) exp(kappa*sqrt(1-y.^2/V^2)) ./ (pi*besseli(0,kappa)*sqrt(V.^2-y.^2)  );

px = densityx(x);
py = densityy(x);

% simul
Nsample=100000;
wop = circ_vmrnd(0, kappa, Nsample);
cv = V*cos(wop); sv = V*sin(wop);

% plot cos and sin
figure;
subplot(1,2,1);
histogram(cv,1000,'Normalization','probability'); hold on;
plot(x(2:end-1),px(2:end-1)./px(end-1).*0.065); title('cos');
subplot(1,2,2);
histogram(sv,1000,'Normalization','probability'); hold on;
plot(x(2:end-1),py(2:end-1)./py(end-1).*0.009); title('sin');
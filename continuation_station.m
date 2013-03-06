clear; close all; clc;

N = 16;
Q = 4;
opts = optimset('display','off');

%% PARAMS

conv = 2.54/100;    % m/in

b       = [.4, .2]*conv; % base (beam1, beam2), m
ht      = [0.03, 0.03].*conv; % thickness (beam1, beam2), m
A       = b.*ht; % area (beam1, beam2), m^2
A1 = A(1); A2 = A(2);
L = 2.2*conv; % beam length
m = [.005, .005]; % magnet mass
M1 = m(1); M2 = m(2);
E = 7.5e10; % youngs of ceramic (estimated)
I = b.*ht.^3/12; % area moment
I1 = I(1); I2 = I(2); 
rho     = 1.2e3; % density, kg/m^3
d = 3/8*conv; % distance between two beams
a = .375/2 * conv; % magnet height / 2
V = [1,1].*(0.375^2)/4*pi* 2*a * conv^2; % magnet volume
V1 = V(1); V2 = V(2);
mu0 = 4*pi*10^(-7); % H/m, perm of free space
sigma_mag = (13200/10000)/mu0; % magnetization strenght per volume
sigma1 = sigma_mag; sigma2 = sigma_mag;

alpha1 =  M1 + 3/8 *rho*A1*L;
alpha2 =  M2 + 3/8 *rho*A2*L;
beta1 =  M1 + 33/140 *rho*A1*L;
beta2 =  M2 + 33/140 *rho*A2*L;
gamma1 = 3*E*I1/L^3;
gamma2 = 3*E*I2/L^3;

zeta = 0.05; % damping coeff
mu1 = 2*zeta*sqrt(gamma1*beta1);
mu2 = 2*zeta*sqrt(gamma2*beta2);

Amp = 15; % m/s^2, excitation amplitude

const1 = mu0.*V1.*V2.*sigma1.*sigma2;

%%

eqhandle = @(t,r,Omega) eqns_for_cont(t,r,Omega,Amp,L, a, d, const1, gamma1, gamma2, alpha1, alpha2, beta1, beta2, mu1, mu2);


om1 = 10*2*pi;
om2 = 10.05*2*pi;

[tt,D,~] = chebdif(N);

tt = (tt + 1)/2;  D = D*2;

D = [D(1:end-1,:);
    1, zeros(1,N-1), -1];

Dk = kron(eye(Q),D);

mask = [ones(N,1); 0];  mk = repmat(mask,Q,1);

resid = @(x,Om) Dk*x  - 3*(2*pi/Om)*mk.*reshape(eqhandle((2*pi/Om)*tt,reshape(x,N+1,Q),Om),Q*(N+1),1);


x0 = zeros(Q*(N+1),1);


x1 = fsolve(@(x) resid(x,om1),x0,opts);

X1 = reshape(x1,N+1,4);

% figure(1)
% subplot(211)
% plot(X1(:,1),X1(:,2))
% subplot(212)
% plot(X1(:,3),X1(:,4))

x2 = fsolve(@(x) resid(x,om2),x1,opts);

X2 = reshape(x2,N+1,4);


% subplot(211)
% hold on
% plot(X2(:,1),X2(:,2))
% subplot(212)
% hold on
% plot(X2(:,3),X2(:,4))

done = false;
xsi = 1;
ctr = 3;

ds = 1;

STORAGE = zeros(Q*(N+1),1000);

STORAGE(:,1) = x1; STORAGE(:,2) = x2;
PARAMS = zeros(1,1000);
PARAMS(1) = om1; PARAMS(2) = om2;

figure(2)
clf
hold on
view(3)
while ~done
    
    xpred = x2 + ds*(x2 - x1);
    ompred = om2 + ds*(om2 - om1);
    
    [OUT,~,~,INFO] = fsolve(@(IN) ...
        [resid(IN(1:end-1),IN(end));
        norm(IN - [x2;om2]) - ds],[xpred;ompred],opts);
    
    
    x1 = x2;
    om1 = om2;
    
    x2 = OUT(1:end-1); om2 = OUT(end);
    
    STORAGE(:,ctr) = x2; PARAMS(ctr) = om2;
    
    AMP = max(x2(1:N+1)) - min(x2(1:N+1));
    
    X2 = reshape(x2,N+1,4);
    
    plot3(om2/2/pi*ones(size(X2(:,1))),X2(:,1),X2(:,2),'k-')
    
    drawnow
                                    
    xsi = 3/INFO.iterations;
    
    if xsi > 2, xsi = 2; elseif xsi < 0.5, xsi = 0.5; end
    
    ds = ds*xsi;
    
    
    if om2 < 8*2*pi || om2 > 125*2*pi
        done = true;
    else
        ctr = ctr + 1;
    end
    
    
    
    
end

STORAGE = STORAGE(:,1:ctr); PARAMS = PARAMS(1:ctr);

%%

sample = randi(ctr);

Y = reshape(STORAGE(:,ctr),N+1,Q);

[TT,YY] = ode45(@(t,x) eqhandle(t,x',PARAMS(sample))',[0 40*2*pi/PARAMS(sample)],Y(1,:));

figure(3)
clf
hold on
plot(Y(:,1),Y(:,2),'k-',YY(end/2:end,1),YY(end/2:end,2),'r-')




clear; close all; clc;

addpath phd_utils
add_packages

%% Settings

N = 128; % N+1: number of cheby. nodes
opts = optimset('display','off'); % fsolve options

%% Specify parameters and dynamics

%-------------------------------
% System Parameters
%-------------------------------
mu = 0.1; % damping
d = 0.1; % height above platform/length of arm
wn2 = 1; % g/l
A = 3; % tilt amplitude (rads)

%-------------------------------
% Equation of motion
%-------------------------------

eqn = @(t,x,eta) [x(:,2), ... x1dot = x2
    -( mu*x(:,2) + ... x2dot = damping +
    ( -d*A*eta^2*sin(eta*t) + (A*eta*cos(eta*t)).^2.*cos(x(:,1)) ...
    -wn2*sin(eta*t) ).*sin(x(:,1)) )]; % restoring force terms

%-------------------------------
% Phase condition \ddot x(T) = 0
%-------------------------------

phasecond = @(t,x,eta) -( mu*x(2) + ...
    ( -d*A*eta^2*sin(eta*t(end)) + (A*eta*cos(eta*t(end))).^2.*cos(x(1)) ...
    -wn2*sin(eta*t(end)) ).*sin(x(1)) );

%% Set bifurcation parameter initial guesses

eta1 = 0.3; % bifurcation param, forcing frequency
eta2 = 0.31; % second value


%% Generate some simulations to look for periodic behavior

Ntests = 100;

ftest = linspace(0.2,1.8,Ntests);

figure(1)
clf

subplot(221)
clf
subplot(223)
clf
subplot(2,2,[2 4])
clf

for i = 1:Ntests
    
   dt = 2*pi/ftest(i)/N;
   
   [T,X] = ode45(@(t,x) eqn(t,x',ftest(i))',0:dt:150*2*pi/ftest(i),[0.1 0.1]);
   
   X = X(100*N:end,:);
   T = T(100*N:end,:);
    
   subplot(221)
   plot(T(30*N:end),X(30*N:end,1),'b-')
   subplot(223)
   plot(T(30*N:end),X(30*N:end,2),'b-')
   
   subplot(2,2,[2 4])
   plot(asin(sin(X(:,1))),X(:,2),'b-',asin(sin(X(1:N:end,1))),X(1:N:end,2),'go','markerfacecolor','g','markersize',4)
   title(['Frequency ' num2str(ftest(i))])
   eval(['export_fig ' 'fig' num2str(floor(100*ftest(i))) '.pdf -transparent']) 
   
end
    
    


%% Generate spectral differentiation matrix

[tt,D,~] = chebdif(N);

tt = (tt + 1)/2; % shift time from (-1,1) to (0,1)
D = 2*D; % and rescale differentiation matrix appropriately

D = [D(1:N,:);
    1, zeros(1,N-1), -1]; % remove last row of D, and replace with x(0) = x(N)

Dk = kron(eye(2),D); % inflate D to number of states

%% Constuct residual

BOUNDARY_CONDITION = [zeros(N,1);
    2*pi; % <------- x(0) + 2*pi = x(1) or v/v
    zeros(N+1,1)]; % change to x(0) + 2*pi = x(1)

mask = [ones(N,1); 0]; % to kill last row of the vector field so that only the boundary condition row of D is active
mk = repmat(mask,2,1); % inflate to # states

resid = @(x,T,eta) (Dk*x - ...
    T*mk.*reshape(eqn(T*tt,reshape(x,N+1,2),eta),2*(N+1),1) - ...
    BOUNDARY_CONDITION);

% residual = D*x - T*f(T*tau,x,eta),W where tau \in (0,1), t = T*tau,
% note that the BOUNDARY_CONDITION term doesn't influence the vector field
% at all, but is tacked on in a position corresponding to the masked rows
% of "eqn" so that the rows of D that were replaced with boundary
% conditions can be enforced.
%
% Also note that T (the period) is a free variable and the phase condition
% "phasecond" from above is added to keep the system square.

%% Get first two solutions

% initial guesses
x0 = ones(2*(N+1),1);
T0 = 2*pi/eta1;

% first solution, states and period
[s1,~,e1,~,jac] = fsolve(@(IN) ...
    [resid(IN(1:end-1),IN(end),eta1);
    phasecond(IN(end),[IN(N+1), IN(end-1)],eta1)],[x0;T0],opts);

x1 = s1(1:end-1); % first theta, dot theta solution
T1 = s1(end); % first period
disp('Condition number, attempt 1:')
cond(jac)

[s2,~,e2,~,jac] = fsolve(@(IN) ...
    [resid(IN(1:end-1),IN(end),eta2);
    phasecond(IN(end),[IN(N+1), IN(end-1)],eta2)],s1,opts);

if ~(e1 && e2)
    error('Initial steps failed')
end


% second soln, states and period
x2 = s2(1:end-1); % states
T2 = s2(end); % period
disp('Condition number, attempt 1:')
cond(jac)

% Note that the convention for "s1" and "s2" is that it holds the states
% and the period as [theta; dot theta; T]

%% plot first two solutions

figure(1)
clf
hold on
subplot(121)
hold on
xlabel('$\theta$','interpreter','latex','fontsize',14)
ylabel('$\dot \theta$','interpreter','latex','fontsize',14)

subplot(122)
hold on
xlabel('$\eta$','interpreter','latex','fontsize',14)
ylabel('$E$','interpreter','latex','fontsize',14)

%% Set loop parameters

done = false;
xsi = 1;
STORAGE = zeros(2*(N+1),500); % states only
PERIOD = zeros(1,500); % periods
PARAM = zeros(1,500); % parameters

STORAGE(:,1:2) = [x1, x2];
PERIOD(1:2) = [T1, T2];
PARAM(1:2) = [eta1; eta2];

ctr = 3;

%% Loop/continue

while ~done
    
    
    %--------------
    % PREDICTOR
    %--------------
    spred = s2 + xsi*(s2 - s1);
    etapred = eta2 + xsi*(eta2 - eta1);
    
    %--------------
    % CORRECTOR
    %--------------
    [OUT,fval,exitflag,info,jac] = fsolve(@(IN) ...
        [resid(IN(1:end-2),IN(end-1),IN(end)); % RESIDUAL
        phasecond(IN(end-1),[IN(N+1), IN(end-2)],IN(end)); % PHASE CONDITION
        (IN - [spred; etapred])'*([s2; eta2] - [s1; eta1])],... % DOT PRODUCT W OLD SOLN
        [spred;etapred],opts);
    
    %--------------
    % STORAGE
    %--------------
    
    s1 = s2; 
    eta1 = eta2;
    
    s2 = OUT(1:end-1);
    eta2 = OUT(end);
    x2 = s2(1:end-1);
    T2 = s2(end);
    
    STORAGE(:,ctr) = x2;
    PERIOD(ctr) = T2;
    PARAM(ctr) = eta2;
    
    %--------------
    % STEP CONTROL
    %--------------
    
    xsi = 3/info.iterations;
    
    if xsi > 2, xsi = 2; elseif xsi < 0.5, xsi = 0.5; end
    
    %--------------
    % PLOTS
    %--------------
    
    subplot(121)
    h = plot(x2(1:N+1),x2(N+2:end));
    subplot(122)
    plot(eta2,max(x2(N+2:end)) - min(x2(N+2:end)),'k.')
    drawnow
    
    %--------------
    % BREAK LOOP?
    %--------------
    
    if eta2 > 3 || eta2 < 0.2
        done = true;
        STORAGE = STORAGE(:,1:ctr);
        PERIOD = PERIOD(1:ctr);
        PARAM = PARAM(1:ctr);
    else
        delete(h)
    end
    
    ctr = ctr + 1;
    
        
end

%% Verify a random solution with ODE45

% choose a random trajectory
randindex = randi(ctr);


% ode45 it
[ts,xs] = ode45(@(t,x) eqn(t,x',PARAM(randindex))',...
    [0 PERIOD(randindex)],...
    [STORAGE(1,randindex); STORAGE(N+2,randindex)]);

% plot it: red is ode45, blue is collocation
figure(2)
clf

hold on
xlabel('$\theta$','interpreter','latex','fontsize',14)
ylabel('$\dot \theta$','interpreter','latex','fontsize',14)


plot(xs(:,1),xs(:,2),'r-o','markersize',4)
plot(STORAGE(1:N+1,randindex),STORAGE(N+2:end,randindex))
legend('ode45','collocation')

%% Simpson-Porco 2015
%% Parameters
% Microgrid's electrical parameters
clear variables
rng(5);

sample = 1;

if sample == 1
    N = [4 6 8]; % number of nodes in each inverter chain
else
    N = 4;
end
%%%%%%%%%%%%%%%%%%%%%%% Microgrid topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverter parameters %
Vdc = 650;
C = 25e-6;      % Filter Capacitance  [F]
Lf = 1.8e-3;    % Filter Inductance   [H]
Lo = 1.8e-3;    % Output Impedance    [H]
sys.Prat = [1.4 .7 .7 1.4]'; % Rated (max) Active Power   [kW]
sys.Qrat = [1.4 .7 .7 1.4]';   % Rated (max) Reactive Power [kVAr]

% Line impedances %
R = .5*[1 1 1];       % Resistance [Ohm]
L = [1 1 1]*5e-3; % Reactance  [H]

% Load impedances %
Rl = [1 1 1 1]'*0;       % Resistance [Ohm]
Ll = [1 1 1 1]'*0e-3;     % Reactance  [H] ????
loadconfig = [1 0 0 1]';
sys.Pl = sys.Prat*.5e3;  % Load Demanded Active Power   [kW]
sys.Ql = sys.Qrat*.5e3;  % Load Demanded Reactive Power [kVAr]

% Operation set points
sys.omegaast = 50*2*pi; % 50Hz
sys.Vast = 325.3;       % 230V RMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocating memory for variables
invChains(length(N)) = struct(...
    'y',0,...
    'P',0,...
    'Q',0,...
    'time',0,...
    'sys',0,...
    'ctrl',0,...
    'colvec', 0);

ii = 1;
for n = N
    sys.n = n;
    
    % Sampled line parameters
    if sample == 1
        R = datasample(R,sys.n-1);
        L = datasample(L,sys.n-1);
    end
    Zij = R + L*1i;   % Line Impedance
    Yij = 1./Zij;                   % Line Admittance
    gij = real(Yij);                % Line Conductance
    bij = imag(Yij);                % Line Susceptance
    sys.Bij = diag(bij,1) + diag(bij,-1);
    sys.Gij = diag(gij,1) + diag(gij,-1);
    
    % Sampled load parameters
    if sample == 1
        Rl = datasample(Rl,sys.n);
        Ll = datasample(Ll,sys.n);
    end
    sys.Zl = diag(Rl + Ll*1i); % Load impedance
    sys.Yl = 1./sys.Zl;                      % Load admittance
    sys.Yl(sys.Yl==Inf)=0;
    
    % Shunt admittance
    sys.Gii = diag(sum(sys.Gij,1)) + real(sys.Yl);
    sys.Bii = diag(sum(sys.Bij,1)) + imag(sys.Yl);
    
    % Communication
    % line
    sys.Link = diag(ones(1,sys.n-1),1) + diag(ones(1,sys.n-1),-1);
    
    % ring
    sys.Link = diag(ones(1,sys.n-1),1) + diag(ones(1,sys.n-1),-1) + ...
        diag(1,sys.n-1) + diag(1,-sys.n+1);
    
    % Sampled Rated Powers
    if sample == 1
        sys.Prat = datasample(sys.Prat,sys.n);
        sys.Qrat = sys.Prat;
    end
    
    % Load demanded power (??????)
    if sample == 1
        loadconfig = datasample(loadconfig,sys.n);
        sys.Pl = loadconfig.*datasample(sys.Pl,sys.n);
        sys.Ql = sys.Pl;
    else
        sys.Pl = loadconfig.*sys.Pl;
        sys.Ql = loadconfig.*sys.Ql;
    end
    
    sys.Plr = sys.Pl.*(1 - .3*rand(sys.n,1));
    sys.Qlr = sys.Ql.*(1 - .3*rand(sys.n,1));
    % checking
%     sys.Qast = (abs(imag(sys.Yl))*sys.Vast^2*[1 1 1 1]');
%     sys.Past = sys.Vast^2./Rl;
%     sys.Past(sys.Past==Inf)=0;
%     Qi = (abs(sys.Bii)*sys.Vast^2 - abs(sys.Bij)*sys.Vast^2)*[1 1 1 1]';
    
    % assert(all(sys.Qast-Qi<0.1),'erro')
    % return
    %% Simulation
    
    % Primary Controller
%     ctrl.etaP = [2.5 2.5 2.5 2.5]*1e-3; % P-omega droop coeff
    ctrl.etaP = 1./sys.Prat./100;
%     ctrl.etaQ = [1.5 1.5 1.5 1.5]*1e-3; % Q-V droop coeff
    ctrl.etaQ = 1./sys.Qrat./100;
    
    ctrl.k = [1.7 1.7 1.7 1.7]; % Int Frequency gain
    ctrl.kappa = [1 1 1 1]; % Int Voltage gain
    
    if sample == 1
        ctrl.k = datasample(ctrl.k,sys.n);
        ctrl.kappa = datasample(ctrl.kappa,sys.n);
    end
    
    % Secondary Controller
    % study 1a
    ctrl.Beta = 0;    % voltage control
    ctrl.h = 10;          % reactive power sharing [V]
    
    % study 1.b
    % ctrl.Beta = 2.2;    % voltage control
    % sys.h = 10;          % reactive power sharing [V]
    
    % study 1.c
    % ctrl.Beta = 1.2;    % voltage control
    % sys.h = 180; % [V]  % reactive power sharing [V]
    
    % study 1.d
    % ctrl.Beta = 4;    % voltage control
    % sys.h = 100; % [V]  % reactive power sharing [V]
    
    %%%%
    ctrl.Beta = 4.*[0 1 0*ones(1,n-2)]';
    ctrl.h = 2;
    
    ctrl.H = ctrl.h*sys.Link;
    
    
    %% Simulation
    %%%%%%%%%%%%%%%%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sys.simtime = 150; % [s]
    % flags
    flag.PriF = 1;
    flag.PriV = 1;
    flag.SecOn = 1;
    flag.SecOnT = 15;
    flag.LoadOff = 0;
    flag.LoadOffT = 30;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% Initial Conditions %%%% TODO: use previously generated dist
    delta0 = zeros(n,1);
    omega0 = sys.omegaast*ones(n,1) -.25 + .75*rand + [1; zeros(n-1,1)];%(n,1);
    V0 = (sys.Vast)*ones(n,1); - 1 + 2*rand(n,1);
    xi0 = zeros(n,1);
    zeta0 = zeros(n,1);
    state0 = [delta0  % volt angle (delta)
        omega0        % freq (omega)
        V0            % volt (V)
        xi0           % int freq (xi)
        zeta0         % int volt (zeta)
        ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% ODE %%%%%%%%%%%%%%%%%%
    [time,y] = ode45(@microgrid,[0 sys.simtime], state0, opts, sys, ctrl, flag);
    
    % u_final = zeros(length(time),sys.n*2);
    P_final = zeros(length(time),sys.n);
    Q_final = zeros(length(time),sys.n);
    for tt = time'
        ttt = find(ismember(time,tt));
        [dx, yy] = microgrid(tt,y(ttt,:)', sys, ctrl, flag);
        %     u_final(ttt,:) = yy.u;
        P_final(ttt,:) = yy.P;
        Q_final(ttt,:) = yy.Q;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % y = y(find(time>10):end,:);
    % P_final = P_final(find(time>10):end,:);
    % Q_final = Q_final(find(time>10):end,:);
    % time = time(time>10);
    % Populating result structure
    invChain = struct(...
        'y',y,...
        'P',P_final,...
        'Q',Q_final,...
        'time',time,...
        'sys',sys,...
        'ctrl',ctrl,...
        'colvec', colormap(hsv(n+1)));
    
    % accumulating final results
    invChains(ii) = invChain;
    ii = ii + 1;
end

%% Plotting
Nplot = 4;
y = invChains(Nplot).y;
n = invChains(Nplot).sys.n;
time = invChains(Nplot).time;
P_final = invChains(Nplot).P;
Q_final = invChains(Nplot).Q;

delta = y(:,1:n);
omega = y(:,n+1:2*n);
V = y(:,2*n+1:3*n);
xi = y(:,3*n+1:4*n);
zeta = y(:,4*n+1:5*n);

set(groot,'DefaultAxesFontSize', 16)
set(groot,'DefaultAxesLineWidth',2)
set(groot,'DefaultAxesBox','on')
set(groot,'DefaultLineLinewidth',2.5)
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextInterpreter','latex')

figure(1)
plot(time,wrapToPi(delta));
plot(time,(delta));
legend('DG1','DG2','DG3','DG4','location','best');
title('Volt angle $\delta$')
xlabel('$t$ [s]')

figure(2)
plot(time,omega);
legend('DG1','DG2','DG3','DG4','location','best');
title('Frequency $\omega$')
xlabel('$t$ [s]')
% axis([-Inf Inf 310 320])

figure(3)
plot(time,V);
legend('DG1','DG2','DG3','DG4','location','best');
title('Voltage $V$')
xlabel('$t$ [s]')
% axis([-Inf Inf 270 380])

figure(4)
plot(time,xi);
legend('DG1','DG2','DG3','DG4','location','best');
title('Secondary freq var $\xi$')
xlabel('$t$ [s]')

figure(5)
plot(time,zeta);
legend('DG1','DG2','DG3','DG4','location','best');
title('Secondary volt var $\zeta$')
xlabel('$t$ [s]')

figure(6)
plot(time,P_final);
legend('DG1','DG2','DG3','DG4','location','best');
title('Active Power $P$')
xlabel('$t$ [s]')

figure(7)
plot(time,Q_final);
legend('DG1','DG2','DG3','DG4','location','best');
title('Reactive Power $Q$')
xlabel('$t$ [s]')

distFig('screen','south')
%% Model
% Active and Reactive Powers
%
% $$P_i = \sum_{j=1}^nb_{ij}V_iV_j\sin(\delta_i-\delta_j) \\Q_i = b_{ii}V_i^2-\sum_{j=1}^nb_{ij}V_iV_j\cos(\delta_i-\delta_j)$$
%
% Frequency Controller
%
% $$\dot\omega_i = -(\omega_i-\omega^*)-\eta_{Pi}P_i+\xi_i \\k_i \dot\xi_i =
% -(\omega_i-\omega^*)-\sum_{j=1}^ng_{ij}(\xi_i-\xi_j)$$
%
% Voltage Controller
%
% $$\dot V_i = -(V_i-V^*) - \eta_{Qi}Q_i+\zeta_i \\\kappa_i\dot\zeta_i = -\beta_i(V_i-V^*)-\sum_{j=1}^nh_{ij}\left(\frac{Q_i}{Q_i^*}-\frac{Q_j}{Q_j^*}\right)$$
%
% Where $b_{ij} < 0\text{, and } b_{ii} = \sum_{j=1}^nb_{ij} + \hat b_{ii}$,
% and $\hat b_{ii} < 0$ is the shunt (load) susceptance .

function [dydt, y] = microgrid(t, state, sys, ctrl, flag)
n = sys.n;
I = eye(n);
Iv = ones(n,1);
O = zeros(n);
Ov = zeros(n,1);

Bij = sys.Bij;
Bii = sys.Bii;
Gij = sys.Gij;
Gii = sys.Gii;
omegaast = sys.omegaast*ones(n,1);
Vast = sys.Vast*ones(n,1);
Qrat = sys.Qrat;
Prat = sys.Prat;
% Qast = sys.Qast;
% Past = sys.Past;
Pl = sys.Pl;
Ql = sys.Ql;
Plr = sys.Plr;
Qlr = sys.Qlr;
Yl = sys.Yl;
Link = sys.Link;

Beta = ctrl.Beta;
h = ctrl.h;
H = ctrl.H;
K = diag(ones(1,n)./ctrl.k);
Kappa = diag(ones(1,n)./ctrl.kappa);

% states
% delta = wrapToPi(state(1:n));
delta = (state(1:n));
omega = state(n+1:2*n);
V = state(2*n+1:3*n);
xi = state(3*n+1:4*n);
zeta = state(4*n+1:5*n);

% Provides Vij in diagonal+1 and Vi^2 in main diagonal
VV = repmat(V,1,n).*repmat(V',n,1);

% theta_i - theta_i+1
delta_ij = repmat(delta,1,n) - repmat(delta',n,1);
% delta_ij(delta_ij<0.001)=0;


if t>flag.SecOnT
    flag.Sec = 1*flag.SecOn;
else
    flag.Sec = 0;
end
if t>flag.LoadOffT && flag.LoadOff
    Qlr(1) = 0; Plr(1) = 0;
    Ql(1) = 0; Pl(1) = 0;
%     VV = VV;
end 

% Powers used in controller
P = (abs(Bij).*sin(delta_ij)).*VV*Iv;
Q = (VV.*abs(Bii) - (abs(Bij).*cos(delta_ij)).*VV)*Iv;

% Powers output "measurement"
y.P = P;
y.Q = Q;
% y.P = (Gii.*VV - ...
%     ((B).*sin(delta_ij)+Gij.*cos(delta_ij)).*VV)*Iv;
% y.Q = (-VV.*(Bii) - ...
%     (-(B).*cos(delta_ij)+Gij.*sin(delta_ij)).*VV)*Iv;
y.P = Plr + (abs(Bij).*sin(delta_ij)).*VV*Iv;
y.Q = Qlr + (VV.*abs(Bii) - abs(Bij).*cos(delta_ij).*VV)*Iv;

% Y = imag(Yl) + Bij;
% y.Q = diag(V)*Y*V;

% Reactive Power Sharing
QDelta = repmat(y.Q./(Prat*.5e3),1,n) - repmat((y.Q./(Prat*.5e3))',n,1);

A = [...
    O I O O O;                      % deltadot = omega
    [O -I O flag.Sec*I O]*flag.PriF;    % omegadot = -omega + xi
    [O O -I O flag.Sec*I]*flag.PriV;    % Vdot = -V + zeta
    [O -K O -Link.*K O]*flag.Sec; % xidot = -omega/K -G*xi/K
    [O O -Beta.*Kappa O O]*flag.Sec];   % zetadot = -Beta*Kappa*V

ref = [omegaast; Vast; sys.Prat*.0e3; sys.Qrat*.0e3];
Aref = [O O O O;
    [I O I.*ctrl.etaP O]*flag.PriF;  % omegadot += omegaast + etaP*Past
    [O I O I.*ctrl.etaQ]*flag.PriV;  % Vdot += Vast + etaQ*Qast
    [K O O O]*flag.Sec;              % xidot += omegast/K
    [O Beta.*Kappa O O]*flag.Sec];    % zetadot += Vast*Beta/Kappa

PQ = [Ov
    -ctrl.etaP.*y.P*flag.PriF          % omegadot -= etaP*P
    -ctrl.etaQ.*y.Q*flag.PriV          % Vdot -= etaQ*Q
    %     -ctrl.etaP'.*Pc
    %     -ctrl.etaQ'.*Qc
    Ov
    %     Qshare
    -H./ctrl.kappa.*QDelta*Iv*flag.Sec;
    ];

dydt = A*state + ... % linear dynamics
    Aref*ref + ...     % ref (ast)
    PQ;           % nonlinear dynamics (powers)
end
%% Simpson-Porco 2015
%% Parameters
% Microgrid's electrical parameters
clearvars -except sol f
rng(1);

sample = 1; % if 1, creates length(N) chains with N inverters

if sample == 1
    N = [4 6 8]; % number of nodes in each inverter chain
else
    N = 4;
end

controllers = ["C1"];
controllers = ["C2"];
controllers = ["C1" "C2"];
% C1 - DSS Controller (Proposed)
% C2 - Simpson2015 (V-Q 'smart' tuning)


%%%%%%%%%%%%%%%%%%%%%%% Microgrid topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverter parameters % (for a 4 inverters chain)
Vdc = 650;
C = 25e-6;      % Filter Capacitance  [F]
Lf = 1.8e-3;    % Filter Inductance   [H]
Lo = 1.8e-3;    % Output Impedance    [H]
conf.Prat = 1.4;%[1.4 1.4 1.4 1.4]'; % Rated (max) Active Power   [kW]
conf.Qrat = .8;%[.8 .8 .8 .8]';   % Rated (max) Reactive Power [kVAr]

% Line impedances %
conf.R = [.8 .4 .7]';   % Resistance [Ohm]
conf.L = [3.6 1.8 1.9]'*1e-3; % Reactance  [H]

% Load impedances %
conf.Rl = [1 1 1 1]'*0;       % Resistance [Ohm]
conf.Ll = [1 1 1 1]'*0e-3;     % Reactance  [H] ????

conf.loadconfig = [1 1 0 1]'; % ??????
conf.invRating = ones(4,1);
conf.invRating(conf.loadconfig==0) = .5;
conf.invRating(conf.loadconfig==1) = 1;

conf.loaddist = (1 - .1*rand(4,1));

conf.Past = conf.invRating.*conf.Prat*.7e3;  % Load Active Power   [kW]
conf.Qast = conf.invRating.*conf.Qrat*.45e3;  % Load Reactive Power [kVAr]

% Operation set points
conf.omegaast = 50*2*pi; % 50Hz
conf.Vast = 325.3;       % 230V RMS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocating memory for variables
invChains(length(N)) = struct(...
    'y',0,...
    'P',0,...
    'Q',0,...
    'norms',0,...
    'time',0,...
    'sys',0,...
    'ctrl',0,...
    'colvec', 0);

ii = 1;
sys = conf;

%% DSS Controller Optimisation
% Initial parameters for the optimisation
Zij = conf.R + conf.L*1i;   % Line Impedance
Yij = 1./Zij;                   % Line Admittance
gij = real(Yij);                % Line Conductance
bij = imag(Yij);                % Line Susceptance
conf.Bij = diag(bij,1) + diag(bij,-1);
conf.Gij = diag(gij,1) + diag(gij,-1);
conf.Zl = diag(conf.Rl + conf.Ll*1i);  % Load impedance
conf.Yl = 1./conf.Zl;                 % Load admittance
conf.Yl(conf.Yl==Inf)=0;
conf.Bii = diag(sum(conf.Bij,1)) + imag(conf.Yl);

if exist('sol','var')
    optctrl = sol;
elseif ~exist('ctrl','var')
    try 
        Optimisation_SP;
    catch exp
        error(exp.message);
    end
    if constraints == 0
        clear ctrl sol
        error('constraints not satisfied during optimisation')
    end
    optctrl = sol
end
    
%% loop
for n = N
    sys.n = n;
    
    % Line parameters
    sys.R = [sys.R; datasample(conf.R,sys.n-length(sys.R)-1)];
    sys.L = [sys.L; datasample(conf.L,sys.n-length(sys.L)-1)];
    % Loads (not used, Pl and Ql instead) ?
    sys.Rl = [sys.Rl; datasample(conf.Rl,sys.n-length(sys.Rl))];
    sys.Ll = [sys.Ll; datasample(conf.Ll,sys.n-length(sys.Ll))];
    
    % Load configuration
    sys.loadconfig = [sys.loadconfig; ...
        datasample(conf.loadconfig,sys.n-length(sys.loadconfig))];
    sys.invRating(sys.loadconfig==0) = .5;
    sys.invRating(sys.loadconfig==1) = 1;
    
    % Rated Powers
    sys.Prat = [sys.Prat; ...
        sys.invRating(length(sys.Prat)+1:end).* ...
        datasample(conf.Prat,sys.n-length(sys.Prat))];
    sys.Qrat = [sys.Qrat; ...
        sys.invRating(length(sys.Qrat)+1:end).* ...
        datasample(conf.Qrat,sys.n-length(sys.Qrat))];
    
    % Power set points
    sys.Past = [sys.Past; ...
        sys.invRating(length(sys.Past)+1:end).* ...
        datasample(conf.Prat*.7e3,sys.n-length(sys.Past))];
    sys.Qast = [sys.Qast; ...
        sys.invRating(length(sys.Qast)+1:end).* ...
        datasample(conf.Qrat*.45e3,sys.n-length(sys.Qast))];
    
    % Load demanded Powers
    sys.Pl = sys.loadconfig.* sys.Past;
    sys.Ql = sys.loadconfig.* sys.Qast;
    sys.loaddist = [sys.loaddist; ...
        datasample(conf.loaddist,sys.n-length(sys.loaddist))];
    sys.Plr = sys.Pl.*sys.loaddist;
    sys.Qlr = sys.Ql.*sys.loaddist;
    
    Zij = sys.R + sys.L*1i;   % Line Impedance
    Yij = 1./Zij;                   % Line Admittance
    gij = real(Yij);                % Line Conductance
    bij = imag(Yij);                % Line Susceptance
    sys.Bij = diag(bij,1) + diag(bij,-1);
    sys.Gij = diag(gij,1) + diag(gij,-1);
    
    sys.Zl = diag(sys.Rl + sys.Ll*1i);  % Load impedance
    sys.Yl = 1./sys.Zl;                 % Load admittance
    sys.Yl(sys.Yl==Inf)=0;
    
    % Shunt admittance
    sys.Gii = diag(sum(sys.Gij,1)) + real(sys.Yl);
    sys.Bii = diag(sum(sys.Bij,1)) + imag(sys.Yl);
    %% %%%%%%%%%%%%%%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    sys.simtime = 150; % [s]
    % flags
    flag.PriF = 1;
    flag.PriV = 1;
    flag.SecOn = 1;
    flag.SecOnT = 0;
    flag.LoadOffOn = 1;
    flag.LoadOffT = 50;
    flag.LoadOnT = 100;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%%%%%% Controllers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for control = controllers
        sys.Link = diag(ones(1,sys.n-1),1) + diag(ones(1,sys.n-1),-1) + ...
                diag(1,sys.n-1) + diag(1,-sys.n+1);
        if control == "C1"            
            ctrl = optctrl;
            % Communication
            gim1 = ctrl.comm(1);
            gip1 = ctrl.comm(2);
            him1 = ctrl.comm(3);
            hip1 = ctrl.comm(4);
            
            gii = gip1 + gim1;
            hii = hip1 + him1;
            ctrl.G = diag([gip1 ; gii*ones(n-2,1); gim1])...
                - gip1*diag([1,ones(n-2,1)'],1)...
                - gim1*diag([ones(n-2,1)',1],-1);
            ctrl.H = diag([hip1 ; hii*ones(n-2,1); him1])...
                - hip1*diag([1,ones(n-2,1)'],1)...
                - him1*diag([ones(n-2,1)',1],-1);            
            
        elseif control == "C2"
            % Communication
            % line
%             sys.Link = diag(ones(1,sys.n-1),1) + diag(ones(1,sys.n-1),-1);
            % ring
            sys.Link = diag(ones(1,sys.n-1),1) + diag(ones(1,sys.n-1),-1) + ...
                diag(1,sys.n-1) + diag(1,-sys.n+1);
            %%% Primary Controller %%%
            ctrl.etaP = 1./sys.Prat./300;
            ctrl.etaP = 2.5e-3*ones(sys.n,1);
            ctrl.etaQ = 1./sys.Qrat./800;
            ctrl.etaQ = 2.5e-3*ones(sys.n,1);

            %%% Secondary Controller %%%
            ctrl.k = 1.7.*ones(1,sys.n); % Int Frequency gain
            ctrl.kappa = 1.*ones(1,sys.n); % Int Voltage gain

            % Simpson2015's Study 1.
            study = '1.c';
            switch study
                case '1.a' % reactive power sharing only
                    ctrl.Beta = 0;          % voltage control
                    ctrl.h = 50;            % reactive power sharing [V]
                case '1.b' % voltage control only
                    ctrl.Beta = 2.2;
                    ctrl.h = 0;
                case '1.c' % compromise V-Q
                    ctrl.Beta = 1.2;
                    ctrl.h = 18.0;
                case '1.d' % 'smart' tuning (leader-follower)
                    ctrl.Beta = zeros(n,1);
                    ctrl.Beta(2) = 4;       % leader
                    ctrl.h = 100;
                    sys.Link(2,:) = zeros(1,n);
                otherwise
                    ctrl.Beta = 0;
                    ctrl.h = 10;
            end
            ctrl.H = ctrl.h*sys.Link;
        end        
        ctrl.control = control;    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %% %% Initial Conditions
        delta0 = zeros(n,1);
        omega0 = sys.omegaast*ones(n,1); + [1; zeros(n-1,1)];
        V0 = (sys.Vast)*ones(n,1);% - 1 + 2*rand(n,1);
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
        [time,y] = ode45(@microgrid,[0 sys.simtime], state0, opts, sys,...
            ctrl, flag);
        
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
        
        %         return
        
        %% % --- Calculating Norms --- %%%
        l2norma = zeros(length(y),n);
        z = reshape(y, length(y),n,(size(y,2))/n);
        for idx_t = time' % every time step
            tt = find(ismember(time,idx_t));
            for idx_a = 1:n % every agent
                ztemp = z(tt,idx_a,:);
                ztemp = reshape(ztemp,1,(size(y,2))/n);
                
                % angle error
%                 if wrap
                ztemp(1) = 0.*wrapToPi(ztemp(1));
                ztemp(2) = ztemp(2) - conf.omegaast;
                ztemp(3) = ztemp(3) - conf.Vast;
%                 end
                
                % l2 norm for each agent
                l2norma(tt,idx_a) = norm(ztemp);
            end
        end
        % max l2 norm among agents
        l2normt = max(l2norma,[],2);                
        % supreme max l2 
        [linfnorm, ~] = max(max(l2norma));
        
        norms.l2 = l2normt;
        norms.linf = linfnorm;
        %%% --- ----------------- --- %%%
        
        % y = y(find(time>10):end,:);
        % P_final = P_final(find(time>10):end,:);
        % Q_final = Q_final(find(time>10):end,:);
        % time = time(time>10);
                
        % - Populating result structure - %
        invChain = struct(...
            'y',y,...
            'P',P_final,...
            'Q',Q_final,...
            'norms',norms,...
            'time',time,...
            'sys',sys,...
            'ctrl',ctrl,...
            'colvec', colormap(hsv(n+1)));
        
        invChains(ii) = invChain;
        ii = ii + 1;
    end
end

%% Plotting
Nplot = 2;
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

% figure(1)
% plot(time,wrapToPi(delta));
% plot(time,(delta));
% legend('DG1','DG2','DG3','DG4','location','best');
% title('Volt angle $\delta$')
% xlabel('$t$ [s]')

figure(2);
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
% axis([-Inf Inf 100 400])

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

% if ~exist('f','var')
    distFig('screen','east')
% end


%%
f = figure(8);
clf
Ns = length(N)*length(controllers);
% leg = string(1:Ns);
for x = 1:Ns
    if mod(x,2)==1 % Controller C1
%         leg(x) = strcat("C1 N=",);
        colvec = [1 0 0];
    else % Controller C2
%         leg(x) = strcat("C2 N=",string(mod(x,length(controllers))));
        colvec = [0 1 0];
    end
    
    plot(invChains(x).time,invChains(x).norms.l2,'color',colvec.*x/Ns)    
    leg(x) = strcat(invChains(x).ctrl.control," N=",...
        num2str(invChains(x).sys.n));
    hold on
end
legend(leg)
%%
%% ODE Func
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
% Qrat = sys.Qrat;
% Prat = sys.Prat;
Qast = sys.Qast;
Past = sys.Past;
Pl = sys.Pl;
Ql = sys.Ql;
Plr = sys.Plr;
Qlr = sys.Qlr;
Yl = sys.Yl;

% states
% delta = wrapToPi(state(1:n));
delta = (state(1:n));
% omega = state(n+1:2*n);
V = state(2*n+1:3*n);
% xi = state(3*n+1:4*n);
% zeta = state(4*n+1:5*n);

% Provides Vij in diagonal+-1 and Vi^2 in main diagonal
VV = repmat(V,1,n).*repmat(V',n,1);

delta_ij = repmat(delta,1,n) - repmat(delta',n,1);
% delta_ij(delta_ij<0.001)=0;


if t>flag.SecOnT
    flag.Sec = 1*flag.SecOn;
else
    flag.Sec = 0;
end
if t>flag.LoadOffT && flag.LoadOffOn && t<flag.LoadOnT
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
y.P = Plr + (abs(Bij).*sin(delta_ij)).*VV*Iv;
y.Q = Qlr + (VV.*abs(Bii) - abs(Bij).*cos(delta_ij).*VV)*Iv;

% Y = imag(Yl) + Bij;
% y.Q = diag(V)*Y*V;

% Reactive Power Sharing
QDelta = repmat(y.Q./(Qast),1,n) - repmat((y.Q./(Qast))',n,1);

ref = [omegaast; Vast; Past; Qast];


if ctrl.control == "C1" %TODO DSS Controller      
    tau = ctrl.ctrlgain(1);
    np = ctrl.ctrlgain(2);
    nq = ctrl.ctrlgain(3);
    kxi = ctrl.ctrlgain(4);
    kzeta = ctrl.ctrlgain(5);
%     tau = 1;
%     np = 2.5e-3;
%     nq = 2.5e-3;
%     kxi = 1.7;    
%     kzeta = 1;
    
    T = diag(ones(1,n)./tau);    
    Np = diag(np);
    Nq = diag(nq);
    K = diag(ones(1,n)./kxi);
    Kappa = diag(ones(1,n)./kzeta);
    
    L = diag(ones(1,n).*ctrl.ctrlxtra(1));
    M = diag(ones(1,n).*ctrl.ctrlxtra(2));
    Og = diag(ones(1,n).*ctrl.ctrlxtra(3));
    X = diag(ones(1,n).*ctrl.ctrlxtra(4));
    Y = diag(ones(1,n).*ctrl.ctrlxtra(5));
    
    A = [...
        O I O O O;                        % deltadot = omega
        [O -T O flag.Sec*T*Y O]*flag.PriF;  % omegadot = -omega + xi
        [O O -T*Og O flag.Sec*T*M]*flag.PriV;  % Vdot = -V + zeta
        [O -K*X O -K*ctrl.G O]*flag.Sec;     % xidot = -omega/K -G*xi/K
        [O O -Kappa*L O -Kappa*ctrl.H]*flag.Sec]; % zetadot = -Beta*Kappa*V
    
    Aref = [O O O O;
        [T O T*Np O]*flag.PriF;  % omegadot += omegaast + etaP*Past
        [O T*Og O T*Nq]*flag.PriV;  % Vdot += Vast + etaQ*Qast
        [K*X O O O]*flag.Sec;              % xidot += omegast/K
        [O Kappa*L O O]*flag.Sec];   % zetadot += Vast*Beta/Kappa
    
    PQ = [Ov
        -Np.*T*y.P*flag.PriF          % omegadot -= etaP*P
        -Nq.*T*y.Q*flag.PriV          % Vdot -= etaQ*Q
        Ov
        Ov%-H./ctrl.kappa.*QDelta*Iv*flag.Sec;
        ];
elseif ctrl.control == "C2" 
    K = diag(ones(1,n)./ctrl.k);
    Kappa = diag(ones(1,n)./ctrl.kappa);
    A = [...
        O I O O O;                      % deltadot = omega
        [O -I O flag.Sec*I O]*flag.PriF;    % omegadot = -omega + xi
        [O O -I O flag.Sec*I]*flag.PriV;    % Vdot = -V + zeta
        [O -K O -sys.Link.*K O]*flag.Sec; % xidot = -omega/K -G*xi/K
        [O O -ctrl.Beta.*Kappa O O]*flag.Sec];   % zetadot = -Beta*Kappa*V
    
    Aref = [O O O O;
        [I O I.*ctrl.etaP O]*flag.PriF;  % omegadot += omegaast + etaP*Past
        [O I O I.*ctrl.etaQ]*flag.PriV;  % Vdot += Vast + etaQ*Qast
        [K O O O]*flag.Sec;              % xidot += omegast/K
        [O ctrl.Beta.*Kappa O O]*flag.Sec];    % zetadot += Vast*Beta/Kappa
    
    PQ = [Ov
        -ctrl.etaP.*y.P*flag.PriF          % omegadot -= etaP*P
        -ctrl.etaQ.*y.Q*flag.PriV          % Vdot -= etaQ*Q
        Ov
        -ctrl.H./ctrl.kappa.*QDelta*Iv*flag.Sec;
        ];
end

dydt = A*state + ... % linear dynamics
    Aref*ref + ...     % ref (ast)
    PQ;           % nonlinear dynamics (powers)
end
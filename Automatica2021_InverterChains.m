%% Simpson-Porco 2015
home
%% Parameters
% Microgrid's electrical parameters
clearvars -except sol f
rng(1);

SStest = 1; % if 1, creates length(N) chains with N inverters
savefigs = 0;
voltctrl = 0;
printfigs = 1;
calcnorms = 1;

if SStest == 1
%     N = 4:32:100; % number of nodes in each inverter chain
    N = 4:6:40; % number of nodes in each inverter chain
else
    N = 60;
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
% conf.R = [.8 .8 .8]';   % Resistance [Ohm]
% conf.L = [3.6 3.6 3.6]'*1e-3; % Reactance  [H]

% Load impedances %
conf.Rl = [1 1 1 1]'*0;       % Resistance [Ohm]
conf.Ll = [1 1 1 1]'*0e-3;     % Reactance  [H] ????

conf.loadconfig = [1 1 1 1]'; % ??????
conf.invRating = ones(4,1);
conf.invRating(conf.loadconfig==0) = 0.5;
conf.invRating(conf.loadconfig==1) = 1;

conf.loaddist = (1 - .0*rand(4,1));

% Load Active [kW] and reactive [kVAr] Power   
conf.Past = conf.invRating.*conf.Prat.*(.9-.0*rand(4,1))*1e3;  
conf.Qast = conf.invRating.*conf.Qrat.*(.8-.0*rand(4,1))*1e3;

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
        
        if voltctrl
            Optimisation_SP;
        else
            Optimisation_Vcte;
        end
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
fprintf('Simulating %d strings...\n',length(N)*length(controllers))
for n = N
    sys.n = n;
    fprintf('\t %d/%d \n',n,N(end))
    
    % Line parameters
    zselect = datasample(1:3,1);
    sys.R = [sys.R; repmat(conf.R(zselect),sys.n-length(sys.R)-1,1)];
    sys.L = [sys.L; repmat(conf.L(zselect),sys.n-length(sys.L)-1,1)];
    % Loads (not used, Pl and Ql instead) ?
    sys.Rl = [sys.Rl; datasample(conf.Rl,sys.n-length(sys.Rl))];
    sys.Ll = [sys.Ll; datasample(conf.Ll,sys.n-length(sys.Ll))];
    
    % Load configuration
    sys.loadconfig = [sys.loadconfig; ...
%   repmat(datasample(conf.loadconfig,1),sys.n-length(sys.loadconfig),1)];
        repmat(1,sys.n-length(sys.loadconfig),1)];
    sys.invRating(sys.loadconfig==0) = 0.5;
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
        datasample(conf.Past,sys.n-length(sys.Past))];
    sys.Qast = [sys.Qast; ...
        sys.invRating(length(sys.Qast)+1:end).* ...
        datasample(conf.Qast,sys.n-length(sys.Qast))];
    
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
    sys.simtime = 8000; % [s]
    % flags
    flag.PriF = 1;
    flag.PriV = voltctrl;
    flag.SecOn = 1;
    flag.SecOnT = 0;
    flag.Dist = 1;
    flag.ic = 0;
    flag.DistOnT = 0;
    flag.DistOffT = 4000;
    opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%%%%%% Controllers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for control = controllers
%       sys.Link = diag(ones(1,sys.n-1),1) + diag(ones(1,sys.n-1),-1) + ...
%               diag(1,sys.n-1) + diag(1,-sys.n+1);
        if control == "C1"            
            ctrl = optctrl;
            if ~voltctrl
                ctrl.ctrlgain(1) = optctrl.tau;     % tau
                ctrl.ctrlgain(2) = optctrl.etai;    % etaP
                ctrl.ctrlgain(3) = 1;               % etaQ
                ctrl.ctrlgain(4) = optctrl.ki;      % kxi
                ctrl.ctrlgain(5) = 1;               % kzeta
                
                ctrl.ctrlxtra(1) = 1;%optctrl.Oo;
                ctrl.ctrlxtra(2) = 1;%optctrl.Oxi;
                ctrl.ctrlxtra(3) = 1;
                ctrl.ctrlxtra(4) = 0;%optctrl.Xid;
                ctrl.ctrlxtra(5) = 1;%optctrl.Xio;
                ctrl.ctrlxtra(6) = 1;
                
                gim1 = optctrl.gim1;
                gip1 = optctrl.gip1;                
                him1 = 0;                
                hip1 = 0;
            else                
                gim1 = ctrl.comm(1);
                gip1 = ctrl.comm(2);
                him1 = ctrl.comm(3);
                hip1 = ctrl.comm(4);
            end
            % Communication
            
            gii = gip1 + gim1;
            hii = hip1 + him1;
            ctrl.G = (diag([gip1 ; gii*ones(n-2,1); gim1]) ...
                - gip1*diag([1,ones(n-2,1)'],1)...
                - gim1*diag([ones(n-2,1)',1],-1));...
%             - gim1.*diag(1,n-1)...   % those lines implement ring 
%             - gip1.*diag(1,-n+1));   % connection if first line is zero
%           ctrl.G = ctrl.G - diag(sum(ctrl.G,2)); 
            ctrl.H = diag([hip1 ; hii*ones(n-2,1); him1])...
                - hip1*diag([1,ones(n-2,1)'],1)...
                - him1*diag([ones(n-2,1)',1],-1);            
            
        elseif control == "C2"
            % Communication
            % line
%           sys.Link = diag(ones(1,sys.n-1),1) + diag(ones(1,sys.n-1),-1);
            
            sys.Link = 1*(diag([1 ; 2*ones(n-2,1); 1]) ...
                - 1*diag([1,ones(n-2,1)'],1)...
                - 1*diag([ones(n-2,1)',1],-1));
            
            % ring
%           sys.Link = diag(ones(1,sys.n-1),1) + ...
%            diag(ones(1,sys.n-1),-1) + diag(1,sys.n-1) + diag(1,-sys.n+1);
%           sys.Link = -sys.Link + diag(sum(sys.Link,2));
            %% Primary Freq Controller %%%
            ctrl.etaP = 6./sys.Past;
%             ctrl.etaP = 6e-4*ones(sys.n,1);
            ctrl.etaQ = 1./sys.Qast./800;
%             ctrl.etaQ = 2.5e-3*ones(sys.n,1);

            %% Secondary Freq Controller %%%
            ctrl.k = .8.*ones(1,sys.n); % Int Frequency gain
            ctrl.kappa = 1.*ones(1,sys.n); % Int Voltage gain

            %% Simpson2015's Study 1. - Voltage control
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
                    ctrl.Beta(2) = 10;       % leader
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
        omega0 = sys.omegaast*ones(n,1) + [.5; 0*rand(n-1,1)]*flag.ic;
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
        sys.dist = 1 - .1*rand(n,1);
        [time,y] = ode45(@microgrid,[0 sys.simtime], state0, opts, sys,...
            ctrl, flag);
        
        % u_final = zeros(length(time),sys.n*2);
        P_final = zeros(length(time),sys.n);
        Q_final = zeros(length(time),sys.n);
        for tt = time'
            idx_t = find(ismember(time,tt));
            [dx, yy] = microgrid(tt,y(idx_t,:)', sys, ctrl, flag);
            %     u_final(ttt,:) = yy.u;
            P_final(idx_t,:) = yy.P;
            Q_final(idx_t,:) = yy.Q;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % - Populating result structure - %
        invChain.y = y;
        invChain.P = P_final;
        invChain.Q = Q_final;
        invChain.norms=[];
        invChain.time = time;
        invChain.sys = sys;
        invChain.ctrl = ctrl;
        invChain.colvec = colormap(hsv(n));
        
        invChains(ii) = invChain;        
        clear invChain
        
%         clear y P_final Q_final norms time ctrl colormap
        
        ii = ii + 1;
    end
end

%% Plotting
fprintf('Plotting...\n')
try
    Nplot = 40;
    Controller = "C2";
    idx_N = find(ismember(N,Nplot));
catch
    Nplot = N(end);
    Controller = "C2";
    idx_N = find(ismember(N,Nplot));
end

if Controller=="C2"
    idx_N = idx_N+1 + 1*(idx_N-1);
else
    idx_N = idx_N + 1*(idx_N-1);
end
    
y = invChains(idx_N).y;
n = invChains(idx_N).sys.n;
time = invChains(idx_N).time;
P_final = invChains(idx_N).P;
Q_final = invChains(idx_N).Q;
colvect = invChains(idx_N).colvec;

delta = (y(:,1:n));
omegastart = (conf.omegaast*time);
omega = y(:,n+1:2*n);
V = y(:,2*n+1:3*n);
xi = y(:,3*n+1:4*n);
zeta = y(:,4*n+1:5*n);

[~, idx_DistOn] = min(abs(time - flag.DistOnT));
[~, idx_DistOff] = min(abs(time - flag.DistOffT));

delta_s = zeros(size(delta));   

delta_s(:) = repmat(delta(end-1,:) - omegastart(end-1),length(delta),1); 
delta_s(idx_DistOn:idx_DistOff,:) = repmat(delta(idx_DistOff-1,:) ...
    - omegastart(idx_DistOff-1),idx_DistOff-idx_DistOn+1,1);

set(groot,'DefaultAxesFontSize', 12)
set(groot,'DefaultAxesLineWidth',1.5)
set(groot,'DefaultAxesBox','on')
set(groot,'DefaultLineLinewidth',2.5)
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultAxesColorOrder',colvect)
set(groot,'DefaultAxesColorMap',colvect)
xlimm = [-Inf Inf];
% xlimm = [0 200];

figure(1); clf;
hold on
plot(time,180./pi.*(delta - omegastart - delta_s)); 
xlabel('$t$ [s]')
ylabel('$\delta_{i} - \delta^*$')
colorbar(gca,'Ticks',(0:1/(Nplot-1):1),'TickLabels',string(1:Nplot),...
    'LineWidth',1,'Direction','reverse','TickLength',0); colormap(colvect)
axis([xlimm -Inf Inf])

figure(2); clf;
plot(time,omega);
% legend('DG1','DG2','DG3','DG4','location','best');
xlabel('Time [s]')
ylabel('Frequency [Hz]')
% axis([-Inf Inf 310 320])
axis([xlimm -Inf Inf])
colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
    'TickLabels',string(1:Nplot),...
    'LineWidth',1,...
    'Direction','reverse',...
    'TickLength',0); 
colormap(colvect)

figure(3); clf;
plot(time,xi);
% legend('DG1','DG2','DG3','DG4','location','best');
% title('Secondary freq var $\xi$')
xlabel('Time [s]')
ylabel('$\xi$')
axis([xlimm -Inf Inf])
colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
    'TickLabels',string(1:Nplot),...
    'LineWidth',1,...
    'Direction','reverse',...
    'TickLength',0); 
colormap(colvect)

% figure(3)
% plot(time,V);
% % legend('DG1','DG2','DG3','DG4','location','best');
% % title('Voltage $V$')
% xlabel('Time [s]')
% ylabel('Voltage [V]')
% % axis([-Inf Inf 100 400])
colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
    'TickLabels',string(1:Nplot),...
    'LineWidth',1,...
    'Direction','reverse',...
    'TickLength',0); 
colormap(colvect)

% figure(4)
% plot(time,zeta);
% % legend('DG1','DG2','DG3','DG4','location','best');
% % title('Secondary volt var $\zeta$')
% xlabel('Time [s]')
% ylabel('$\zeta$')
colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
    'TickLabels',string(1:Nplot),...
    'LineWidth',1,...
    'Direction','reverse',...
    'TickLength',0); 
colormap(colvect)

figure(4); clf;
plot(time,P_final); 
% legend('DG1','DG2','DG3','DG4','location','best');
% title('Active Power $P$')
xlabel('Time [s]')
ylabel('$P$ [W]')
colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
    'TickLabels',string(1:Nplot),...
    'LineWidth',1,...
    'Direction','reverse',...
    'TickLength',0); 
colormap(colvect)

% figure(6)
% plot(time,Q_final);
% % legend('DG1','DG2','DG3','DG4','location','best');
% % title('Reactive Power $Q$')
% xlabel('Time [s]')
% ylabel('$Q$ [VAr]')
colorbar(gca,'Ticks',(0:1/(Nplot-1):1),...
    'TickLabels',string(1:Nplot),...
    'LineWidth',1,...
    'Direction','reverse',...
    'TickLength',0); 
colormap(colvect)

%% % --- Calculating Norms --- %%%
if calcnorms == 1
    fprintf('Calculating norms...\n')
for idx = 1:length(invChains)
    fprintf('\t%d/%d\n',idx,length(N)*length(controllers))
    n = invChains(idx).sys.n;
    y = invChains(idx).y;
%     y(:,1:n) = wrapTo2Pi(y(:,1:n));
    time = invChains(idx).time;    
    delta = y(:,1:n);
%     delta = wrapTo2Pi(y(:,1:n));
    omegastart = conf.omegaast*time;
%     omegastart = wrapTo2Pi(conf.omegaast*time);
    xi = y(:,3*n+1:4*n);
    l2norma = zeros(size(y,1),n);                   % ok
    z = reshape(y, size(y,1),n,(size(y,2))/n);      % ok
    [~, idx_DistOff] = min(abs(time - flag.DistOffT));
    for tt = time' % every time step
        idx_t = find(ismember(time,tt));
        for idx_a = 1:n % every agent
            ztemp = z(idx_t,idx_a,:);                  % ok
            ztemp = reshape(ztemp,1,(size(y,2))/n);

            ztemp(2) = ztemp(2) - conf.omegaast;
            ztemp(3) = ztemp(3) - conf.Vast;            

            % disturbance on?                
            if flag.Dist && tt > flag.DistOnT && tt < flag.DistOffT 
%                 [~, idx_DistOn] = min(abs(time - flag.DistOnT));
                ztemp(1) = ztemp(1) - omegastart(idx_t) - ...
                  (delta(idx_DistOff-1,idx_a) - omegastart(idx_DistOff-1));
                ztemp(1) = ztemp(1).*180/pi;
                ztemp(4) = ztemp(4) - xi(idx_DistOff-1,idx_a); 
            else
                ztemp(1) = ztemp(1) - omegastart(idx_t) - ...
                    (delta(end,idx_a) - omegastart(end));
                ztemp(1) = ztemp(1).*180/pi;
                ztemp(4) = ztemp(4) - xi(end,idx_a);
            end
            
            ztemp(5) = 0;
            % l2 norm for each agent
            l2norma(idx_t,idx_a) = norm(ztemp);
        end
    end
    % max l2 norm among agents
    l2normt = max(l2norma,[],2);                
    % supreme max l2 
    [linfnorm, ~] = max(max(l2norma));

    invChains(idx).norms.l2 = l2normt;
    invChains(idx).norms.linf = linfnorm;
%         norms=0;
    %%% --- ----------------- --- %%%

    % y = y(find(time>10):end,:);
    % P_final = P_final(find(time>10):end,:);
    % Q_final = Q_final(find(time>10):end,:);
    % time = time(time>10);
end
end
%%
f = figure(5);
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
    if invChains(x).ctrl.control == "C1"
        plot(invChains(x).time,invChains(x).norms.l2,...
            '-','color',colvec.*x/Ns)    
    elseif invChains(x).ctrl.control == "C2"
        plot(invChains(x).time,invChains(x).norms.l2,...
            '--','color',colvec.*x/Ns)    
    end
    leg(x) = strcat("N=",num2str(invChains(x).sys.n)," ",...
        invChains(x).ctrl.control);
    hold on
end

xlabel('Time [s]')
ylabel('$\sup_i\left|x_i-x_i^*\right|_2$')
% legend(leg)

%%
norms_ours_v = [invChains(1:2:end).norms];
norms_naiv_v = [invChains(2:2:end).norms];
figure(6); clf
stem(N-.1,[norms_ours_v.linf],'--*','r','linewidth',2); hold on
stem(N+.1,[norms_naiv_v.linf],'--*','g','linewidth',2)
legend('DSS Controller','Standard Controller','location','best')
xlabel('Number of agents $n$')
ylabel('$\sup_i\left|\left|x_i-x_i^*\right|\right|_\infty$')
xlim([N(1)-2 N(end)+2])
xticks(N)

%% Printing/Saving figures
% clc
if savefigs == 1
    
    figs = findobj(0, 'type', 'figure');
    figs = figs(end:-1:1);
    datetime_char = string(datetime('now','format','yyyyMMMdd_hhmm'));
    for k=1:length(figs)
        % print each figure in figs to a separate .eps file
        print(figs(k), '-depsc', ...
            sprintf("/figs/N%d_%s_%d.eps", Nplot, datetime_char, k))
        %TODO SAVE .fig
    end
end
%%
% distFig('screen','east')
%% ODE Func
function [dydt, y] = microgrid(t, state, sys, ctrl, flag)

persistent tlast
if ~isempty(tlast)
    tlast = 0;
end

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

if t>flag.SecOnT
    flag.Sec = 1*flag.SecOn;
else
    flag.Sec = 0;
end

% Powers used in controller
P = (abs(Bij).*sin(delta_ij)).*VV*Iv;
Q = (VV.*abs(Bii) - (abs(Bij).*cos(delta_ij)).*VV)*Iv;

if t>flag.DistOnT && flag.Dist && t<flag.DistOffT
    Plr(1) = Plr(1).*0.8;%sys.dist(1:4:n);
%     Plr(ceil(end/2)) = Plr(ceil(end/2)).*0.6;%sys.dist(1:4:n);
    dist = zeros(n,1);
else
    dist = zeros(n,1);
end

% Powers output "measurement"
y.P = P;
y.Q = Q;
y.P = Plr + P;
y.Q = Qlr + Q;

% Y = imag(Yl) + Bij;
% y.Q = diag(V)*Y*V;

% Reactive Power Sharing
QDelta = repmat(y.Q./(Qast),1,n) - repmat((y.Q./(Qast))',n,1);

ref = [omegaast; Vast; Past; Qast];

if ctrl.control == "C1" 
    tau =   ctrl.ctrlgain(1);  % + makes the system faster
    np =    ctrl.ctrlgain(2);  % + improves power sharing
    nq =    ctrl.ctrlgain(3);
    kxi =   ctrl.ctrlgain(4);  % + makes system slower
    kzeta = ctrl.ctrlgain(5);    
    
    Oo      = diag(ones(1,n).*ctrl.ctrlxtra(1));    % omega -> omega
    Oxi     = diag(ones(1,n).*ctrl.ctrlxtra(2));    % omega -> xi
    Vzeta   = diag(ones(1,n).*ctrl.ctrlxtra(3));    % V -> zeta    
    Xid     = diag(ones(1,n).*ctrl.ctrlxtra(4));    % xi -> delta
    Xio     = diag(ones(1,n).*ctrl.ctrlxtra(5));    % xi -> omega
    ZetaV   = diag(ones(1,n).*ctrl.ctrlxtra(6));    % zeta -> V
    
    
    T = diag(ones(1,n)./tau);    
    Np = diag(ones(1,n)*np);
    Nq = diag(nq);
    K = diag(ones(1,n)./kxi);
    Kappa = diag(ones(1,n)./kzeta);
    
    A = [...
        O I O O O;                        
        [O -T*Oo O flag.Sec*T*Oxi O]*flag.PriF;  
        [O O -T O flag.Sec*T*Vzeta]*flag.PriV;  
        [O -K*Xio O -K*ctrl.G O]*flag.Sec;     
        [O O -Kappa*ZetaV O -Kappa*ctrl.H]*flag.Sec*flag.PriV]; 
    
    Aref = [O O O O;
        [T*Oo O T*Np O]*flag.PriF;  
        [O T O T*Nq]*flag.PriV;  
        [K*Xio O O O]*flag.Sec;              
        [O Kappa*ZetaV O O]*flag.Sec*flag.PriV];
    
    PQ = [Ov
        -Np.*T*y.P*flag.PriF          % omegadot -= etaP*P
        -Nq.*T*y.Q*flag.PriV          % Vdot -= etaQ*Q
        Ov
%         -Np*K*y.P*flag.Sec
        Ov%-H./ctrl.kappa.*QDelta*Iv*flag.Sec;
        ];
    
    Dist = [Ov; Np*dist; Ov; Ov; Ov];
elseif ctrl.control == "C2" 
    if t>30
     a = 1;
    end
%     
    K = diag(ones(1,n)./ctrl.k);
    Kappa = diag(ones(1,n)./ctrl.kappa);
    A = [...
        O I O O O;                      % deltadot = omega
        [O -I O flag.Sec*I O]*flag.PriF;    % omegadot = -omega + xi
        [O O -I O flag.Sec*I]*flag.PriV;    % Vdot = -V + zeta
        [O -K O -sys.Link*K O]*flag.Sec; % xidot = -omega/K -G*xi/K
        [O O -ctrl.Beta.*Kappa O O]*flag.Sec*flag.PriV];   
    
    Aref = [O O O O;
        [I O I.*ctrl.etaP O]*flag.PriF;  % omegadot += omegaast + etaP*Past
        [O I O I.*ctrl.etaQ]*flag.PriV;  % Vdot += Vast + etaQ*Qast
        [K O O O]*flag.Sec;              % xidot += omegast/K
        [O ctrl.Beta.*Kappa O O]*flag.Sec*flag.PriV];    
    
    PQ = [Ov
        -ctrl.etaP.*y.P*flag.PriF          % omegadot -= etaP*P
        -ctrl.etaQ.*y.Q*flag.PriV          % Vdot -= etaQ*Q
        Ov
        -ctrl.H./ctrl.kappa.*QDelta*Iv*flag.Sec*flag.PriV;
        ];
        
    Dist = [Ov; ctrl.etaP.*dist; Ov; Ov; Ov];
end

dydt = A*state + ... % linear dynamics
    Aref*ref + ...     % ref (ast)
    PQ + ...        % nonlinear dynamics (powers)
    Dist;           % Disturbance (power)

if t - tlast > 50
    tlast = t;
    fprintf('%d/%d',t,sys.simtime);
end
end
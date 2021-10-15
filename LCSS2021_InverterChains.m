%% String Stability in Microgrids using Frequency Controlled Inverter Chains
% 
%%% Paper: 
% G. F. Silva, A. Donaire, M. M. Seron, A. McFadyen and J. Ford, 
% "String Stability in Microgrids Using Frequency Controlled Inverter Chains," 
% in IEEE Control Systems Letters, vol. 6, pp. 1484-1489, 2022, doi: 10.1109/LCSYS.2021.3114143.
%
%%% Code:
% 

home
%% Parameters
% Microgrid's electrical parameters
clearvars -except temp sol
rng(1);
load temp

% Toggles
SStest = 1; % if 1, creates length(N) chains with N inverters
savefigs = 0;
voltctrl = 0;
printfigs = 1;
calcnorms = 1;

controllers = ["C1" "C2"];
% C1 - DSS Controller (Proposed)
% C2 - SimpsonPorco2015 (V-Q 'smart' tuning)

ctrl_opt = 23; % Controller structure 2.3 (used in the optimisation function)

if SStest == 1
    N = [4:2:20]; % number of nodes in each inverter chain
else
    N = 8;
end

%%%%%%%%%%%%%%%%%%%%%%% Microgrid topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverter parameters % (for a 4 inverters chain)
Vdc = 650;
C = 25e-6;      % Filter Capacitance  [F]
Lf = 1.8e-3;    % Filter Inductance   [H]
Lo = 1.8e-3;    % Output Impedance    [H]
conf.Prat =[1.4 1.4 1.4 1.4]'; % Rated (max) Active Power   [kW]
conf.Qrat = .8;%[.8 .8 .8 .8]';   % Rated (max) Reactive Power [kVAr]

% Line impedances %
conf.R = [.8 .4 .7]';   % Resistance [Ohm]
conf.L = [3.6 1.8 1.9]'*1e-3; % Reactance  [H]

% Load impedances %
conf.Rl = [1 1 1 1]'*0;       % Resistance [Ohm]
conf.Ll = [1 1 1 1]'*0e-3;     % Reactance  [H] ????

conf.loadconfig = [1 1 1 1]'; 
conf.invRating = ones(4,1);
conf.invRating(conf.loadconfig==0) = 0.5;
conf.invRating(conf.loadconfig==1) = 1;

conf.loaddist = (1 - .0*rand(4,1));

% Load Active [kW] and reactive [kVAr] Power   
conf.Past = conf.invRating.*conf.Prat.*(.9-.0*rand(4,1))*1e3;  
conf.Qast = conf.invRating.*conf.Qrat.*(.9-.0*rand(4,1))*1e3;

% Disturbance
conf.dist = [1.1 1.1 1.1 1.1]';

% Operation set points
conf.omegaast = 50*2*pi; % 50Hz
conf.Vast = 325.3;       % 230V RMS 
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocating memory for variables
invChains(length(N)) = struct(...
    'y',0,...
    'P',0,...
    'Q',0,...
    'Pi',0,...
    'Qi',0,...
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
Yij = 1./Zij;               % Line Admittance
gij = real(Yij);            % Line Conductance
conf.bij = imag(Yij);       % Line Susceptance
conf.Bij = diag(conf.bij,1) + diag(conf.bij,-1);
conf.Gij = diag(gij,1) + diag(gij,-1);
conf.Zl = diag(conf.Rl + conf.Ll*1i);  % Load impedance
conf.Yl = 1./conf.Zl;                  % Load admittance
conf.Yl(conf.Yl==Inf)=0;
conf.Bii = diag(sum(conf.Bij,1)) + imag(conf.Yl);

%%
if exist('sol','var')
    optctrl = sol;
elseif ~exist('ctrl','var')
    try 
        disp('Running optimisation...')
        Optimisation_Vcte;
    catch expt
        error(expt.message);
    end
    if constraints == 0
        clear ctrl sol
        error('constraints not satisfied during optimisation')
    end
    optctrl = sol;
end
    
%% loop
fprintf('Simulating %d strings...\n',length(N)*length(controllers))
for n = N
    sys.n = n;
    fprintf('\t %d/%d (%d) ',ii,length(N)*length(controllers),n)
    
    % Line parameters
    zselect = datasample(1:3,1);
    sys.R = [sys.R; repmat(conf.R(zselect),sys.n-length(sys.R)-1,1)];
    sys.L = [sys.L; repmat(conf.L(zselect),sys.n-length(sys.L)-1,1)];
    % Loads (not currently used, using Pl and Ql instead)
    sys.Rl = [sys.Rl; datasample(conf.Rl,sys.n-length(sys.Rl))];
    sys.Ll = [sys.Ll; datasample(conf.Ll,sys.n-length(sys.Ll))];
    
    % Load configuration
    sys.loadconfig = [sys.loadconfig; ...
  repmat(datasample(conf.loadconfig,1),sys.n-length(sys.loadconfig),1)];
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
    
    % Disturbance    
    sys.dist = [sys.dist; datasample(conf.dist,sys.n-length(sys.dist))];
      
    Zij = sys.R + sys.L*1i;   % Line Impedance
    Yij = 1./Zij;             % Line Admittance
    gij = real(Yij);          % Line Conductance
    bij = imag(Yij);          % Line Susceptance %TODO
    sys.Bij = diag(bij,1) + diag(bij,-1);
    sys.Gij = diag(gij,1) + diag(gij,-1);
    
    sys.Zl = diag(sys.Rl + sys.Ll*1i);  % Load impedance
    sys.Yl = 1./sys.Zl;                 % Load admittance
    sys.Yl(sys.Yl==Inf)=0;
    
    % Shunt admittance
    sys.Gii = diag(sum(sys.Gij,1)) + real(sys.Yl);
    sys.Bii = diag(sum(sys.Bij,1)) + imag(sys.Yl);
    %% %%%%%%%%%%%%%%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    sys.simtime = 3000; % [s]
    % flags
    flag.PriF = 1;
    flag.PriV = voltctrl;
    flag.SecOn = 1;
    flag.SecOnT = 0;
    flag.Dist = 1;
    flag.ic = 0;
    flag.DistOnT = 0;
    flag.DistOffT = 1500;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%%%%%% Controllers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for control = controllers
        fprintf("%s ",control)
        
        % line        
        sys.Comm = diag([1,ones(n-2,1)'],1) + diag([ones(n-2,1)',1],-1);
        
        % ring
        sys.Comm = (diag(ones(1,sys.n-1),1) + ...
            diag(ones(1,sys.n-1),-1) + diag(1,sys.n-1) + diag(1,-sys.n+1));
        
        if control == "C1"            
            ctrl = optctrl;
            if ~voltctrl % Assuming constant voltage
                ctrl.ctrlgain(1) = optctrl.tau;     % tau
                ctrl.ctrlgain(2) = optctrl.etai;    % etaP
                ctrl.ctrlgain(3) = 1;               % (future use)
                ctrl.ctrlgain(4) = optctrl.ki;      % kxi
                ctrl.ctrlgain(5) = 1;               % (future use)
                
                % if additional gains are being used
                try ctrl.ctrlxtra(1) = optctrl.Oo;
                catch 
                    ctrl.ctrlxtra(1) = 1;
                end
                
                try ctrl.ctrlxtra(2) = optctrl.Oxi;
                catch
                    ctrl.ctrlxtra(2) = 1;
                end
                
                try ctrl.ctrlxtra(3) = 1;
                catch
                    ctrl.ctrlxtra(3) = 1;
                end
                
                try ctrl.ctrlxtra(4) = optctrl.Xid;
                catch
                    ctrl.ctrlxtra(4) = 0;
                end
                
                try ctrl.ctrlxtra(5) = optctrl.Xio;
                catch
                    ctrl.ctrlxtra(5) = 1;
                end
                
                ctrl.ctrlxtra(6) = 1;
                
                gim1 = optctrl.comm(1);
%                 gip1 = optctrl.comm(2);
                
                try
                lim1 = optctrl.comm(3);
%                 lip1 = optctrl.comm(4);
                catch
                end
                
                him1 = 0;                
                hip1 = 0;
            else
                gim1 = ctrl.comm(1); 
                gip1 = ctrl.comm(2);                
                him1 = ctrl.comm(3);
                hip1 = ctrl.comm(4);
            end
            
            % Communication
            % - for \xi
            gip1 = gim1;
            lip1 = lim1;
            
            gii = (gip1 + gim1);
            ctrl.G = (diag([gip1 ; gii*ones(n-2,1); gim1])*0 ...
                - gip1*diag([1,ones(n-2,1)'],1)...
                - gim1*diag([ones(n-2,1)',1],-1)...
                - gim1.*diag(1,n-1)...   % those lines implement ring
                - gip1.*diag(1,-n+1));   % connection if first line is zero
            ctrl.G = ctrl.G - diag(sum(ctrl.G,2));
            
            % - for \delta
            try
            lii = lip1 + lim1;
            ctrl.L = (diag([lip1 ; lii*ones(n-2,1); lim1])*0 ...
                - lip1*diag([1,ones(n-2,1)'],1)...
                - lim1*diag([ones(n-2,1)',1],-1)...
                - lim1.*diag(1,n-1)...   % those lines implement ring
                - lip1.*diag(1,-n+1));   % connection if first line is zero
            ctrl.L = ctrl.L - diag(sum(ctrl.L,2));            
               
            catch
            end
            % - for \zeta (if volt control)
            hii = hip1 + him1;
            ctrl.H = diag([hip1 ; hii*ones(n-2,1); him1])...
                - hip1*diag([1,ones(n-2,1)'],1)...
                - him1*diag([ones(n-2,1)',1],-1);
            
        elseif control == "C2"
            %% Communication
            % line
            sys.Link = (diag([1 ; 2*ones(n-2,1); 1]) ...
                - 1*diag([1,ones(n-2,1)'],1)...
                - 1*diag([ones(n-2,1)',1],-1));
            
            % ring
            sys.Link = (diag(ones(1,sys.n-1),1) + ...
                diag(ones(1,sys.n-1),-1) + diag(1,sys.n-1) + diag(1,-sys.n+1));
            sys.Link = -sys.Link + diag(sum(sys.Link,2));
            %% Primary Freq Controller %%%
            ctrl.etaP = 3e-3;%1./sys.Past/2;
            ctrl.etaQ = 1./sys.Qast./800;

            %% Secondary Freq Controller %%%
            ctrl.k = 1.*ones(1,sys.n); % Int Frequency gain
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
        delta0 = [-.1; (-.1+.1*rand(n-1,1))]*flag.ic; zeros(n,1);
        omega0 = sys.omegaast*ones(n,1) + [.5; 0.1*rand(n-1,1)]*flag.ic;
        V0 = (sys.Vast)*ones(n,1);% - 1 + 2*rand(n,1);
        xi0 = zeros(n,1);
        zeta0 = 0*ones(n,1);
        state0 = [delta0  % volt angle (delta)
            omega0        % freq (omega)
            V0            % volt (V)
            xi0           % int freq (xi)
            zeta0         % int volt (zeta)
            ];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%       
        ctrl.opt = ctrl_opt;
        %% %% ODE %%%%%%%%%%%%%%%%%%        
        [time,y] = ode15s(@microgrid,[0 sys.simtime], state0, opts, sys,...
            ctrl, flag);
        
        % u_final = zeros(length(time),sys.n*2);
        P_final = zeros(length(time),sys.n);
        Q_final = zeros(length(time),sys.n);
        Pi_final = zeros(length(time),sys.n);
        Qi_final = zeros(length(time),sys.n);
        for tt = time'
            idx_t = find(ismember(time,tt));
            [dx, yy] = microgrid(tt,y(idx_t,:)', sys, ctrl, flag);
            %     u_final(ttt,:) = yy.u;
            P_final(idx_t,:) = yy.P;
            Q_final(idx_t,:) = yy.Q;
            Pi_final(idx_t,:) = yy.Pi;
            Qi_final(idx_t,:) = yy.Qi;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % - Populating result structure - %
        invChain.y = y;
        invChain.P = P_final;
        invChain.Q = Q_final;
        invChain.Pi = Pi_final;
        invChain.Qi = Qi_final;
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
    fprintf("\n")
end

%% Plotting
plots.powerinj = 1;
plots.delta = 0;
plots.omega = 1;
plots.xi = 0;
plots.P = 1;
plots.Q = 0;
plots.Pfactor = 0;

Cplot = controllers(2);
Nplot = 4;

plotting

%% % --- Calculating Norms --- %%%
normscalc

normsplot

%% Printing/Saving figures
% clc
if savefigs == 1    
    figs = findobj(0, 'type', 'figure');
    figs = figs(end:-1:1);
    datetime_char = string(datetime('now','format','yyyyMMMdd_hhmm'));
    for k=1:length(figs)
        % print each figure in figs to a separate .eps file
        print(figs(k), '-depsc', ...
            sprintf("figs/N%d%s_%s_%s.eps", Nplot, Cplot, datetime_char, figs(k).Name))        
    end
    savefig(figs,sprintf("figs/N%d%s_%s.fig", Nplot, Cplot, datetime_char))
    disp('Figures saved.')
end

%%
disp('*** End of Script ***') 
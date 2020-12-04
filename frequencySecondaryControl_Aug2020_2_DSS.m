%% Simulation script for inverter strings
% All values are in Volt, Watt, ohm, second, etc

%% Settings
clear invChains
rng(1);

save_results = 0;
desc = 'v1';

plot_results = 1;
save_plots = 0;
wrap = 0;

% Number of inverters
N = [4];
% N=10;
    
vars = 'absolute';
vars = 'incremental';
controllers = ["DSS"];
% controllers = [controllers "DSS"];
control_plot = "Org"; 
control_plot = "DSS";
simtime = containers.Map({'DSS','Org'},{1, 10});
PlDist = containers.Map({'DSS','Org'},{0, 0}); % should be the same!

% Set to 1 if random values of rated powers and loads for different tests
% are sought, r = 0 uses Simpson-Porco (2015)
r = 0;
if r == 0
    load inverters
end

%% Allocating memory for variables
invChains(length(N)) = struct(...    
            'y',0,...
            'u',0,...
            'P',0,...
            'Q',0,...
            'time',0,...
            'sys',0,...         
            'ctrl',0,...             'z_normt',z_l2normt,...            'z_infnorm',z_infnorm,...
            'control',0,...
            'colvec', 0);

%         'y',0,...
%         'u',0,...
%         'time',0,...
%         'sys',0,...
%         'ctrl',0,...
%         'z_normt',0,...
%         'z_infnorm',0,...
%         'control',0,...
%         'colvec', 0);

%% Parameters 
% Per unit bases
Pbase = 1; % in kilovolt-ampere [kVA]
Qbase = 1; % in kilovolt-ampere [kVAr] (?)
Vbase = 325.3; % in volt [V]
Zbase = Vbase*Vbase/(Pbase*1000);
Ybase = 1/Zbase;

% The desired reference frequency
freq = 50; %[Hz]

% The rated power of the inverters (Pos?)
PastPos = [1.4, .7, .7, 1.4]; % from Simpson-Porco (2015)
QastPos = [.8, .4, .4, .8]; % from Simpson-Porco (2015)

% The active power of loads (this corresponds to the error P_i^*-P_{Li}
PlPos = 0*[1, 0.5, 0.5, 1]; % Simpson-Porco (2015) does not report Pl...
% The reactive power of loads (corresponds to Q_i^* - Q_L)
QlPos = 0*[.1, 0.05, 0.05, .1]; % from Majumder (2009)

% The inductive part of the lines is a fraction of 
bij = -1/0.6; % from Majumder (2009)
% per unit
bij = bij/Ybase; %-1/3.488;  

% A particular set of values for the line impedances
%bb=[-0.4073   -8.7218  -59.8396 -130.5942  -36.5307]';



%% Simulating
i = 1;
for n = N    
    %% Getting the controller
    if exist('sol','var')
        ctrl = sol;
    elseif ~exist('ctrl','var')
        try % Optimisation 
            Optimisation;
        catch exp
            error(exp.message);
        end 
        if constraints == 0
            clear ctrl sol
            error('constraints not satisfied during optimisation')
        end
        ctrl = sol
    end
    
    if r==1
        Past = datasample(PastPos,n)';
        Qast = datasample(QastPos,n)';
        Pl = datasample(PlPos,n)';
        Ql = datasample(QlPos,n)';
    else
        Past = PastPos';
        Qast = QastPos';
        Pl = PlPos';
        Ql = QlPos';
    end
    
    
    % per unit
    Past = Past/Pbase;
    Qast = Qast/Qbase;
%     Past = ones(n,1);

    % Choosing appropriate values for the controllers (from Simpson-Porco (2015))
    if r == 1
        tau = 10*ones(1,n);
        etap = 1./Past/100;
    else
        tau = 10*ones(1,n);
        etap = [2.5e-3 5e-3 5e-3 2.5e-3];
        etaq = [1.5e-3 3e-3 3e-3 1.5e-3];
    end
    
    k = 1.7 * ones(1,n);
    kappai = ones(1,n);

    % per unit
    Pl = Pl/Pbase*1;
    Ql = Ql/Pbase*1;
    
%     Pl(end) = 0;
    
    % String connection of inverters assuming no self connection
    if r==1
        bb=bij*rand(1,n-1);
    end

    B = diag(bb,1) + diag(bb,-1); 
    
    % The desired reference frequency
    wref = freq*2*pi* ones(n,1);

    % The voltages are assumed to remain constant at the rated voltage of the
    % distribution line
    % absolute voltage
    voltage = 325.3;
    % per unit
    voltage = voltage/Vbase; %find_voltages(B,Past,Pl,1,n);
    VV = voltage*voltage*ones(n,n);     % alpha

    %% Choose the initial conditions
    if r == 1
        dw = [-1; (rand(n-1,1)-ones(n-1,1))];
    end
    
    xi0 = 0*rand(n,1); % secondary control variable (freq)
    zeta0 = 0*rand(n,1); % secondary control variable (volt)
    if strcmp(vars,'incremental')
        theta0 = 0*ones(n-1,1);
        omega0 = 0*rand(n,1);
        v0 = 1*rand(n,1);
        initial_state = [theta0; omega0; v0; xi0; zeta0];
    else
        theta0 = ones(n-1,1);
        omega0 = wref + wd;
        v0 = voltage*rand(n,1);
        initial_state = [delta0; omega0; v0; xi0; zeta0];
    end
    
    % Defining system struct
    sys = struct(...
        'n', n,...
        'V', voltage,...
        'wref', wref,...
        'B', B,...
        'Pl', Pl,...
        'Ql', Ql);
    
    opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
    for control = controllers
        sys.simtime = simtime(control);        
        sys.PlDist = PlDist(control);
        
        % Communication control in the same string
        if control == "DSS"
            gim1 = ctrl.comm(1);
            gip1 = ctrl.comm(2);
            him1 = ctrl.comm(3);
            hip1 = ctrl.comm(4);
            
            gii = gip1 + gim1;
            hii = hip1 + him1;
            G = diag([gip1 ; gii*ones(n-2,1); gim1])...
                - gip1*diag([1,ones(n-2,1)'],1)...
                - gim1*diag([ones(n-2,1)',1],-1);
            H = diag([hip1 ; hii*ones(n-2,1); him1])...
                - hip1*diag([1,ones(n-2,1)'],1)...
                - him1*diag([ones(n-2,1)',1],-1);
            
        else
            G = diag([.5 ; ones(n-2,1);.5])...
                - 0.5*diag([1,ones(n-2,1)'],1)...
                - 0.5*diag([ones(n-2,1)',1],-1);
            H = diag([.5 ; ones(n-2,1);.5])...
                - 0.5*diag([1,ones(n-2,1)'],1)...
                - 0.5*diag([ones(n-2,1)',1],-1);
        end
        sys.G = G;
        sys.H = H;
        
        if strcmp(vars,'incremental')
            if control == "DSS"
                [time,y] = ode45(@inverters_theta_dss,[0 sys.simtime], ...
                    initial_state, opts, sys, ctrl);
                % outputing the controller
                u_final = zeros(length(time),n*2);
                P_final = zeros(length(time),n);
                Q_final = zeros(length(time),n);
                for tt = time'
                    ttt = find(ismember(time,tt));
                    [dx, u, P, Q] = inverters_theta_dss(tt,y(ttt,:)',sys,ctrl);
                    u_final(ttt,:) = u;
                    P_final(ttt,:) = P;
                    Q_final(ttt,:) = Q;
                end
            else
%                 [time,y] = ode45(@inverters_theta,[0 sys.simtime],...
%                     initial_state,opts,sys,tau,etap,k);
%                 % outputing the controller
%                 u_final = zeros(length(time),n);
%                 for tt = time'
%                     ttt = find(ismember(time,tt));
%                     [dx, u] = inverters_theta(tt,y(ttt,:)',sys,tau,etap,k);
%                     u_final(ttt,:) = u;
%                 end
            end
        else 
%             if control == "DSS"
%                 [time,y] = ode45(@inverters_dss,[0 sys.simtime], ...
%                     initial_state, opts,sys,ctrl);
%                 % outputing the controller
%                 u_final = zeros(length(time),n);
%                 for tt = time'
%                     ttt = find(ismember(time,tt));
%                     [dx, u] = inverters_dss(tt,y(ttt,:)',sys,ctrl);
%                     u_final(ttt,:) = u;
%                 end
%             else
%                 [time,y] = ode45(@inverters,[0 sys.simtime],...
%                     initial_state,opts,sys,tau,etap,k);
%                 % outputing the controller
%                 u_final = zeros(length(time),n);
%                 for tt = time'
%                     ttt = find(ismember(time,tt));
%                     [dx, u] = inverters(tt,y(ttt,:)',sys,tau,etap,k);
%                     u_final(ttt,:) = u;
%                 end
%             end
        end

        Frequency_error=y(end,[n+1:2*n])-wref';
        Secondary_steady_state_value=y(end,end-n+1:end);

        V = [eye(n-1) -ones(n-1,1)];

        % Calculating norms
%         znorm = zeros(length(y),n);
%         if strcmp(vars,'incremental')
%             z = reshape([zeros(length(time),1) y], length(time),n,(size(y,2)+1)/n);        
%             for idx_t = time'
%                 tt = find(ismember(time,idx_t));
%                 for idx_a = 1:n
%                     ztemp = z(tt,idx_a,:);
%                     ztemp = reshape(ztemp,1,3);
% 
%                     % angle error
%                     if wrap
%                         ztemp(1) = wrapToPi(ztemp(1));
%                     end
%                     znorm(tt,idx_a) = norm(ztemp);
%                 end
%             end
%             z_l2normt = max(znorm,[],2);
%             [z_infnorm, ~] = max(max(znorm));
%         else
%             z = reshape(y, length(time),n,(size(y,2))/n);
% 
%             for idx_t = time'
%                 tt = find(ismember(time,idx_t));
%                 for idx_a = 1:n
%                     ztemp = z(tt,idx_a,:);
%                     ztemp = reshape(ztemp,1,3);
% 
%                     % angle error
% %                     if wrap
%                         ztemp(1) = wrapToPi(ztemp(1));
% %                     end
% 
%                     % freq error
%                     ztemp(2) = ztemp(2) - wref(1);
%                     znorm(tt,idx_a) = norm(ztemp);
%                 end
%             end
%             z_l2normt = max(znorm,[],2);
%             [z_infnorm, ~] = max(max(znorm));            
%         end
        

        % Populating result structure
        invChain = struct(...
            'y',y,...
            'u',u_final,...
            'P',P_final,...
            'Q',Q_final,...
            'time',time,...
            'sys',sys,...         
            'ctrl',ctrl,...             'z_normt',z_l2normt,...            'z_infnorm',z_infnorm,...
            'control',control,...
            'colvec', colormap(hsv(n+1)));

        % accumulating final results
        invChains(i) = invChain;
        i = i + 1;
    end
end

%% Saving results
if save_results == 1
    for control = controllers
        datetime_char = char(datetime('now','format','yyMMdd_HHmm'));
        save(['results/' datetime_char '_' control{:} '_' desc],'invChains','N')
        save(['results/' datetime_char '_' control{:} '_' desc '_ctrl'],'ctrl')
    end
end

%% Plotting
% load results\200918_1523_DSS_result.mat
% load results\200918_1523_Org_result.mat
if plot_results == 1
    ChainLength = N(1);
    ChainToPlot = find(ismember(N,ChainLength));    
%     control_plot = "DSS";
    if length(controllers)>1&&control_plot == "DSS"
        ChainToPlot = ChainToPlot + 1;
    end
    plotResults(invChains,ChainToPlot,N,control_plot,vars,wrap);
    if save_plots
        for i = 1:6
            figure(i); 
            print([control_plot{:} '_' desc '_' num2str(i)],'-depsc');
            savefig([control_plot{:} '_' desc '_' num2str(i)])
        end
    end
end

%% ------------------------------------------------------------------------
function [dydt, u, P, Q] = inverters(t, state, sys, tau, eta, k)
% Implements the inverter system in the absolute variables
n = sys.n;
wref = sys.wref;
simtime = sys.simtime;

B = sys.B;
Lambda = sys.Lambda;

Pl = sys.Pl;
PlDist = sys.PlDist;

I = eye(n);
O = zeros(n,n);
Ii= ones(n,1);

Ti = diag(ones(1,n)./tau);
N = diag(eta);
Ki = 1*diag(ones(1,n)./k);

Pldist = 0;
if t>simtime/2 % constant disturbance in Pl after half simulation time
    Pldist = PlDist;
end
Pll = Pl(1:end-1) + Pldist;
% Pll(end+1) = 0;
Pll(end+1) = Pl(end) + Pldist;

A = [O, I, O;
    O, -Ti, Ti;
    O, -Ki, -Ki*Lambda];   
Aref =[O, O;
    Ti, Ti*N;
    Ki, O];    

%state(1:n) = wrapToPi(state(1:n));

Delta = repmat(state(1:n),1,n) - repmat(state(1:n)',n,1); 
% Delta is defined just for F, as F = bij*V^2*sin(delta_i-delta_j)
F=-(B.*sin(Delta))*Ii; % 

dydt = A * state + ...  % free term 
    Aref * [wref;Pll] ...reference (accounts for -w* in (21)) + disturbance
       + [zeros(n,1);Ti*N*F;zeros(n,1)]; % non-linear term
%dydt(1:n-1)= wrapToPi(dydt(1:n-1));

% outputing controller
Au = [zeros(n), -Ti, Ti];
u = Au*state + Ti*N*F;
end

function [dydt, u, P, Q] = inverters_dss(t, state, sys, ctrl)
% Implements the inverter system in the absolute variables with extra
% controller gains to satisfy Monteil et al's DSS conditions.

n = sys.n;
wref = sys.wref;
simtime = sys.simtime;

B = sys.B;
Lambda = sys.Lambda;

Pl = sys.Pl;
PlDist = sys.PlDist;

I = eye(n);
O = zeros(n,n);
Ii= ones(n,1);

Ti = diag(ones(1,n)./ctrl.tau);
N = diag(ctrl.etai);
Ki = 1*diag(ones(1,n)./ctrl.ki);

Pldist = 0;
if t>simtime/2 % constant disturbance in Pl after half simulation time
    Pldist = PlDist;
end
Pll = Pl(1:end-1) + Pldist;
% Pll(end+1) = 0;
Pll(end+1) = Pl(end) + Pldist;

A = [O, I, O;
    O, -Ti*ctrl.Oi, Ti*ctrl.Zi;
    O, -Ki*ctrl.Li, -Ki*Lambda];   

Aref =[O, O;
    Ti*ctrl.Oi, Ti*N;
    Ki*ctrl.Li, O];    

%state(1:n) = wrapToPi(state(1:n));

Delta = repmat(state(1:n),1,n) - repmat(state(1:n)',n,1); 
% Delta is defined just for F, as F = bij*V^2*sin(delta_i-delta_j)
F=-(B.*sin(Delta))*Ii; 

dydt = A * state + ...  % free term 
    Aref * [wref;Pll/ctrl.Oi] ...reference (accounts for -w* in (21)) + disturbance
       + [zeros(n,1);Ti*N*F;zeros(n,1)]; % non-linear term
%dydt(1:n-1)= wrapToPi(dydt(1:n-1));

% outputing controller
Au = [zeros(n,n), -Ti*ctrl.Oi, Ti*ctrl.Zi];
u = Au*state + Ti*N*F;
end

function [dydt, u, P, Q] = inverters_theta(t, state, sys, tau, eta, k)
% Implements the inverter system in the incremental variables as per
% Schiffer et al's transformations

n = sys.n;
wref = sys.wref;
simtime = sys.simtime;

B = sys.B;
Lambda = sys.Lambda;

Pl = sys.Pl;
PlDist = sys.PlDist;

I = eye(n);
O = zeros(n,n);
Ii= ones(n,1);
V=[eye(n-1) -ones(n-1,1)];

Ti = diag(ones(1,n)./tau);
N = diag(eta);
Ki = 1*diag(ones(1,n)./k);

Pldist = 0;
if t>150 % constant disturbance in Pl after half simulation time
    Pldist = PlDist;
end
Pll = Pl(1:end-1) + Pldist;
% Pll(end+1) = 0;
Pll(end+1) = Pl(end) + Pldist;

A = [zeros(n-1,n-1), V, 0*V;
    zeros(n,n-1), -Ti*1, Ti;
    zeros(n,n-1), -Ki, -Ki*Lambda];

Delta = repmat(state(1:n-1),1,n-1) - repmat(state(1:n-1)',n-1,1);
FF=[(-B(1:n-1,1:n-1).*sin(Delta(1:n-1,:)))*Ii(1:n-1); -B(n,n-1)*sin(state(n-1))];

dydt = A * state...
    + [zeros(n-1,1);Ti*(N*(Pll));zeros(n,1)] ... % disturbance (Gill: ci)
    + [zeros(n-1,1);Ti*N*FF;zeros(n,1)]; % non-linear term

Au = [zeros(n,n-1), -Ti, Ti];

u = Au*state + Ti*N*FF;
end

function [dydt, u, P, Q] = inverters_theta_dss(t, state, sys, ctrl)
% Implements the inverter system in the incremental variables as per
% Schiffer et al's transformations, with extra controller gains to satisfy
% Monteil et al's DSS conditions. 
tau = ctrl.ctrlgain(1);
np = ctrl.ctrlgain(2);
nq = ctrl.ctrlgain(3);
kxi = ctrl.ctrlgain(4);
kzeta = ctrl.ctrlgain(5);
l = ctrl.ctrlxtra(1);
m = ctrl.ctrlxtra(2);
o = ctrl.ctrlxtra(3);
x = ctrl.ctrlxtra(4);
y = ctrl.ctrlxtra(5);

n = sys.n;
simtime = sys.simtime;

B = sys.B;
G = sys.G;
H = sys.H;

Pl = sys.Pl;
Ql = sys.Ql;

PlDist = sys.PlDist;

I = eye(n);
Os = zeros(n,n);
Ii = ones(n,1);
R = [eye(n-1) -ones(n-1,1)];

T = diag(ones(1,n)./tau);
Np = diag(np);
Nq = diag(nq);
Kxi = diag(ones(1,n)./kxi);
Kzeta = diag(ones(1,n)./kzeta);
L = diag(ones(1,n)*l);
M = diag(ones(1,n)*m);
O = diag(ones(1,n)*o);
X = diag(ones(1,n)*x);
Y = diag(ones(1,n)*y);

% constant disturbance in Pl 
Pldist = 0;
if t>150 
    Pldist = PlDist;
end
Pll = Pl(1:end-1) + Pldist;
% Pll(end+1) = 0;
Pll(end+1) = Pl(end) + Pldist;

A = [
    zeros(n-1,n-1) R      0*R      0*R   0*R;
    zeros(n,n-1)  -T      Os       T*Y   Os;
    zeros(n,n-1)   Os    -T*O      Os    T*M;
    zeros(n,n-1)  -Kxi*X  Os      -Kxi*G Os;
    zeros(n,n-1)   Os    -Kzeta*L  Os   -Kzeta*H;    
    ];

% Delta below is theta_i - theta_i+1
Delta = repmat(state(1:n-1),1,n-1) - repmat(state(1:n-1)',n-1,1);
VV = repmat(state(2*n:3*n-1),1,n).*repmat(state(2*n:3*n-1)',n,1);

FFp=[(-B(1:n-1,1:n-1).*sin(Delta(1:n-1,:)))*Ii(1:n-1).*diag(VV,1);
    -B(n,n-1)*sin(state(n-1))].*VV(n,n-1);

FFq=[diag(VV(1:n-1,1:n-1)).^2 - (B(1:n-1,1:n-1).*cos(Delta(1:n-1,:)))*Ii(1:n-1).*diag(VV,1);
    VV(n,n) - B(n,n-1)*cos(-state(n-1)).*VV(n,n-1)];
% FFq=[- (B(1:n-1,1:n-1).*cos(Delta(1:n-1,:)))*Ii(1:n-1);
%     - B(n,n-1)*cos(state(n-1))];

dydt = A * state...
    + [zeros(n-1,1); T*(Np*(Pll)); T*(Nq*Ql); zeros(n,1); zeros(n,1)] ... % disturbance
    + [zeros(n-1,1); T*Np*FFp; T*Nq*FFq; zeros(n,1); zeros(n,1)];


Au = [zeros(n,n-1)  -T      Os       T*Y   Os;
    zeros(n,n-1)   Os    -T*O      Os    T*M];

u = Au*state + [T*Np*FFp; T*Nq*FFq];

P = Pll + FFp;
Q = Ql + FFq;
end
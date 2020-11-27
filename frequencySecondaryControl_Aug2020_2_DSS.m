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
N = [5];
% N=10;
    
vars = 'absolute';
vars = 'incremental';
controllers = ["Org"];
controllers = [controllers "DSS"];
control_plot = "Org"; 
control_plot = "DSS";
simtime = containers.Map({'DSS','Org'},{300, 300});
PlDist = containers.Map({'DSS','Org'},{500000, 50}); % should be the same!

% Set to 1 if random values of rated powers and loads for different tests
% are sought, r = 0 uses a saved example with n = 6 
r = 1;
if r == 0
    load inverters
end

%% Allocating memory for variables
invChains(length(N)) = struct(...
        'y',0,...
        'u',0,...
        'time',0,...
        'sys',0,...
        'ctrl',0,...
        'z_normt',0,...
        'z_infnorm',0,...
        'control',0,...
        'colvec', 0);

%% Parameters 
% Per unit bases
Pbase = 1; % in kilovolt-ampere [kVA]
Vbase = 415; % in volt [V]
Zbase = Vbase*Vbase/(Pbase*1000);
Ybase = 1/Zbase;

% The desired reference frequency
freq = 50; %[Hz]

% The rated power of the inverters (Pos?)
% PastPos = [1, 2, 1.5]; % from Majumder (2009) (original)
PastPos = [1, 2, 1.5, 1.5]; % from Majumder (2009)
% The active power of loads (this corresponds to the error P_i^*-P_{Li}
% PlPos = [1.8, 0.8, 0,3]; % original
PlPos = [1.8, 0.8, 0.8, 0.8, 1.8];

% The reactive power of loads (corresponds to Q_i^* - Q_L)
QlPos = [1.6, 0.6, 0.6, 0.6, 1.6]; % from Majumder (2009)

% The inductive part of the lines is a fraction of 
bij = -1/0.6; % from Majumder (2009)
% per unit
bij = bij/Ybase; %-1/3.488;  

% A particular set of values for the line impedances
%bb=[-0.4073   -8.7218  -59.8396 -130.5942  -36.5307]';

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
% Removing the additional term in the secondary variable dynamics
ctrl.Yi = 0; 
ctrl.Oi = 1;

%% Simulating
i = 1;
for n = N    
    if r == 1
        Past = datasample(PastPos,n)';
    end
    % per unit
    Past = Past/Pbase;
%     Past = ones(n,1);

    % Choosing appropriate values for the controllers
    tau = 10*ones(1,n);
    eta = 1./Past/100;
    k = 1.7 * ones(1,n);

    % A particular load condition
    %Pl=[1.8000    1.8000    0.8000    1.8000    3.0000  0]';

    if r == 1
        Pl = datasample(PlPos,n)';
        Ql = datasample(QlPos,n)';
    end

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
    voltage = 415;
    % per unit
    voltage = voltage/Vbase; %find_voltages(B,Past,Pl,1,n);
    VV = voltage*voltage*ones(n,n);     % alpha

    %%% Choose the initial conditions
    xi0 = 0*rand(n,1); % secondary control variable

    if r == 1
        dw = [-1; (rand(n-1,1)-ones(n-1,1))];
        save inverters Past bb dw Pl
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
            ctrl.gii = ctrl.gip1 + ctrl.gim1;
            Lambda = diag([ctrl.gip1 ; ctrl.gii*ones(n-2,1); ctrl.gim1])...
                - ctrl.gip1*diag([1,ones(n-2,1)'],1)...
                - ctrl.gim1*diag([ones(n-2,1)',1],-1);
        else
            Lambda = diag([.5 ; ones(n-2,1);.5])...
                - 0.5*diag([1,ones(n-2,1)'],1)...
                - 0.5*diag([ones(n-2,1)',1],-1);
        end
        sys.Lambda = Lambda;
        
        if strcmp(vars,'incremental')
            initial_state = [ones(n-1,1);1*rand(n,1);xi0];  %% for incremental
            if control == "DSS"
                [time,y] = ode23s(@inverters_theta_dss,[0 sys.simtime], ...
                    initial_state, opts, sys, ctrl);
                % outputing the controller
                u_final = zeros(length(time),n);
                for tt = time'
                    ttt = find(ismember(time,tt));
                    [dx, u] = inverters_theta_dss(tt,y(ttt,:)',sys,ctrl);
                    u_final(ttt,:) = u;
                end
            else
                [time,y] = ode23s(@inverters_theta,[0 sys.simtime],...
                    initial_state,opts,sys,tau,eta,k);
                % outputing the controller
                u_final = zeros(length(time),n);
                for tt = time'
                    ttt = find(ismember(time,tt));
                    [dx, u] = inverters_theta(tt,y(ttt,:)',sys,tau,eta,k);
                    u_final(ttt,:) = u;
                end
            end
        else 
            initial_state = [1*ones(n,1);wref+dw;xi0];           %% for absolute
            if control == "DSS"
                [time,y] = ode45(@inverters_dss,[0 sys.simtime], ...
                    initial_state, opts,sys,ctrl);
                % outputing the controller
                u_final = zeros(length(time),n);
                for tt = time'
                    ttt = find(ismember(time,tt));
                    [dx, u] = inverters_dss(tt,y(ttt,:)',sys,ctrl);
                    u_final(ttt,:) = u;
                end
            else
                [time,y] = ode45(@inverters,[0 sys.simtime],...
                    initial_state,opts,sys,tau,eta,k);
                % outputing the controller
                u_final = zeros(length(time),n);
                for tt = time'
                    ttt = find(ismember(time,tt));
                    [dx, u] = inverters(tt,y(ttt,:)',sys,tau,eta,k);
                    u_final(ttt,:) = u;
                end
            end
        end

        Frequency_error=y(end,[n+1:2*n])-wref';
        Secondary_steady_state_value=y(end,end-n+1:end);

        V = [eye(n-1) -ones(n-1,1)];

        % Calculating norms
        znorm = zeros(length(y),n);
        if strcmp(vars,'incremental')
            z = reshape([zeros(length(time),1) y], length(time),n,(size(y,2)+1)/n);        
            for idx_t = time'
                tt = find(ismember(time,idx_t));
                for idx_a = 1:n
                    ztemp = z(tt,idx_a,:);
                    ztemp = reshape(ztemp,1,3);

                    % angle error
                    if wrap
                        ztemp(1) = wrapToPi(ztemp(1));
                    end
                    znorm(tt,idx_a) = norm(ztemp);
                end
            end
            z_l2normt = max(znorm,[],2);
            [z_infnorm, ~] = max(max(znorm));
        else
            z = reshape(y, length(time),n,(size(y,2))/n);

            for idx_t = time'
                tt = find(ismember(time,idx_t));
                for idx_a = 1:n
                    ztemp = z(tt,idx_a,:);
                    ztemp = reshape(ztemp,1,3);

                    % angle error
%                     if wrap
                        ztemp(1) = wrapToPi(ztemp(1));
%                     end

                    % freq error
                    ztemp(2) = ztemp(2) - wref(1);
                    znorm(tt,idx_a) = norm(ztemp);
                end
            end
            z_l2normt = max(znorm,[],2);
            [z_infnorm, ~] = max(max(znorm));            
        end
        

        % Populating result structure
        invChain = struct(...
            'y',y,...
            'u',u_final,...
            'time',time,...
            'sys',sys,...         
            'ctrl',ctrl,...
            'z_normt',z_l2normt,...
            'z_infnorm',z_infnorm,...
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

n = sys.n;
wref = sys.wref;
simtime = sys.simtime;

B = sys.B;
Lambda = sys.Lambda;

Pl = sys.Pl;
PlDist = sys.PlDist;

I = eye(n);
O = zeros(n,n);
Ii = ones(n,1);
V = [eye(n-1) -ones(n-1,1)];

Ti = diag(ones(1,n)./ctrl.tau);
N = diag(ctrl.etai);
Ki = 1*diag(ones(1,n)./ctrl.ki);

Pldist = 0;
if t>150 % constant disturbance in Pl 
    Pldist = PlDist;
end
Pll = Pl(1:end-1) + Pldist;
% Pll(end+1) = 0;
Pll(end+1) = Pl(end) + Pldist;

A = [zeros(n-1,n-1), V, 0*V;
    zeros(n,n-1), -Ti*ctrl.Oi, Ti*ctrl.Zi;
    zeros(n,n-1), -Ki*ctrl.Li, -Ki*Lambda];

Delta = repmat(state(1:n-1),1,n-1) - repmat(state(1:n-1)',n-1,1);
FF=[(-B(1:n-1,1:n-1).*sin(Delta(1:n-1,:)))*Ii(1:n-1);
    -B(n,n-1)*sin(state(n-1))];

dydt = A * state...
    + [zeros(n-1,1); Ti*(N*(Pll)); zeros(n,1)] ... % disturbance
    + [zeros(n-1,1); Ti*N*FF; Ki*ctrl.Yi*FF];
       
Au = [zeros(n,n-1), -Ti*ctrl.Oi, Ti*ctrl.Zi];

u = Au*state + Ti*N*FF;

% P = Pll + V^2*sin(Delta);
% Q = Qll + V^2*cos(Delta);
end
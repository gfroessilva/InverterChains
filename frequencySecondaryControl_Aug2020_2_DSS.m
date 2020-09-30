%% Simulation script for inverter strings
% All values are in Volt, Watt, ohm, second, etc

%% Settings
clear invChains
rng(1);

save_results = 1;
plot_results = 1;

% Allocating memory for variables
invChains(length(N)) = struct(...
        'y',0,...
        'u',0,...
        'time',0,...
        'n',0,...         'ctrl','',...
        'z_normt',0,...
        'z_infnorm',0,...
        'control',0,...
        'colvec', 0);
    
% Number of inverters
% N = [50 100 150 200 250 300 350 400 450 500];
N = [50 100 150 200];
vars = 'absolute';
vars = 'incremental';
controllers = ["Org" "DSS"];
simulation_time = 0:1000;

% Set to 1 if random values of rated powers and loads for different tests
% are sought, r = 0 uses a saved example with n = 6 
r = 1;
if r == 0
    load inverters
end

%% Parameters 
% Per unit bases
Pbase = 1; % in kilovoltaampere [kVA]
Vbase = 415; % in volt [V]
Zbase = Vbase*Vbase/(Pbase*1000);
Ybase = 1/Zbase;

% The rated power of the inverters
PastPos = [1,2,1.5]; % from Majumder (2009)
% The active power of loads (this corresponds to the error P_i^*-P_{Li}
PlPos = [1.8,0.8,0,3];

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
        error('constraints not satisfied during optimisation')
    end
    ctrl = sol;
end
% Removing the additional term in the secondary variable dynamics
ctrl.Yi = 0; 

%% Simulating
i = 1;
for n = N    
    if r == 1
        Past = datasample(PastPos,n)';
    end
    % per unit
    Past = Past/Pbase;

    % Choosing appropriate values for the controllers
    tau = 10*ones(1,n);
    eta = 1./Past/100;
    k = 1.7 * ones(1,n);

    % A particular load condition
    %Pl=[1.8000    1.8000    0.8000    1.8000    3.0000  0]';

    if r == 1
        Pl = datasample(PlPos,n)';
    end

    % per unit
    Pl = Pl/Pbase;

    % String connection of inverters assuming no self connection
    if r==1
        bb=bij*rand(1,n-1);
    end

    B = diag(bb,1) + diag(bb,-1); 

    % The desired reference frequency
    wref = 50*2*pi* ones(n,1);

    % The voltages are assumed to remain constant at the rated voltage of the
    % distribution line
    % absolute voltage
    voltage = 415;
    % per unit
    voltage = voltage/Vbase; %find_voltages(B,Past,Pl,1,n);
    VV = voltage*voltage*ones(n,n);     % alpha

    %%% Choose the initial conditions
    z0 = 0*rand(n,1);

    if r == 1
        dw = [-1; (rand(n-1,1)-ones(n-1,1))];
        save inverters Past bb dw Pl
    end
    
    opts = odeset('RelTol',1e-6,'AbsTol',1e-5);
    for control = controllers
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
        if strcmp(vars,'incremental')
            initial_state = [0*ones(n-1,1);rand(n,1);z0];  %% for incremental
            if control == "DSS"
                [time,y] = ode45(@inverters_theta_dss,[0 simulation_time(end)], ...
                    initial_state, opts,n,ctrl,Pl,B,Lambda);
                % outputing the controller
                u_final = zeros(length(time),n);
                for tt = time'
                    ttt = find(ismember(time,tt));
                    [dx, u] = inverters_theta_dss(tt,y(ttt,:)',n,ctrl,Pl,B,Lambda);
                    u_final(ttt,:) = u;
                end
            else
                [time,y] = ode45(@inverters_theta,[0 simulation_time(end)],...
                    initial_state,opts,n,tau,eta,k,Pl,B,Lambda);
                % outputing the controller
                u_final = zeros(length(time),n*3-1);
                for tt = time'
                    ttt = find(ismember(time,tt));
                    [dx, u] = inverters_theta(tt,y(ttt,:)',n,tau,eta,k,Pl,B,Lambda);
                    u_final(ttt,:) = u;
                end
            end
        else
            initial_state = [0*ones(n,1);wref+dw;z0];           %% for absolute
            [time,y] = ode45(@(tt,yy)inverters(yy,n,tau,eta,k,Pl,B,wref,Lambda),...
            simulation_time,initial_state);
        end

        Frequency_error=y(end,[n+1:2*n])-wref';
        Secondary_steady_state_value=y(end,end-n+1:end);


        V = [eye(n-1) -ones(n-1,1)];

        % Calculating norms
        znorm = zeros(length(y),n);
        if strcmp(vars,'incremental')
            z = reshape([zeros(length(time),1) y], length(time),n,(size(y,2)+1)/n);        
        else
            z = reshape(y, length(time),n,(size(y,2))/n);
        end
        for idx_t = time'
            tt = find(ismember(time,idx_t));
            for idx_a = 1:n
                ztemp = z(tt,idx_a,:);
                ztemp = reshape(ztemp,1,3);

                % angle error
                ztemp(1) = wrapToPi(ztemp(1));
                znorm(tt,idx_a) = norm(ztemp);
            end
        end
        z_l2normt = max(znorm,[],2);
        [z_infnorm, ~] = max(max(znorm));

        % Populating result structure
        invChain = struct(...
            'y',y,...
            'u',u_final,...
            'time',time,...
            'n',n,...         'ctrl',cc,...
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
        save(['results/' datetime_char '_' control{:} '_result'],'invChains')
        save(['results/' datetime_char '_' control{:} '_ctrl'],'ctrl')
    end
end

%% Plotting
% load results\200918_1523_DSS_result.mat
% load results\200918_1523_Org_result.mat
if plot_results == 1
    ChainLength = N(1);
    ChainToPlot = find(ismember(N,ChainLength));
    control = "Org";
    if control == "DSS"
        ChainToPlot = ChainToPlot + 1;
    end
    plotResults(invChains,ChainToPlot,N,control,vars);
end

%% ------------------------------------------------------------------------
function dydt = inverters(state,n,tau,eta,k,Pl,B,wref,Lambda)
   % Implements the inverter system in the absolute variables
    % setting up the simulation
    I = eye(n);
    O = zeros(n,n);
    Ii= ones(n,1);
    
    Ti = diag(ones(1,n)./tau);
    N = diag(eta);
    Ki = 1*diag(ones(1,n)./k);
    
    A = [O, I, O; O, -Ti, Ti; O, -Ki, -Ki*Lambda];   
    Aref =[O, O; Ti, Ti*N; Ki, O];    
       
    %state(1:n) = wrapToPi(state(1:n));
    Delta = repmat(state(1:n),1,n) - repmat(state(1:n)',n,1); 
    % Delta is defined just for F, as F = bij*V^2*sin(delta_i-delta_j)
    F=(B.*sin(Delta))*Ii; % 
     
    dydt = A * state + Aref * [wref;Pl] ... % free term + disturbance
           + [zeros(n,1);Ti*N*F;zeros(n,1)]; % non-linear term
    %dydt(1:n-1)= wrapToPi(dydt(1:n-1));
end


%% ------------------------------------------------------------------------
function [dydt, u] = inverters_theta(t,state,n,tau,eta,k,Pl,B,Lambda)
% Implements the inverter system in the incremental variables as per
% Schiffer et al's transformations
    % setting up the simulation
    I = eye(n);
    O = zeros(n,n);
    Ii= ones(n,1);
    V=[eye(n-1) -ones(n-1,1)];
    
    Ti = diag(ones(1,n)./tau);
    N = diag(eta);
    Ki = 1*diag(ones(1,n)./k);
    
    A = [zeros(n-1,n-1), V, 0*V;
        zeros(n,n-1), -Ti, Ti;
        zeros(n,n-1), -Ki, -Ki*Lambda];
    
    Au = [zeros(n-1,n-1), 0*V, 0*V;
        zeros(n,n-1), -Ti, Ti;
        zeros(n,n-1), -Ki, -Ki*Lambda];
    
    Delta = repmat(state(1:n),1,n) - repmat(state(1:n)',n,1);
    FF=[(-B(1:n-1,:).*sin(Delta(1:n-1,:)))*Ii; -B(n,n-1)*sin(state(n))];
 
    dydt = A * state + [zeros(n-1,1);Ti*(N*Pl);zeros(n,1)] ... % disturbance (Gill: ci)
           + [zeros(n-1,1);Ti*N*FF;zeros(n,1)]; % non-linear term  
       
    u = Au*state + [zeros(n-1,1);Ti*N*FF;zeros(n,1)];
end

function [dydt, u] = inverters_theta_dss(t,state,n,ctrl,Pl,B,Lambda)
% Implements the inverter system in the incremental variables as per
% Schiffer et al's transformations, with extra controller gains to satisfy
% Monteil et al's DSS conditions. 
%
I = eye(n);
O = zeros(n,n);
Ii= ones(n,1);
V=[eye(n-1) -ones(n-1,1)];

Ti = diag(ones(1,n)./ctrl.tau);
N = diag(ctrl.etai);
Ki = 1*diag(ones(1,n)./ctrl.ki);

A = [zeros(n-1,n-1), V, 0*V;
    zeros(n,n-1), -Ti*ctrl.Oi, Ti*ctrl.Zi;
    zeros(n,n-1), -Ki*ctrl.Li, -Ki*Lambda];

Delta = repmat(state(1:n),1,n) - repmat(state(1:n)',n,1);
FF=[(-B(1:n-1,:).*sin(Delta(1:n-1,:)))*Ii; -B(n,n-1)*sin(state(n))];

dydt = A * state + [zeros(n-1,1); Ti*(N*Pl); zeros(n,1)] ... % disturbance
           + [zeros(n-1,1); Ti*N*FF; Ki*ctrl.Yi*FF];
       
Au = [%zeros(n-1,n-1), 0*V, 0*V;
    zeros(n,n-1), -Ti*ctrl.Oi, Ti*ctrl.Zi];
    %zeros(n,n-1), -Ki*ctrl.Li, -Ki*Lambda];

u = Au*state + Ti*N*FF;%[zeros(n-1,1); Ti*N*FF; Ki*ctrl.Yi*FF];

% u = u(2);
end
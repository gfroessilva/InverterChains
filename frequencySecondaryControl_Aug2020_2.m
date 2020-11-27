%% Simulation script for inverter strings
%% All values are in Volt, Watt, ohm, second, etc
clf
clear
% Number of inverters
n = 3;
rng(1);
%% Set to 1 if random values of rated powers and loads for different tests are sought
%% r = 0 uses a saved example with n = 6 
r = 1;

if r == 0
    load inverters
end

% Per unit bases
Pbase = 1; % in kilovoltaampere [kVA]
Vbase = 415; % in volt [V]
Zbase = Vbase*Vbase/(Pbase*1000);
Ybase = 1/Zbase;

% The rated power of the inverters
PastPos = [1,2,1.5]; % from Majumder (2009)

if r == 1
    Past = datasample(PastPos,n)';
end
% per unit
Past = Past/Pbase;

% Choosing appropriate values for the controllers
tau = 10*ones(1,n);
eta = 1./Past/100;
k = 1.7 * ones(1,n);

% The active power of loads (this corresponds to the error P_i^*-P_{Li}
PlPos = [1.8,0.8,0,3]; 

%% A particular load condition
%Pl=[1.8000    1.8000    0.8000    1.8000    3.0000  0]';

if r == 1
    Pl = datasample(PlPos,n)';
end

% per unit
Pl = Pl/Pbase;

% The inductive part of the lines is a fraction of 
bij = -1/0.6; % from Majumder (2009)
% per unit
bij = bij/Ybase; %-1/3.488;  

% A particular set of values for the line impedances
%bb=[-0.4073   -8.7218  -59.8396 -130.5942  -36.5307]';

% String connection of inverters assuming no self connection
if r==1
    bb=bij*rand(1,n-1);
end

B = diag(bb,1) + diag(bb,-1); 

% Communication control in the same string
Lambda = diag([.5 ; ones(n-2,1);.5]) - 0.5*diag([1,ones(n-2,1)'],1) - 0.5*diag([ones(n-2,1)',1],-1);

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
    dw = (rand(n,1)-ones(n,1));
    save inverters Past bb dw Pl
end


% initial_state = [0*ones(n-1,1);1*rand(n,1);z0];  %% for incremental
initial_state = [0*ones(n,1);wref+dw;z0];           %% for absolute

simulation_time = 0:100;

[time,y] = ode45(@(tt,yy)inverters(yy,n,tau,eta,k,Pl,B,wref,Lambda),...
    simulation_time,initial_state);
% [time,y] = ode45(@(tt,yy)inverters_theta(yy,n,tau,eta,k,Pl,B,Lambda),...
%     simulation_time,initial_state);

Frequency_error=y(end,n+1:2*n)-wref';
Secondary_steady_state_value=y(end,end-n+1:end);

S = diag(ones(1,n))-diag((ones(1,n-1)),-1); % to calculate ang differences
%S=S(1:n-1,1:n-1);
Del = S'*y(:,1:n)'; % delta (angle differences)
V = [eye(n-1) -ones(n-1,1)];

%%
figure(1)
plot(time,wrapToPi(y(:,1:n-1)),'LineWidth',2);
plot(time,(y(:,1:n-1)),'-','LineWidth',2)
title('Angles \delta_i')
set(gca,'FontSize',16,'LineWidth',3,'Box','on')
grid on

%%
figure(2)
plot(time,wrapToPi(Del(1:n-1,:))','LineWidth',2)
plot(time,Del(1:n-1,:)','LineWidth',2)
title('Angle differences \delta_{ij}')
set(gca,'FontSize',16,'LineWidth',3,'Box','on')
grid on

%%
figure(3)
plot(time,y(:,n:2*n-1),'LineWidth',2) %% for incremental
plot(time,y(:,n+1:2*n),'LineWidth',2) %% for absolute
title('Frequency')
set(gca,'FontSize',16,'LineWidth',3,'Box','on')
grid on

%%
figure(4)
plot(time,y(:,2*n:end),'LineWidth',2) %% for incremental
plot(time,y(:,2*n+1:end),'LineWidth',2) %% for absolute
title('Secondary variable')
set(gca,'FontSize',16,'LineWidth',3,'Box','on')
grid on

%% %% ---------------------------------------------------------------------
% Implements the inverter system in the absolute variables
function dydt = inverters(state,n,tau,eta,k,Pl,B,wref,Lambda)
   
    % setting up the simulation
    I = eye(n);
    O = zeros(n,n);
    Ii= ones(n,1);
    
    Ti = diag(ones(1,n)./tau);
    N = diag(eta);
    Ki = 1*diag(ones(1,n)./k);
    
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
     
    dydt = A * state + Aref * [wref;Pl] ... % free term + disturbance
           + [zeros(n,1);Ti*N*F;zeros(n,1)]; % non-linear term
    %dydt(1:n-1)= wrapToPi(dydt(1:n-1));
end


%%%% ----------------------------------------------------------------------
% Implements the inverter system in the incremental variables as per
% Schiffer et al's transformations
function dydt = inverters_theta(state,n,tau,eta,k,Pl,B,Lambda,t)

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

    Delta = repmat(state(1:n-1),1,n-1) - repmat(state(1:n-1)',n-1,1);
    FF=[(-B(1:n-1,1:n-1).*sin(Delta(1:n-1,:)))*Ii(1:n-1); -B(n,n-1)*sin(state(n-1))];
 
    dydt = A * state + [zeros(n-1,1);Ti*(N*Pl);zeros(n,1)] ... % disturbance
           + [zeros(n-1,1);Ti*N*FF;zeros(n,1)]; % non-linear term  
end

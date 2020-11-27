syms t              % time
syms z0             % initial condition
syms k1 k2 k3 k4    % constants based on desired voltage and reactive power
                    % rating and time constant of power sensor

n = 3;
z0 = 12245;
tau = 10;
etap = 1e-3;
k1 = 1/tau;
Vd = 500*ones(n,1); Qd = 5*ones(n,1);
k2 = sum(Vd/etap+Qd);
k3 = 1e-2;
k4 = sqrt(4*k2*k3+k1^2);
% solution from Schiffer 2014                    
num = (2*k2*(-1+exp(k4*t))+z0*(k1*(1-exp(k4*t))+k4*(1+exp(k4*t))));
den = (k1*(-1+exp(k4*t))+k4*(1+exp(k4*t))+2*k3*z0*(-1+exp(k4*t)));
zt = num/den;

figure(1)
fplot(zt,[0 .3])
limzt = (2*k2+z0*(-k1+k4)) / ...
        (k1+k4+2*k3*z0);

% time derivative of solution
dtzt = diff(zt,t);

figure(2)
fplot(dtzt,[0 0.3])
% solution when dtzt equals 0
% t_sol = solve(dtzt<=0.1,t);

% evaluating zt when t=t_sol
% zt_sol = subs(zt,t,t_sol.z0);


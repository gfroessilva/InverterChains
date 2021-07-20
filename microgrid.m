function [dydt, y] = microgrid(t, state, sys, ctrl, flag)

persistent tlast msg
if isempty(tlast)
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
% delta = wrapTo2Pi(state(1:n));
delta = (state(1:n));
% state(1:n) = wrapTo2Pi(state(1:n));
% omega = state(n+1:2*n);
V = state(2*n+1:3*n);
% xi = state(3*n+1:4*n);
% zeta = state(4*n+1:5*n);

% Provides Vij in diagonal+-1 and Vi^2 in main diagonal
VV = repmat(V,1,n).*repmat(V',n,1);

delta_ij = repmat(delta,1,n) - repmat(delta',n,1);
% delta_ij = wrapTo2Pi(repmat(delta,1,n) - repmat(delta',n,1));

if t>flag.SecOnT
    flag.Sec = 1*flag.SecOn;
else
    flag.Sec = 0;
end

% Powers (injections) used in controller
P = (abs(Bij).*sin(delta_ij)).*VV*Iv;
Q = (VV.*abs(Bii) - (abs(Bij).*cos(delta_ij)).*VV)*Iv;

% Plr(1) = Plr(1).*1.05;

if t>flag.DistOnT && flag.Dist && t<flag.DistOffT
    Plr(1) = Plr(1).*1.2; %
    Plr(2:2:end) = Plr(2:2:end).*sys.dist(2:2:end);
    dist = zeros(n,1);
else
    dist = zeros(n,1);
end

% Powers output "measurement"
y.Pi = P;
y.Qi = Q;
y.P = Plr + P; % Pi
y.Q = Qlr + Q;

% Y = imag(Yl) + Bij;
% y.Q = diag(V)*Y*V;

% Reactive Power Sharing
QDelta = repmat(y.Q./(Qast),1,n) - repmat((y.Q./(Qast))',n,1);

ref = [omegaast; Vast; Past; Qast];

if ctrl.control == "C1"
    tau =   ctrl.ctrlgain(1);  % + makes the system faster???
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
    
    switch ctrl.opt
        case 0
        case 1
        case 2
            A = [...
                O I O O O;
                [O -T*Oo O flag.Sec*T*Oxi O]*flag.PriF;
                [O O -T O flag.Sec*T*Vzeta]*flag.PriV;
                [O -K*Xio O -K*ctrl.G O]*flag.Sec; 
                [O O -Kappa*ZetaV O -Kappa*ctrl.H]*flag.Sec*flag.PriV];
            
            Aref = [O O O O;
                [T*Oo O T*Np O]*flag.PriF;
                [O T O T*Nq]*flag.PriV;
                [K*Xio O K*Xid O]*flag.Sec;
                [O Kappa*ZetaV O O]*flag.Sec*flag.PriV];
            
            PQ = [Ov % deltadot
                -Np.*T*(Plr + P)*flag.PriF          % omegadot -= etaP*P
                -Nq.*T*(Qlr + Q)*flag.PriV          % Vdot -= etaQ*Q
                -Xid.*K*(sin(delta_ij)*Iv)*flag.Sec %(Alterante 2)
                -ctrl.H./kzeta.*QDelta*Iv*flag.Sec*flag.PriV;
                ];
        case 21
        case 22
            A = [...
                O I O O O;
                [O -T*Oo O flag.Sec*T*Oxi O]*flag.PriF;
                [O O -T O flag.Sec*T*Vzeta]*flag.PriV;
                [O -K*Xio O -K*ctrl.G O]*flag.Sec; 
                [O O -Kappa*ZetaV O -Kappa*ctrl.H]*flag.Sec*flag.PriV];
            
            Aref = [O O O O;
                [T*Oo O T*Np O]*flag.PriF;
                [O T O T*Nq]*flag.PriV;
                [K*Xio O K*Xid O]*flag.Sec;
                [O Kappa*ZetaV O O]*flag.Sec*flag.PriV];
            
            PQ = [Ov % deltadot
                -Np.*T*(Plr + P)*flag.PriF          % omegadot -= etaP*P
                -Nq.*T*(Qlr + Q)*flag.PriV          % Vdot -= etaQ*Q
                -Xid.*K*((Plr + P)+sys.Comm.*sin(delta_ij)*Iv)*flag.Sec %(Alterante 2)
                -ctrl.H./kzeta.*QDelta*Iv*flag.Sec*flag.PriV;
                ];
        case 23
            A = [...
                O I O O O;
                [O -T*Oo O flag.Sec*T*Oxi O]*flag.PriF;
                [O O -T O flag.Sec*T*Vzeta]*flag.PriV;
                [-K*ctrl.L -K*Xio O -K*ctrl.G O]*flag.Sec; % Alternate CL 2.3 (ctlr.L)
                %                 [O -K*Xio O -K*ctrl.G O]*flag.Sec; % Alternate CL 2.4
                [O O -Kappa*ZetaV O -Kappa*ctrl.H]*flag.Sec*flag.PriV];
            Aref = [O O O O;
                [T*Oo O T*Np O]*flag.PriF;
                [O T O T*Nq]*flag.PriV;
                [K*Xio O O O]*flag.Sec;
                [O Kappa*ZetaV O O]*flag.Sec*flag.PriV];
            
            PQ = [Ov % deltadot
                -Np.*T*(Plr + P)*flag.PriF          % omegadot -= etaP*P
                -Nq.*T*(Qlr + Q)*flag.PriV          % Vdot -= etaQ*Q
                % xidot:
                %                 -Xid.*K*(sys.Comm.*sin(delta_ij)*Iv)*flag.Sec %(Alterante 2)
                %         (-Xid.*K*(Plr + P))*flag.Sec % (Alternate 2.1)
                %         (-Xid.*K*(sys.Comm.*delta_ij*Iv)-Xid.*K*(Plr + P))*flag.Sec % (Alternate 2.2)
                Ov; % (Alterante CL 2.3)
                %         -ctrl.L*K*(Plr + P); % (Alternate CL 2.4)
                -ctrl.H./kzeta.*QDelta*Iv*flag.Sec*flag.PriV;
                ];
        case 24
            A = [...
                O I O O O; %OK
                [O -T*Oo O flag.Sec*T*Oxi O]*flag.PriF; %OK
                [O O -T O flag.Sec*T*Vzeta]*flag.PriV; 
                [O -K*Xio O -K*ctrl.G O]*flag.Sec; % Alternate CL 2.4 %OK
                [O O -Kappa*ZetaV O -Kappa*ctrl.H]*flag.Sec*flag.PriV];
            Aref = [O O O O;
                [T*Oo O T*Np O]*flag.PriF; %OK
                [O T O T*Nq]*flag.PriV;
                [K*Xio O O O]*flag.Sec;
                [O Kappa*ZetaV O O]*flag.Sec*flag.PriV];
            
            PQ = [Ov % deltadot
                -Np.*T*(Plr + P)*flag.PriF          % omegadot -= etaP*P
                -Nq.*T*(Qlr + Q)*flag.PriV          % Vdot -= etaQ*Q                
                -ctrl.L*K*(Plr + P); % (Alternate CL 2.4)
                -ctrl.H./kzeta.*QDelta*Iv*flag.Sec*flag.PriV;
                ];
    end
    
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
        -ctrl.etaP.*(Plr + P)*flag.PriF          % omegadot -= etaP*P
        -ctrl.etaQ.*(Qlr + Q)*flag.PriV          % Vdot -= etaQ*Q
        Ov
        -ctrl.H./ctrl.kappa.*QDelta*Iv*flag.Sec*flag.PriV;
        ];
    
    Dist = [Ov; ctrl.etaP.*dist; Ov; Ov; Ov];
end

dydt = A*state + ... % linear dynamics
    Aref*ref + ...     % ref (ast)
    PQ + ...        % nonlinear dynamics (powers)
    Dist;           % Disturbance (power)

end
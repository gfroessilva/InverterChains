function Jii_L_tilde = JiiL_Vcte(syspar,T,comm,tau,etai,ki,Xio,Xid,Oo,Oxi,ctrl_opt)
gim1 = comm(1);
gip1 = comm(2);
try
    lim1 = comm(3);
    lip1 = comm(4);
catch
end
alpha = syspar.alpha;

switch ctrl_opt
    case 1
        Jii_L = [0 1 0;
            -etai*(2*alpha)/tau -Oo/tau Oxi/tau;
            0 -Xio/ki -(gim1 + gip1)/ki];
    case 2
        Jii_L = [0 1 0;
            -etai*(2*alpha)/tau -Oo/tau Oxi/tau;
            -Xid*(2*alpha)/ki -Xio/ki -(gim1 + gip1)/ki];        
    case 21
        Jii_L = [0 1 0;
            -etai*(2*alpha)/tau -Oo/tau Oxi/tau;
            -Xid*(2*alpha)/ki -Xio/ki -(gim1 + gip1)/ki];                
    case 22
        Jii_L = [0 1 0;
            -etai*(2*alpha)/tau -Oo/tau Oxi/tau;
            -Xid*(2*alpha)/ki -Xio/ki -(gim1 + gip1)/ki];                
    case 23
        Jii_L = [0 1 0;
            -etai*(2*alpha)/tau -Oo/tau Oxi/tau;
            -(lip1 + lim1)/ki  -Xio/ki -(gim1 + gip1)/ki];
    case 24
        Jii_L = [0 1 0;
            -etai*(2*alpha)/tau -Oo/tau Oxi/tau;
            -3*alpha*(lip1 + lim1)/ki -Xio/ki -(gim1 + gip1)/ki];
end

Jii_L_tilde = T*Jii_L/T;
end
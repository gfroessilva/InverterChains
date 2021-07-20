function Jip1_U_tilde = Jip1U_Vcte(syspar,T,comm,tau,etai,ki,Xio,Xid,Oo,Oxi,ctrl_opt)
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
        Jip1_U = [0 0 0;
            etai*alpha/tau 0 0;
            0 0 gip1/ki];
    case 2
        Jip1_U = [0 0 0;
            etai*alpha/tau 0 0;
            Xid*alpha/ki 0 gip1/ki];
    case 21
        Jip1_U = [0 0 0;
            etai*alpha/tau 0 0;
            Xid*alpha/ki 0 gip1/ki];        
    case 22
        Jip1_U = [0 0 0;
            etai*alpha/tau 0 0;
            Xid*alpha/ki 0 gip1/ki];        
    case 23
        Jip1_U = [0 0 0;
            etai*alpha/tau 0 0;
            (Xid*alpha + lip1)/ki 0 gip1/ki];
    case 24
        Jip1_U = [0 0 0;
            etai*alpha/tau 0 0;
            2*(lip1*alpha)/ki 0 gip1/ki];
end

Jip1_U_tilde = T*Jip1_U/T;
end
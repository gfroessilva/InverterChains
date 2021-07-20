function Jii_U_tilde = JiiU_Vcte(syspar,T,comm,tau,etai,ki,Xio,Xid,Oo,Oxi,ctrl_opt)
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
        Jii_U = [0 1 0;
            0 -Oo/tau Oxi/tau;
            0 -Xio/ki -(gim1 + gip1)/ki];
    case 2
        Jii_U = [0 1 0;
            0 -Oo/tau Oxi/tau;
            0 -Xio/ki -(gim1 + gip1)/ki];        
    case 21
        Jii_U = [0 1 0;
            0 -Oo/tau Oxi/tau;
            0 -Xio/ki -(gim1 + gip1)/ki];         
    case 22
        Jii_U = [0 1 0;
            0 -Oo/tau Oxi/tau;
            0 -Xio/ki -(gim1 + gip1)/ki];         
    case 23
        Jii_U = [0 1 0;
            0 -Oo/tau Oxi/tau;
            -(lim1+lip1)/ki -Xio/ki -(gim1 + gip1)/ki];
    case 24
        Jii_U = [0 1 0;
            0 -Oo/tau Oxi/tau;
            0 -Xio/ki -(gim1 + gip1)/ki];
end

Jii_U_tilde = T*Jii_U/T;
end
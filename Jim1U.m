function Jim1_U_tilde = Jim1U(Gamma,syspar,comm,ctrlgain,ctrlxtra,Vbar)
tau = ctrlgain(1);
np = ctrlgain(2);
nq = ctrlgain(3);
kxi = ctrlgain(4);
kzeta = ctrlgain(5);
l = ctrlxtra(1);
m = ctrlxtra(2);
x = ctrlxtra(3);
y = ctrlxtra(4);
gim1 = comm(1);
gip1 = comm(2);
him1 = comm(3);
hip1 = comm(4);
B = syspar.B;

Jim1_U = [0 0 0 0 0;
    0 0 0 0 0;
    nq*B*Vbar^2/tau 0 0 0 0;
    0 0 0 gim1/kxi 0;
    0 0 0 0 him1/kzeta];
Jim1_U_tilde = Gamma*Jim1_U/Gamma;
end
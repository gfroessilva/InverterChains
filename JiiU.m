function Jii_U_tilde = JiiU(Gamma,syspar,comm,ctrlgain,ctrlxtra,Vbar)
tau = ctrlgain(1);
np = ctrlgain(2);
nq = ctrlgain(3);
kxi = ctrlgain(4);
kzeta = ctrlgain(5);
l = ctrlxtra(1);
m = ctrlxtra(2);
o = ctrlxtra(3);
x = ctrlxtra(4);
y = ctrlxtra(5);
gim1 = comm(1);
gip1 = comm(2);
him1 = comm(3);
hip1 = comm(4);
B = syspar.B;

Jii_U = [0              1       0       0               0;
    0                   -1/tau  0       y/tau           0;
    nq/tau*2*B*Vbar^2   0       -o/tau  0               m/tau;
    0                   -x/kxi  0       -(gip1+gim1)/kxi 0;
    0                   0       -l/kzeta 0          -(hip1+him1)/kzeta];
Jii_U_tilde = Gamma*Jii_U/Gamma;
end
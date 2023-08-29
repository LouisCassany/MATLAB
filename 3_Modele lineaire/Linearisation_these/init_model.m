function [meal_signal, params_model, CI, basal] = init_model(patient,sim_time, meal)
% meal = [meal_time, meal_duration, meal_size in grams];

load patients_data

D = meal(3)*1000;

Lstruttura = patients(patient);

basal = Lstruttura.u2ss*Lstruttura.BW; % Patient basal in pmol/min

%% Constantes du modèle

kgri = Lstruttura.kmax;
kabs = Lstruttura.kabs;
BW = Lstruttura.BW;
kmin = Lstruttura.kmin;
kmax = Lstruttura.kmax;
b = Lstruttura.b;
c = Lstruttura.d;
alpha = 5/(2*D*(1-b));      %Dalla man 2006
beta = 5/(2*D*c);           %Dalla man 2006
f = Lstruttura.f;
k1 = Lstruttura.k1;
k2 = Lstruttura.k2;
Vg = Lstruttura.Vg;
ke1 = Lstruttura.ke1;
ke2 = Lstruttura.ke2;
m1 = Lstruttura.m1;
m2 = Lstruttura.m2;
m4 = Lstruttura.m4;
m5 = Lstruttura.m5;
m30 = Lstruttura.m30;
Ilb = Lstruttura.Ilb;
Ipb = Lstruttura.Ipb;
Sb = m30*Ilb*m4*Ipb;
Heb = Lstruttura.HEb;       
m6 = m5*Sb+Heb;
m3 = m30;
VI = Lstruttura.Vi;
Fcns = Lstruttura.Fsnc;
Vmo = Lstruttura.Vm0;
Vmx = Lstruttura.Vmx;
r1 = Lstruttura.r1;
r2 = Lstruttura.r2;
Kmo = Lstruttura.Km0;
p2U = Lstruttura.p2u;
Ib = Lstruttura.Ib;
Gb = Lstruttura.Gb;
Gth = Lstruttura.Gth;
kp1 = Lstruttura.kp1;
kp2 = Lstruttura.kp2;
kp3 = Lstruttura.kp3;
kh = Lstruttura.kXGn;
ki = Lstruttura.ki;
Hb =  Lstruttura.Gnb;
zeta = Lstruttura.kcounter;
ka1 = Lstruttura.ka1;
ka2 = Lstruttura.ka2;
kd = Lstruttura.kd;
Ts = 1/Lstruttura.ksc;
rho = Lstruttura.alfaG;
sigma2 = Lstruttura.kGSRs;     
sigma = 6.67;                % Analyse paramétrique Louis
delta = Lstruttura.kGSRd;
n = Lstruttura.k01g;
kh1 = Lstruttura.SQgluc_k1;
kh2 = Lstruttura.SQgluc_kc1;
Hinf = 0;
kh3 = Lstruttura.SQgluc_k2;
SRHb = Hb*n;                 % From 2014

params_model=[kgri kabs BW kmin kmax b c alpha beta f k1 k2 Vg ke1 ke2 m1 m2 m6 m3 m4 VI Fcns Vmo Vmx r1 r2 Kmo p2U Ib Gb Gth kp1 kp2 kp3 kh ki Hb zeta ka1 ka2 kd Ts rho sigma2 sigma delta n kh1 kh2 Hinf kh3 SRHb D].';

%% Conditions initiales

Gpb = Lstruttura.Gpb;
Gtb = Lstruttura.Gtb;
Ipb = Lstruttura.Ipb;
Ilb = Lstruttura.Ilb;
Isc1ss = Lstruttura.isc1ss;
Isc2ss = Lstruttura.isc2ss;
Gsb = Lstruttura.Gpb;
Hb =  Lstruttura.Gnb;
SRHb = Hb*n;  
Hsc1b = 0;
Hsc2b = 0;

CI = [0 0 0 Gpb Gtb Ipb 0 Ib Ib Ilb Isc1ss Isc2ss Gsb Hb 0 SRHb Hsc1b Hsc2b].';
CI = CI(1:13);

%% Meal

meal_duration = meal(2);
meal_time = meal(1);
k_meal = 1/15;

meal_signal = zeros(1,sim_time+1).';

for index=1:meal_duration
    meal_signal(meal_time+index)=D;
end

meal_signal = [(0:sim_time).' meal_signal*k_meal];

end


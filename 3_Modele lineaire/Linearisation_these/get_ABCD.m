function [sys] = get_ABCD(params_model, basal, x, x_equil,index)

syms Qsto1 Qsto2 Qgut Gp Gt Ip X I_ XL Il Isc1 Isc2 Gs

syms meal_symbolic Ins

syms risk E XH_dot SRHs_dot EGP Ra Uii Uid Rai G_dot SRHd SRH Qsto kempt

%% Constantes du modèle
kgri = params_model(1);
kabs = params_model(2);
BW = params_model(3);
kmin = params_model(4);
kmax = params_model(5);
b = params_model(6);
c = params_model(7);
alpha = params_model(8);
beta = params_model(9);
f = params_model(10);
k1 = params_model(11);
k2 = params_model(12);
Vg = params_model(13);
ke1 = params_model(14);
ke2 = params_model(15);
m1 = params_model(16);
m2 = params_model(17);
m6 = params_model(18);
m3 = params_model(19);
m4 = params_model(20);
VI = params_model(21);
Fcns = params_model(22);
Vmo =params_model(23);
Vmx = params_model(24);
r1 = params_model(25);
r2 = params_model(26);
Kmo = params_model(27);
p2U = params_model(28);
Ib = params_model(29);
Gb = params_model(30);
Gth = params_model(31);
kp1 =params_model(32);
kp2 =params_model(33);
kp3 = params_model(34);
kh = params_model(35);  
ki = params_model(36);              
Hb =  params_model(37);            
zeta = params_model(38);           
ka1 = params_model(39);
ka2 = params_model(40);
kd = params_model(41);
Ts = params_model(42);              
rho = params_model(43);          
sigma2 = params_model(44);        
sigma = params_model(45);      
delta =params_model(46);     
n = params_model(47);             
kh1 = params_model(48);
kh2 = params_model(49);
Hinf = params_model(50);
kh3 = params_model(51);
SRHb = params_model(52);
D = params_model(53);

%% Modèle NL ordre 13

G = Gp/Vg;
I = Ip/VI;

if x(4)/Vg >= Gb
    risk = 0;
elseif ((Gth <= x(4)/Vg) && (x(4)/Vg < Gb))
    risk = 10*(log10(G/Gb)^r2)^2;
else
    risk = 10*(log10(Gth/Gb)^r2)^2;
end


if x(4) > ke2
    E = ke1*(Gp-ke2);
else
    E = 0;
end

if(kp1-kp2*(x(4)/Vg)-kp3*x(9) < 0)
    EGP = 0;
else
    EGP = kp1-kp2*Gp-kp3*XL;
end

Ra = (f*kabs*Qgut)/(BW);
Uii = Fcns;
Uid = ((Vmo+Vmx*X*(1+r1*risk))*Gt)/(Kmo+Gt);
Rai = ka1*Isc1+ka2*Isc2;

Qsto = Qsto1+Qsto2; 

% kempt = kmin + ((kmax-kmin)/2) * (tanh(alpha*(Qsto-b*D))-tanh(beta*(Qsto-c*D))+2);
kempt = kmin + ((kmax-kmin)/2) * (tanh(alpha*((x(1)+x(2))-b*D))-tanh(beta*((x(1)+x(2))-c*D))+2);

x_dot = [-kgri*Qsto1+meal_symbolic ;... % Qsto1_dot
    -kempt*Qsto2+kgri*Qsto1 ;...        % Qsto2_dot
    -kabs*Qgut + kempt*Qsto2 ;...       % Qgut_dot
    EGP+Ra-Uii-E-k1*Gp+k2*Gt;...        % Gp_dot
    -Uid+k1*Gp-k2*Gt;...                % Gt_dot
    -(m2+m4)*Ip+m1*Il+Rai;...           % Ip_dot
    -p2U*X+p2U*(I-Ib); ...              % X_dot
    -ki*(I_-I); ...                     % I__dot
    -ki*(XL-I_); ...                    % XL_dot
    -(m1+m3)*Il+m2*Ip ;...              % Il_dot
    -(kd+ka1)*Isc1+Ins/BW               % Isc1_dot
    kd*Isc1-ka2*Isc2;...                % Isc2_dot
    -(1/Ts)*Gs+(1/Ts)*Gp]; ...          % Gs_dot
    
x_symbolic = [Qsto1 Qsto2 Qgut Gp Gt Ip X I_ XL Il Isc1 Isc2 Gs].';
u = [meal_symbolic Ins].';
% y = G; % x4/Vg
y = Gs/Vg; % x13/Vg

x_dot = eval(x_dot);

A = jacobian(x_dot,x_symbolic);
B = jacobian(x_dot,u);
C = jacobian(y,x_symbolic);
D = jacobian(y,u);

% x_equil = x;

Qsto1 = x_equil(1);
Qsto2 = x_equil(2);
Qgut = x_equil(3);
Gp = x_equil(4);
Gt = x_equil(5);
Ip = x_equil(6);
X = x_equil(7);
I_ = x_equil(8);
XL = x_equil(9);
Il = x_equil(10);
Isc1 = x_equil(11);
Isc2 = x_equil(12);
Gs = x_equil(13);

meal_symbolic = 0;

if(index == 1)
    Ins = 0;
else
    Ins = basal;
end

sys.A = eval(A);
sys.B = eval(B);
sys.C = eval(C);
sys.D = eval(D);

end


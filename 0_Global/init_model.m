function [r, params] = init_model()
    r = 50; % cable length (m)
    m = 4; % kite's mass (kg)
    g_k = 0.16652; % turn rate
    M_k = -3.6438; % mass distribution effect (+ : "piqueur", -:"redresseur")
    A_k = 15; % kite surface (m^2)
    rho_air = 1.225; % air density (kg/m^3)
    g = 9.81; % gravity acceleration (m/s^2)
    K_i = 0.4474; % border/choquer coeficient
    alpha_i0 = 30*pi/180; % initial value for alpha (°)
    C_L = 1.2; % lift force coef
    C_D = 0.2; % drag force coef
    C_T = 0.2; % trans force coef
    params = [r; m; g_k; M_k; A_k; rho_air; g; K_i; alpha_i0; C_D; C_L; C_T];
end

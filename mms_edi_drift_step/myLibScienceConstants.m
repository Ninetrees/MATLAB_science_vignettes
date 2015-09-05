% myLibScienceConstants

c      =  2.99792458e+08;  % speed of light (SI)
mass_e =  9.10938188e-31;  % electron mass, kg (SI)
q      =  1.60217650e-19;  % coulomb (SI)
mV2V   =  1.0e-3;          % mV > V (SI)
nT2T   =  1.0e-9;          % nT > T (SI)
C2V2T  =  mV2V / nT2T;     % Combining constants to save flops; potentially used 100Ks of times

twoPi     = 2.0 * pi;
halfPi    = pi / 2.0;
quarterPi = pi / 4.0;

deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;

Re = 6371; % kilometers

q_over_m = -q / mass_e;

% Non-relativistic: KE = mv^2 / 2 ~> v_nr = sqrt (2.0 * 500 eV / mass_e)
v_500eV_electron = 13262052.; % m/s
v_1keV_electron  = 18755373.; % m/s
% Relativistic
v_500eV_electron = 13252328.; % m/s, difference of 0.073% ~> 0.99926678 * nr
v_1keV_electron  = 18727897.; % m/s, difference of 0.147% ~> 0.99853503 * nr

% myLibScienceConstants

% All these values are assumed to be precise (no uncertainty), unless otherwise noted.
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

Re = 6371000.0; % meters
% m_Earth =
% G =
GM = 3.986e+14; % m^3 / s^2

q_over_m = -q / mass_e;
q_over_mass_e_nT2T = q * nT2T / mass_e;

% See 'Gyrotime from B calculations * M15.pdf'.
% Non-relativistic: KE = mv^2 / 2 ~> v_nr = sqrt (2.0 * 500 eV / mass_e)
v_500eV_electron = 13262052.; % m/s
v_1keV_electron  = 18755373.; % m/s
% Relativistic
vr_500eV_electron = 13252328.; % m/s, difference of 0.073% ~> 0.99926678 * nr
vr_1keV_electron  = 18727897.; % m/s, difference of 0.147% ~> 0.99853503 * nr

e_gamma_500eV     = 1 / sqrt (1 - (vr_500eV_electron / c)^2);
e_gamma_1keV      = 1 / sqrt (1 - (vr_1keV_electron  / c)^2);
nT2sr_500eV       = e_gamma_500eV * twoPi * mass_e / q * 1.0e9; % T ~> nT
nT2sr_1keV        = e_gamma_1keV  * twoPi * mass_e / q * 1.0e9; % T ~> nT
nT2usr_500eV      = e_gamma_500eV * twoPi * mass_e / q * 1.0e9 * 1.0e6; % T ~> nT, s ~> us
nT2usr_1keV       = e_gamma_1keV  * twoPi * mass_e / q * 1.0e9 * 1.0e6; % T ~> nT, s ~> us

%written by M. Jabas 08/2009 (approx.)
%this database contains the material properties of the
%selected materials and calculates the Q-matrix of the desired
%material
%input: laminatdescription: alph, mat_name, fibrevolume pfi_f

function[Q_lam_t,Q_lam_c,T,rho,zul_t,zul_c,th_skin] = mat_data(angle)
phi_f=0.35;
%density
rho_f_cfk = 1750;       %carbon fibre density
rho_f_gfk = 2550;       %glas fiber density
rho_m = 1140;       %matrix density
rho_mylar = 1390;   %mylar density
th_skin = 12e-6;
rho(1:2) = rho_f_cfk * phi_f + (1-phi_f)* rho_m;
rho(3:4) = rho_f_gfk * phi_f + (1-phi_f)* rho_m;
rho(5) =30;
rho(6) = rho_mylar;
%material properties for tension
%stiffnes for phi_f = 0.35
phi_f_old = 0.35;
%[MPA]         cfk-gew cfk-ud gfk-gew gfk-ud  core1(conticellCC31) MylarADS
%              from paper 48(gauge) 12mu_m
E1t_o  = 1e6.* [35853; 73257; 16850;  31435;  15.5; 4894;];
E2t_o  = 1e6.* [35853; 9744;  16850;  8949;   15.5; 4894;];
v12_o  =       [0.08;  0.25;   0.12;  0.26;   0.22; 0;];
G12_o  = 1e6.* [1294;  1570;   1678;  1299;   3.7; 4894/2;];

%strength patameters
sigma_1_t_zul_o = 1e6.* [376;  914; 291;  554; 0.24; 103;];
sigma_2_t_zul_o = 1e6.* [376;  115; 291;  30;  0.24; 103;];
tau_12_zul_o    = 1e6.* [63;   67;   48;  46;  1.17; NaN;];


%change in fibrevolume:
E1t = phi_f/phi_f_old.* E1t_o(1:4);
E2t = phi_f/phi_f_old.* E2t_o(1:4);
G12 = phi_f/phi_f_old.* G12_o(1:4);
v12 = v12_o;
sigma_1_t_zul = phi_f/phi_f_old.* sigma_1_t_zul_o(1:4);
sigma_2_t_zul = phi_f/phi_f_old.* sigma_2_t_zul_o(1:4);
tau_12_zul    = phi_f/phi_f_old.* tau_12_zul_o(1:4);


zul_t = [sigma_1_t_zul sigma_2_t_zul tau_12_zul;...
    sigma_1_t_zul_o(5) sigma_2_t_zul_o(5) tau_12_zul_o(5);...
    sigma_1_t_zul_o(6) sigma_2_t_zul_o(6) tau_12_zul_o(6)];


for j=1:1:length(E1t)

%calculation of Q-matrix in cordinate system of the material

Q_lam_t(:,:,j) = [E1t(j)/(1-v12(j)^2) v12(j)*E2t(j)/(1-v12(j)^2) 0;...
	 v12(j)*E2t(j)/(1-v12(j)^2) E2t(j)/(1-v12(j)^2) 0; 0 0 G12(j);];

%Q-matrix in global cordinatesystem
%Q_glob_t1 = T.*Q_lam_t1.*T';

end

Q_lam_t(:,:,5) = [E1t_o(5) 0 0;...
	 0 E2t_o(5) 0; 0 0 G12_o(5);];
Q_lam_t(:,:,6) = [E1t_o(6) 0 0;...
	 0 E2t_o(6) 0; 0 0 G12_o(6);];

%===============================================
%===============================================
%material properties for compression
%stiffnes for phi_f = 0.35
phi_f_old = 0.35;
      %        cfk-gew cdk-ud gfk-gew gfk-ud  core1 MylarADS
E1c_o  = 1e6.* [34683; 79265; 16976; 24235;   15.5; 4894;];
E2c_o  = 1e6.* [34683; 5271;  16976; 8526;    15.5; 4894;];
v12    =       [0.08;  0.25;   0.12;  0.26;   0.22; NaN;];
G12_o  = 1e6.* [1294;  1570;   1678;  1299;   3.7;  NaN;];

%strenth patameters
sigma_1_c_zul_o = 1e6.* [263;  368; 224;  337; 0.24; 103;];
sigma_2_c_zul_o = 1e6.* [263;  128; 224;  120; 0.24; 103;];
tau_12_zul_o    = 1e6.* [63;   67;   48;  46;  1.17; NaN;];





%change in fibrevolume:
E1c = phi_f/phi_f_old.* E1c_o(1:4);
E2c = phi_f/phi_f_old.* E2c_o(1:4);
G12 = phi_f/phi_f_old.* G12_o(1:4);
v12 = v12;
sigma_1_c_zul = phi_f/phi_f_old.* sigma_1_c_zul_o(1:4);
sigma_2_c_zul = phi_f/phi_f_old.* sigma_2_c_zul_o(1:4);
tau_12_zul    = phi_f/phi_f_old.* tau_12_zul_o(1:4);

zul_c = [sigma_1_c_zul sigma_2_c_zul tau_12_zul;...
    sigma_1_c_zul_o(5) sigma_2_c_zul_o(5) tau_12_zul_o(5);...
    sigma_1_c_zul_o(6) sigma_2_c_zul_o(6) tau_12_zul_o(6)];



for j=1:1:length(E1c)

%calculation of Q-matrix in cordinatesystem of the material
Q_lam_c(:,:,j) = [E1c(j)/(1-v12(j)^2) v12(j)*E2c(j)/(1-v12(j)^2) 0;...
	 v12(j)*E2c(j)/(1-v12(j)^2) E2c(j)/(1-v12(j)^2) 0; 0 0 G12(j);];



%Q-matrix in global cordinatesystem
%Q_glob_c1 = T.*Q_lam_c1.*T';

end

Q_lam_c(:,:,5) = [E1c_o(5) 0 0;...
	 0 E2c_o(5) 0; 0 0 G12_o(5);];
Q_lam_c(:,:,6) = [E1c_o(6) 0 0;...
	 0 E2c_o(6) 0; 0 0 G12_o(6);];

T = [cos(angle)^2 sin(angle)^2 -2*sin(angle)*cos(angle);...
	sin(angle)^2 cos(angle)^2 2*sin(angle)*cos(angle);...
	sin(angle)*cos(angle) -sin(angle)*cos(angle) cos(angle)^2-sin(angle)^2;];

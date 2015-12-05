% Structure Design Module
% Stefan Leutenegger
% 12/2009
%
% Inputs: 
%    b:        [m]  wing span
%    AR:       [-]  aspect ratio
%    m_pld:    [kg] point mass
%    m_distr:  [kg] distributed mass (e.g. battery)
%    n:        [-]  number of distributed propulsion units
% =========================================================================
function [m_struct,masses,thicknesses,velocities,polar]=StructureDesigner(b,AR,m_pld,m_distr,n,ismars)
mult = 0.95; %set to meet Siebenmanns calculations for 5.4m/AR19.6/mbat4kg plane with
%an additional 20% margin(meeting of 20120607). However, this calculation
%assumed a rib-Structure!!

N  = 24;  % discretization
n_propulsion = 0.6; %propulsion group overall efficiency
if n==0 % e.g. Human Powered, no propulsion applied
    n=1;
    k_prop=0;
else
    k_prop = 1.1e-3;
end

c  = b/AR;
A  = c*b;

m_propulsion_old = 0.2; %[kg] propulsion mass
if(ismars~=1)  
    g=9.806;
    rho=1.225;
    mu=16e-6;
else
    g=3.711;
    rho=.699 / (.1921 * 200); %sort of worst case max
    mu=15e-6;
end

% thicknesses - obsolete
% t_f_wing  = 0.001;  % wing flange thickness
% t_s_wing  = 0.0001;  % wing shell thickness (overall)
% t_sp_wing = 0.0001; % wing spar thickness (overall)
% t_spsw_wing = 0.001; % wing spar sandwich thickness
% mat_s_wing = 1; % 1=f_cfrp, 3=f_gfrp
% t_sw_wing = 0.002; % wing sandwich thickness
% t_fus     = 0.001; % fuselage shell thickness
% t_f_hor  = 0.01*0.001;  % hor flange thickness
% t_s_hor  = 0.0003;  % hor shell thickness (overall)
% t_sp_hor = 0.0001; % hor spar thickness (overall)
% t_spsw_hor = 0.001; % hor spar sandwich thickness
% mat_s_hor = 1; % 1=f_cfrp, 3=f_gfrp
% t_sw_hor = 0.002; % hor sandwich thickness
% t_f_fin  = 0.01*0.001;  % fin flange thickness
% t_s_fin  = 0.0003;  % fin shell thickness (overall)
% t_sp_fin = 0.0001; % fin spar thickness (overall)
% t_spsw_fin = 0.001; % fin spar sandwich thickness
% mat_s_fin = 1; % 1=f_cfrp, 3=f_gfrp
% t_sw_fin = 0.002; % fin sandwich thickness

% aerodynamics look-up
[ReList,wingPolar,wingPolarp20,wingPolarm20,tailPolarp30,A0,A0_tail,s0, s0_tail,s0_LE,s0_TE,Iy_LE,Ix_LE,Iy_TE,Ix_TE]=init;

% get the material data 
% cfk-gew cfk-ud gfk-gew gfk-ud core1(conticellCC31) MylarADS
[Q_lam_t,Q_lam_c,T45,rho_mat,zul_t,zul_c,th_skin] = mat_data(45/180*pi);
Q_mat_crfp=T45*Q_lam_c(:,:,1)*T45';
Q_mat_grfp=T45*Q_lam_c(:,:,3)*T45';

%% geometry
% discretization
y   = b/2/(N/2)/2:b/2/(N/2):b/2-b/2/(N/2)/2;
y2  = 0:b/2/(N/2):b/2;
delta = b/2/(N/2);
for j=1:size(y,2)
    Y(j,:)=y-y2(j); %distance matrix for collocation points
end
Y=Y.*(sign(Y)+1)/2;
% parameters
L = 4*c;           % length of fusealge
D = 0.12*c;        % mean fuselage diameter
D_min = 0.10*c;    % fuselage diameter at tail
D_max = 0.14*c;    % fuselage diameter at wing
AR_hor = 5;        % hor. tail AR
AR_fin = 1.2;      % tail fin AR
b_f    = 0.05;     % flange width (fraction of chord)
h_prof = c*0.125;  % airfoil/profile height at 0.3c
S_hor = c*A*0.35/L;
S_fin = c*A*0.045/L;
c_hor = sqrt(S_hor/AR_hor);
b_hor = sqrt(S_hor*AR_hor);
h_prof_hor = c_hor*0.08;  % airfoil/profile height at 0.3c
c_fin = sqrt(S_fin/AR_fin);
b_fin = sqrt(S_fin*AR_fin);
h_prof_fin = c_fin*0.08;  % airfoil/profile height at 0.3c

%% initial structure mass guess based on A. Noth's model
m_struct=1/g*0.44*b^3.1*AR^-0.25;
A_tot=A+b_fin*c_fin+b_hor*c_hor+0.2*L*D;
wing_m_distr =ones(size(y,2),1)*m_struct*A/A_tot/N;
fin_m = m_struct*b_fin*c_fin/A_tot;
hor_m = m_struct*b_hor*c_hor/A_tot;
fus_m = m_struct*0.2*L*D/A_tot;

% also guess the propulsion group mass - disregarding itself
c_angle = (20-b*0.18)/180*pi;
v_cn = sqrt((m_struct+m_distr+m_pld)*g/(0.5*rho*A*0.8)); 
P_needed=v_cn*(m_struct+m_distr+m_pld)*g*sin(c_angle)...
    +sqrt((m_struct+m_distr+m_pld*g)^3/(0.5*rho*A))/20;
m_propulsion = P_needed/n_propulsion*k_prop;% all motors!

% and now the total
m_tot = m_struct+m_distr+m_pld+m_propulsion;

%% main while
rel_delta_m=1;
while rel_delta_m>0.001
    %% polar
    c_L = -1.5:0.01:1.8;
    c_D_ind = 1.3*c_L.^2/(pi*AR);       % induced drag
    %c_D_inter = 0.07*(c*0.12/c_hor)^2; % interference drag

    for i=1:1:size(ReList,2)
        c_l_distr{i}=0.5*(1+4/pi*sqrt(1-(2*y'/b).^2))*c_L; % c_l distribution
        c_d_distr{i}=interp1(wingPolar{i}(:,2),wingPolar{i}(:,3),c_l_distr{i});
        %alpha_distr{i}=pi/180*interp1(wingPolar{i}(:,2),...
            %wingPolar{i}(:,1),c_l_distr{i});
        c_D_fus=0.074*(ReList(i)/c*L)^(-0.2)...
            *(L*(D/2)^2*pi/A); % fuselage drag (flat plate)
        c_D_hor = 0.074*(ReList(i)/c*c_hor)^(-0.2)*S_hor/A; % hor. tail drag
        c_D_fin = 0.074*(ReList(i)/c*c_fin)^(-0.2)*S_fin/A; % hor. tail drag
        c_D{i}=mean(c_d_distr{i})+c_D_ind+c_D_hor+c_D_fin;  %TODO: fuselage drag missing?
        [gr_max(i),j_gr_max(i)]=max(c_L./c_D{i});
        c_L_gr_max(i)=c_L(j_gr_max(i));
        c_D_gr_max(i)=c_D{i}(j_gr_max(i));
        [cn_max(i),j_cn_max(i)]=max(c_L.^(3/2)./c_D{i});
        c_L_cn_max(i)=c_L(j_cn_max(i));
        c_D_cn_max(i)=c_D{i}(j_cn_max(i));
        [c_L_max(i),j_c_L_max]=max(c_L./c_D{i}.*c_D{i});
        [c_L_min(i),j_c_L_min]=min(c_L./c_D{i}.*c_D{i});
    end

    %% load cases / aerodynamic distributions

    % lift increase approximation:
    % according to McCormick
    dcl_dalpha=2*pi*AR/(AR+(2*(AR+4)/(AR+2)));

    % maximum positive load factor
    n_max_i = 2.1+10900/(m_tot+4536);

    % velocities:

    % manoeuvre speed
    delta_v_m=1;
    v_m_old=0;
    c_L_max_m=mean(c_L_max);
    v_m = sqrt(n_max_i*m_tot*g/(0.5*rho*A*c_L_max_m)); 
    while abs(delta_v_m)>0.01
        Re=min([rho*c*v_m/mu,ReList(length(ReList))]);
        v_m_old = v_m;
        c_L_max_m = interp1(log(ReList),c_L_max,log(Re));
        v_m = sqrt(n_max_i*m_tot*g/(0.5*rho*A*c_L_max_m)); % logarithmic interp.
        delta_v_m=(v_m-v_m_old)/v_m;
    end

    % most negative load factor
    c_L_min_m=interp1(log(ReList),c_L_min,log(Re));
    n_min_i=0.5*rho*A*v_m^2*c_L_min_m/(m_tot*g);

    v_d = 1.5*v_m;  % dive speed

    % loads at manoeuvre speed:
    mu_g=2*m_tot/(rho*A*c*dcl_dalpha);
    Kg=0.88*mu_g/(5.3+mu_g);

    delta_L=min(c_L_max_m,atan(15.25/v_m)*dcl_dalpha)*0.5*rho*v_m^2*A;
    delta_L2=max(c_L_min_m,-atan(15.25/v_m)*dcl_dalpha)*0.5*rho*v_m^2*A;
    n_max_m=max(n_max_i,1+Kg*delta_L/(m_tot*g));
    n_min_m=min(n_min_i,1+Kg*delta_L2/(m_tot*g));

    % loads at dive speed:
    delta_L=min(c_L_max_m,atan(7.6/v_d)*dcl_dalpha)*0.5*rho*v_d^2*A;
    delta_L2=max(c_L_min_m,-atan(7.6/v_d)*dcl_dalpha)*0.5*rho*v_d^2*A;
    n_max_d=max(n_max_i,1+Kg*delta_L/(m_tot*g));
    n_min_d=min(n_min_i,1+Kg*delta_L2/(m_tot*g));
    
    % calculate the propulsion group mass
    delta_v_cn=1;
    v_cn_old=0;
    c_L_cn=mean(c_L_cn_max);
    v_cn = sqrt(m_tot*g/(0.5*rho*A*c_L_cn)); 
    while abs(delta_v_cn)>0.01
        Re_cn=min([rho*c*v_cn/mu,ReList(length(ReList))]);
        v_cn_old = v_cn;
        c_L_cn = interp1(log(ReList),c_L_cn_max,log(Re_cn));
        v_cn = sqrt(m_tot*g/(0.5*rho*A*c_L_cn)); % logarithmic interp.
        delta_v_cn=(v_cn-v_cn_old)/v_cn;
    end
    c_angle = (20-b*0.18)/180*pi;
    P_needed=v_cn*m_tot*g*sin(c_angle)...
        +sqrt((m_tot*g)^3/(0.5*rho*A))/interp1(log(ReList),cn_max,log(Re));
    m_propulsion = P_needed/n_propulsion*k_prop;% all motors!
    
    % At v_d and v_m
    %---------------
    % v_d
    Re_d=min([rho*c*v_d/mu,ReList(length(ReList))]);
    %
    c_L_d_max=2*n_max_d*m_tot*g/(rho*A*v_d^2);
    c_l_d_max=0.5*(1+4/pi*sqrt(1-(2*y'/b).^2))*c_L_d_max;
    ReU=0;
    ReL=0;
    for i=1:size(ReList,2)-1
        if ReList(i+1)>Re_d && ReList(i)<=Re_d
            ReU=ReList(i+1);
            ReL=ReList(i);
            break;
        end
    end
    if Re_d <= ReList(1)
        ReU=ReList(2);
        ReL=ReList(1);
    end
    if Re_d >= ReList(size(ReList,2))
        ReU=ReList(size(ReList,2));
        ReL=ReList(size(ReList,2)-1);
    end

    c_d_d_max_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,3),c_l_d_max);
    c_d_d_max_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,3),c_l_d_max);
    c_m_d_max_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,4),c_l_d_max);
    c_m_d_max_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,4),c_l_d_max);
    r=(log(Re_d)-log(ReL))/(log(ReU)-log(ReL));
    c_d_d_max=c_d_d_max_L+r*(c_d_d_max_U-c_d_d_max_L);
    alpha_d_max_U=pi/180*interp1(wingPolar{i+1}(:,2),...
        wingPolar{i+1}(:,1),c_l_d_max);
    alpha_d_max_L=pi/180*interp1(wingPolar{i}(:,2),...
        wingPolar{i}(:,1),c_l_d_max);
    alpha_d_max=alpha_d_max_L+r*(alpha_d_max_U-alpha_d_max_L);
    c_m_d_max=c_m_d_max_L+r*(c_m_d_max_U-c_m_d_max_L);

    c_L_d_min=2*n_min_d*m_tot*g/(rho*A*v_d^2);
    c_l_d_min=0.5*(1+4/pi*sqrt(1-(2*y'/b).^2))*c_L_d_min;

    c_d_d_min_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,3),c_l_d_min);
    c_d_d_min_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,3),c_l_d_min);
    c_m_d_min_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,4),c_l_d_min);
    c_m_d_min_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,4),c_l_d_min);
    c_d_d_min=c_d_d_min_L+r*(c_d_d_min_U-c_d_d_min_L);
    alpha_d_min_U=pi/180*interp1(wingPolar{i+1}(:,2),...
        wingPolar{i+1}(:,1),c_l_d_min);
    alpha_d_min_L=pi/180*interp1(wingPolar{i}(:,2),...
        wingPolar{i}(:,1),c_l_d_min);
    alpha_d_min=alpha_d_min_L+r*(alpha_d_min_U-alpha_d_min_L);
    c_m_d_min=c_m_d_min_L+r*(c_m_d_min_U-c_m_d_min_L);

    % v_m
    Re_m=rho*v_m*c/mu;
    ReU=0;
    ReL=0;
    for i=1:size(ReList,2)-1
        if ReList(i+1)>Re_m && ReList(i)<=Re_m
            ReU=ReList(i+1);
            ReL=ReList(i);
            break;
        end
    end
    if Re_m <= ReList(1)
        ReU=ReList(2);
        ReL=ReList(1);
    end
    if Re_m >= ReList(size(ReList,2))
        ReU=ReList(size(ReList,2));
        ReL=ReList(size(ReList,2)-1);
    end
    c_L_m_max=min(c_L_max(i),c_L_max(i+1));
    c_l_m_max=0.5*(1+4/pi*sqrt(1-(2*y'/b).^2))*c_L_m_max;


    c_d_m_max_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,3),c_l_m_max);
    c_d_m_max_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,3),c_l_m_max);
    c_m_m_max_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,4),c_l_m_max);
    c_m_m_max_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,4),c_l_m_max);
    r=(log(Re_m)-log(ReL))/(log(ReU)-log(ReL));
    c_d_m_max=c_d_m_max_L+r*(c_d_m_max_U-c_d_m_max_L);
    alpha_m_max_U=pi/180*interp1(wingPolar{i+1}(:,2),...
        wingPolar{i+1}(:,1),c_l_m_max);
    alpha_m_max_L=pi/180*interp1(wingPolar{i}(:,2),...
        wingPolar{i}(:,1),c_l_m_max);
    alpha_m_max=alpha_m_max_L+r*(alpha_m_max_U-alpha_m_max_L);
    c_m_m_max=c_m_m_max_L+r*(c_m_m_max_U-c_m_m_max_L);

    c_L_m_min=max(c_L_min(i),c_L_min(i+1));
    c_l_m_min=0.5*(1+4/pi*sqrt(1-(2*y'/b).^2))*c_L_m_min;

    c_d_m_min_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,3),c_l_m_min);
    c_d_m_min_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,3),c_l_m_min);
    c_m_m_min_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,4),c_l_m_min);
    c_m_m_min_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,4),c_l_m_min);
    c_d_m_min=c_d_m_min_L+r*(c_d_m_min_U-c_d_m_min_L);
    alpha_m_min_U=pi/180*interp1(wingPolar{i+1}(:,2),...
        wingPolar{i+1}(:,1),c_l_m_min);
    alpha_m_min_L=pi/180*interp1(wingPolar{i}(:,2),...
        wingPolar{i}(:,1),c_l_m_min);
    alpha_m_min=alpha_m_min_L+r*(alpha_m_min_U-alpha_m_min_L);
    c_m_m_min=c_m_m_min_L+r*(c_m_m_min_U-c_m_m_min_L);

    % lift distribution for roll at v_m
    c_L_m_r=2*m_tot*g/(rho*A*v_m^2);
    c_d_m_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,3),c_L_m_r);
    c_d_m_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,3),c_L_m_r);
    c_m_m_U=interp1(wingPolar{i+1}(:,2),wingPolar{i+1}(:,4),c_L_m_r);
    c_m_m_L=interp1(wingPolar{i}(:,2),wingPolar{i}(:,4),c_L_m_r);
    c_d_m=c_d_m_L+r*(c_d_m_U-c_d_m_L);
    c_m_m=c_m_m_L+r*(c_m_m_U-c_m_m_L);
    c_l_m_rp=[ones(1,floor(size(y,2)/2))*c_L_m_r,...
        wingPolarp20(2)*ones(1,size(y,2)-floor(size(y,2)/2))]';
    c_d_m_rp=[ones(1,floor(size(y,2)/2))*c_d_m,...
        wingPolarp20(3)*ones(1,size(y,2)-floor(size(y,2)/2))]';
    c_m_m_rp=[ones(1,floor(size(y,2)/2))*c_m_m,...
        wingPolarp20(4)*ones(1,size(y,2)-floor(size(y,2)/2))]';
    c_l_m_rm=[ones(1,floor(size(y,2)/2))*c_L_m_r,...
        wingPolarm20(2)*ones(1,size(y,2)-floor(size(y,2)/2))]';
    c_d_m_rm=[ones(1,floor(size(y,2)/2))*c_d_m,...
        wingPolarm20(3)*ones(1,size(y,2)-floor(size(y,2)/2))]';
    c_m_m_rm=[ones(1,floor(size(y,2)/2))*c_m_m,...
        wingPolarm20(4)*ones(1,size(y,2)-floor(size(y,2)/2))]';

    %% inertia tensor recalculation
    Ixx=sum(Y.*Y*(m_distr/b*delta+wing_m_distr));
    for k=1:floor(n/2)
        Ixx=Ixx+2*y(k*round(N/(n+1)))*y(k*round(N/(n+1)))*m_propulsion;
    end
    Ixx=Ixx+(b_hor/4)^2*hor_m+(b_fin/2)^2*fin_m;
    Izz=Ixx+L^2*(hor_m+fin_m)+(L/2)^2*fus_m;
    Iyy=L^2*(hor_m+fin_m)+(L/2)^2*fus_m;

    %% structure resizing

    %minimum thicknesses:
    t_cfrp_min=0.055/(rho_mat(1)-0.65*1140);
    t_gfrp_min=0.025/(rho_mat(3)-0.65*1140);
    t_core_min=0.5e-3; % fixed to minimum half a millimeter

    %wing flanges 
    %------------

    % roll:
    M_r=sum(y.*(0.5*rho*v_m^2*c*delta*c_l_m_rp)'-y.*(0.5*rho*v_m^2*c*delta*c_l_m_rm)'); %roll moment
    phi_dd=M_r/Ixx;
    n_loc_p=1+y*phi_dd/g; % local roll induced acceleration exxaggeration - positive
    n_loc_m=1-y*phi_dd/g; % local roll induced acceleration exxaggeration - negative
    fz_rp=-0.5*rho*v_m^2*c*delta*c_l_m_rp...
         +n_loc_p'.*(g*(m_distr/b*delta+wing_m_distr));
    for k=1:floor(n/2)
        fz_rp(k*round(N/(n+1)))=fz_rp(k*round(N/(n+1)))...
            +m_propulsion/n*n_loc_p(k*round(N/(n+1)))*g;
    end
    fz_rm=-0.5*rho*v_m^2*c*delta*c_l_m_rm...
         +n_loc_m'.*(g*(m_distr/b*delta+wing_m_distr));
    for k=1:floor(n/2)
        fz_rp(k*round(N/(n+1)))=fz_rp(k*round(N/(n+1)))...
            +m_propulsion/n*n_loc_m(k*round(N/(n+1)))*g;
    end
    fx_rp=0.5*rho*v_m^2*c*delta*c_d_m_rp;
    fx_rm=0.5*rho*v_m^2*c*delta*c_d_m_rm;
    Mbx_rp=Y*fz_rp;
    Mbz_rp=Y*fx_rp;
    Mbx_rm=Y*fz_rm;
    Mbz_rm=Y*fx_rm;

    % maximum load factors -> at v_d
    fz_d_min=0.5*rho*v_d^2*c*delta*...
         (-cos(alpha_m_min).*c_l_d_min-sin(alpha_d_min).*c_d_d_min)...
         +n_min_d*g*(m_distr/b*delta+wing_m_distr);
    for k=1:floor(n/2)
        fz_d_min(k*round(N/(n+1)))=fz_d_min(k*round(N/(n+1)))...
            +m_propulsion/n*n_min_d*g;
    end
    fx_d_min=0.5*rho*v_d^2*c*delta*...
         (sin(alpha_d_min).*c_l_d_min-cos(alpha_d_min).*c_d_d_min);
    fz_d_max=0.5*rho*v_d^2*c*delta*...
         (-cos(alpha_m_max).*c_l_d_max-sin(alpha_d_max).*c_d_d_max)...
         +n_max_d*g*(m_distr/b*delta+wing_m_distr);
    for k=1:floor(n/2)
        fz_d_max(k*round(N/(n+1)))=fz_d_max(k*round(N/(n+1)))...
            +m_propulsion/n*n_max_d*g;
    end
    fx_d_max=0.5*rho*v_d^2*c*delta*...
         (sin(alpha_d_max).*c_l_d_max-cos(alpha_d_max).*c_d_d_max);
    Mbx_d_min=Y*fz_d_min;
    Mbz_d_min=Y*fx_d_min;
    Mbx_d_max=Y*fz_d_max;
    Mbz_d_max=Y*fx_d_max;

    sigma_max_times_t_f=max(max([[abs(Mbz_d_max)*3/((c*b_f)^2)...
        +abs(Mbx_d_max)/(2*h_prof*c*b_f)],...
        [abs(Mbz_d_min)*3/((c*b_f)^2)...
        +abs(Mbx_d_min)/(2*h_prof*c*b_f)],...
        [abs(Mbz_rp)*3/((c*b_f)^2)...
        +abs(Mbx_rp)/(2*h_prof*c*b_f)],...
        [abs(Mbz_rm)*3/((c*b_f)^2)...
        +abs(Mbx_rm)/(2*h_prof*c*b_f)]]));

    % against stress
    t_f_wing_s = sigma_max_times_t_f/zul_c(2,1)*1.5;

    % against buckling (refer to Hertel)
    t_f_wing_b = (sigma_max_times_t_f/(0.4*Q_lam_c(1,1,1))*(c*b_f/2)^2*1.5)^(1/3);
    
    % bending not more than 10°
    M_integral=max([abs(sum(Mbx_d_max*delta)),abs(sum(Mbx_d_min*delta)),...
        abs(sum(Mbx_rp*delta)),abs(sum(Mbx_rm*delta))]);
    t_f_wing_f = M_integral/(tan(10/180*pi)*Q_lam_c(1,1,1)*c*b_f*h_prof^2/2);
    
    t_f_wing=max([t_f_wing_s,t_f_wing_b,t_f_wing_f]);

    %wing spar
    %---------

    % critical shear:
    Q=max(max([abs(cumsum(flipdim(fz_d_max,1)));...
        abs(cumsum(flipdim(fz_d_min,1)));...
        abs(cumsum(flipdim(fz_rm,1)));...
        abs(cumsum(flipdim(fz_rp,1)))])); % max. shear force
    s_t=T45*[0;0;Q/(h_prof)];
    t_sp_wing_s=abs(s_t(1))/(zul_c(1,1)/1.5);

    % buckling
    t_sp_wing_b = (2*Q*h_prof*1.5*rho_mat(5)/(3*4.8*rho_mat(1)*Q_lam_c(1,1,1)))^(1/3);
    t_spsw_wing = Q*h_prof*1.5/(3*4.8*Q_lam_c(1,1,1)*t_sp_wing_b^2)-t_sp_wing_b;
    t_spsw_wing = max([t_core_min,t_spsw_wing]);

    t_sp_wing=max([t_sp_wing_s,t_sp_wing_b,t_cfrp_min*2]);

    % simplified: no sandwich
    t_sp_wing_b_simple = max([(abs(s_t(1))*h_prof^2/(4.8*Q_lam_c(1,1,1)))^(1/3),t_cfrp_min,t_sp_wing_s]);
    if(t_sp_wing_b_simple*rho_mat(1)<t_sp_wing*rho_mat(1)+t_spsw_wing*rho_mat(5))
        t_sp_wing=t_sp_wing_b_simple;
        t_spsw_wing=0;
    end

    % wing shell: torsion
    % -------------------

    %torsion:
    T_m_max=c_m_m_max*c*delta*c*0.5*v_m^2*rho+0.1*c*c_l_m_max*delta*c*0.5*v_m^2*rho;
    T_m_min=c_m_m_min*c*delta*c*0.5*v_m^2*rho+0.1*c*c_l_m_min*delta*c*0.5*v_m^2*rho;
    T_m_rp=c_m_m_rp*c*delta*c*0.5*v_m^2*rho+0.1*c*c_l_m_rp*delta*c*0.5*v_m^2*rho;
    T_m_rm=c_m_m_rm*c*delta*c*0.5*v_m^2*rho+0.1*c*c_l_m_rm*delta*c*0.5*v_m^2*rho;
    T=max(max([abs(cumsum(flipdim(T_m_max,1)));...
        abs(cumsum(flipdim(T_m_min,1)));...
        abs(cumsum(flipdim(T_m_rp,1)));...
        abs(cumsum(flipdim(T_m_rm,1)))]));
    %stress:
    s_t=T45*[0;0;T/(A0*c^2)];
    t_s_wing_s_crfp=abs(s_t(1))/(zul_c(1,1)/1.5);
    t_s_wing_s_grfp=abs(s_t(1))/(zul_c(3,1)/1.5);
    %buckling:
    t_s_wing_b_crfp = (2*abs(T/(A0*c^2)*0.65*c)*0.65*c*1.5*rho_mat(5)/(3*4.8*rho_mat(1)*Q_lam_c(1,1,1)))^(1/3);
    t_s_wing_b_grfp = (2*abs(T/(A0*c^2)*0.65*c)*0.65*c*1.5*rho_mat(5)/(3*4.8*rho_mat(3)*Q_lam_c(1,1,3)))^(1/3);
    %twist (max. 3° overall)
    theta_t_crfp=max([max(abs(cumsum(flipdim(delta*T_m_max*s0/(4*(A0*c^2)^2*Q_mat_crfp(3,3)),2))));...
        max(abs(cumsum(flipdim(delta*T_m_min*s0/(4*(A0*c^2)^2*Q_mat_crfp(3,3)),2))));...
        max(abs(cumsum(flipdim(delta*T_m_rp*s0/(4*(A0*c^2)^2*Q_mat_crfp(3,3)),2))));...
        max(abs(cumsum(flipdim(delta*T_m_rm*s0/(4*(A0*c^2)^2*Q_mat_crfp(3,3)),2))))]);
    theta_t_grfp=max([max(abs(cumsum(flipdim(delta*T_m_max*s0/(4*(A0*c^2)^2*Q_mat_grfp(3,3)),2))));...
        max(abs(cumsum(flipdim(delta*T_m_min*s0/(4*(A0*c^2)^2*Q_mat_grfp(3,3)),2))));...
        max(abs(cumsum(flipdim(delta*T_m_rp*s0/(4*(A0*c^2)^2*Q_mat_grfp(3,3)),2))));...
        max(abs(cumsum(flipdim(delta*T_m_rm*s0/(4*(A0*c^2)^2*Q_mat_grfp(3,3)),2))))]);
    t_s_wing_t_crfp=theta_t_crfp/(3/180*pi);
    t_s_wing_t_grfp=theta_t_grfp/(3/180*pi);
    t_s_wing_crfp=max([t_s_wing_t_crfp;t_s_wing_b_crfp;t_s_wing_s_crfp]); 
    t_s_wing_grfp=max([t_s_wing_t_grfp;t_s_wing_b_grfp;t_s_wing_s_grfp]);
    t_sw_wing_crfp = abs(T/(A0*c^2)*0.65*c)*0.65*c*1.5/(3*4.8*Q_lam_c(1,1,1)*t_s_wing_crfp^2)-t_s_wing_crfp;
    t_sw_wing_grfp = abs(T/(A0*c^2)*0.65*c)*0.65*c*1.5/(3*4.8*Q_lam_c(1,1,3)*t_s_wing_grfp^2)-t_s_wing_grfp;
    % check if minimum/maximum thickness is respected
    t_sw_wing_crfp=min([A0*c^2/(s0*c);max([t_sw_wing_crfp, t_core_min])]);
    t_sw_wing_grfp=min([A0*c^2/(s0*c);max([t_sw_wing_grfp, t_core_min])]);
    t_s_wing_crfp=max([t_s_wing_crfp;t_s_wing_b_crfp/2+t_cfrp_min]); % resp. min theckness
    t_s_wing_grfp=max([t_s_wing_grfp;t_s_wing_b_grfp/2+t_gfrp_min]); % resp. min theckness
    if(t_s_wing_crfp*rho_mat(1)+t_sw_wing_crfp*rho_mat(5)<=...
            t_s_wing_grfp*rho_mat(3)+t_sw_wing_grfp*rho_mat(5)) % then select carbon
        t_s_wing=t_s_wing_crfp;
        t_sw_wing=t_sw_wing_crfp;
        mat_s_wing=1;
    else
        t_s_wing=t_s_wing_grfp;
        t_sw_wing=t_sw_wing_grfp;
        mat_s_wing=2;
    end
    % 

    % Hor. Stabilizer - pull up
    % -------------------------
    Fz=tailPolarp30(2)*b_hor*c_hor*0.5*rho*v_m^2;
    theta_dd=Fz*L/Iyy;
    n_hor=1-theta_dd*L/g;
    Mbx_hor=1/2/2*b_hor*hor_m*n_hor*g-1/2/2*tailPolarp30(2)*b_hor*c_hor*0.5*rho*v_m^2;
    Q_hor=hor_m*n_hor*g-tailPolarp30(2)*b_hor*c_hor*0.5*rho*v_m^2/2;
    T_hor=b_hor/2*tailPolarp30(4)*b_hor*c_hor*c_hor*0.5*rho*v_m^2/2+0.1*c_hor*tailPolarp30(2)*b_hor*c_hor*0.5*rho*v_m^2/2;
    % flanges against stress
    t_f_hor_s = abs(Mbx_hor)/(2*h_prof_hor*c_hor*b_f)/zul_c(2,1)*1.5;
    % flanges against buckling (refer to Hertel)
    t_f_hor_b = (abs(Mbx_hor)/(2*h_prof_hor*c_hor*b_f)/(0.4*Q_lam_c(1,1,1))*(c_hor*b_f/2)^2*1.5)^(1/3);
    t_f_hor=max([t_f_hor_s,t_f_hor_b]);
    % spar against stress
    s_t=T45*[0;0;Q_hor/(h_prof_hor)];
    t_sp_hor_s=abs(s_t(1))/(zul_c(1,1)/1.5);
    % spar against buckling
    t_sp_hor_b = (2*abs(Q_hor)*h_prof_hor*1.5*rho_mat(5)/(3*4.8*rho_mat(1)*Q_lam_c(1,1,1)))^(1/3);
    t_spsw_hor = abs(Q_hor)*h_prof_hor*1.5/(3*4.8*Q_lam_c(1,1,1)*t_sp_wing_b^2)-t_sp_wing_b;
    t_spsw_hor = max([t_core_min,t_spsw_hor]);
    t_sp_hor=max([t_sp_hor_s,t_sp_hor_b,2*t_cfrp_min]);
    % spar simplified: no sandwich
    t_sp_hor_b_simple = max([(abs(s_t(1))*h_prof_hor^2/(4.8*Q_lam_c(1,1,1)))^(1/3),t_cfrp_min,t_sp_hor_s]);
    if(t_sp_hor_b_simple*rho_mat(1)<t_sp_hor*rho_mat(1)+t_spsw_hor*rho_mat(5))
        t_sp_hor=t_sp_hor_b_simple;
        t_spsw_hor=0;
    end
    % shell stress:
    s_t=T45*[0;0;T_hor/(A0*c_hor^2)];
    t_s_hor_s_crfp=abs(s_t(1))/(zul_c(1,1)/1.5);
    t_s_hor_s_grfp=abs(s_t(1))/(zul_c(3,1)/1.5);
    % shell buckling: (shell not simplified: inner shell thickness goes to zero anyways, core material can use full space)
    t_s_hor_b_crfp = (2*abs(T_hor/(A0_tail*c_hor^2)*0.65*c_hor)*0.65*c_hor*1.5*rho_mat(5)/(3*4.8*rho_mat(1)*Q_lam_c(1,1,1)))^(1/3);
    t_s_hor_b_grfp = (2*abs(T_hor/(A0_tail*c_hor^2)*0.65*c_hor)*0.65*c_hor*1.5*rho_mat(5)/(3*4.8*rho_mat(3)*Q_lam_c(1,1,3)))^(1/3);
    t_s_hor_crfp=max([t_s_wing_b_crfp;t_s_hor_s_crfp]); 
    t_s_hor_grfp=max([t_s_wing_b_grfp;t_s_hor_s_grfp]);
    t_sw_hor_crfp = abs(T/(A0_tail*c_hor^2)*0.65*c_hor)*0.65*c_hor*1.5/(3*4.8*Q_lam_c(1,1,1)*t_s_hor_crfp^2)-t_s_hor_crfp;
    t_sw_hor_grfp = abs(T/(A0_tail*c_hor^2)*0.65*c_hor)*0.65*c_hor*1.5/(3*4.8*Q_lam_c(1,1,3)*t_s_hor_grfp^2)-t_s_hor_grfp;
    t_sw_hor_crfp=min([A0_tail*c_hor^2/(s0_tail*c_hor);max([t_sw_hor_crfp, t_core_min])]);
    t_sw_hor_grfp=min([A0_tail*c_hor^2/(s0_tail*c_hor);max([t_sw_hor_grfp, t_core_min])]);
    t_s_hor_crfp=max([t_s_hor_crfp;t_s_hor_b_crfp/2+t_cfrp_min]); % resp. min theckness
    t_s_hor_grfp=max([t_s_hor_grfp;t_s_hor_b_grfp/2+t_gfrp_min]); % resp. min theckness
    if(t_s_hor_crfp*rho_mat(1)+t_sw_hor_crfp*rho_mat(5)<=...
            t_s_hor_grfp*rho_mat(3)+t_sw_hor_grfp*rho_mat(5)) % then select carbon
        t_s_hor=t_s_hor_crfp;
        t_sw_hor=t_sw_hor_crfp;
        mat_s_hor=1;
    else
        t_s_hor=t_s_hor_grfp;
        t_sw_hor=t_sw_hor_grfp;
        mat_s_hor=2;
    end

    % Fin - yaw
    % -------------------------
    Fy=tailPolarp30(2)*b_fin*c_fin*0.5*rho*v_m^2;
    psi_dd=Fy*L/Izz;
    n_fin=psi_dd*L/g;
    Mbx_fin=-1/2*b_fin*fin_m*n_fin*g+1/2*tailPolarp30(2)*b_fin*c_fin*0.5*rho*v_m^2;
    Q_fin=-fin_m*n_fin*g+tailPolarp30(2)*b_fin*c_fin*0.5*rho*v_m^2;
    T_fin=b_fin*tailPolarp30(4)*b_fin*c_fin*c_fin*0.5*rho*v_m^2+0.1*c_fin*tailPolarp30(2)*b_fin*c_fin*0.5*rho*v_m^2;
    % flanges against stress
    t_f_fin_s = abs(Mbx_fin)/(2*h_prof_fin*c_fin*b_f)/zul_c(2,1)*1.5;
    % flanges against buckling (refer to Hertel)
    t_f_fin_b = (abs(Mbx_fin)/(2*h_prof_fin*c_fin*b_f)/(0.4*Q_lam_c(1,1,1))*(c_fin*b_f/2)^2*1.5)^(1/3);
    t_f_fin=max([t_f_fin_s,t_f_fin_b]);
    % spar against stress
    s_t=T45*[0;0;Q_fin/(h_prof_fin)];
    t_sp_fin_s=abs(s_t(1))/(zul_c(1,1)/1.5);
    % spar against buckling
    t_sp_fin_b = (2*abs(Q_fin)*h_prof_fin*1.5*rho_mat(5)/(3*4.8*rho_mat(1)*Q_lam_c(1,1,1)))^(1/3);
    t_spsw_fin = abs(Q_fin)*h_prof_fin*1.5/(3*4.8*Q_lam_c(1,1,1)*t_sp_wing_b^2)-t_sp_wing_b;
    t_spsw_fin = max([t_core_min,t_spsw_fin]);
    t_sp_fin=max([t_sp_fin_s,t_sp_fin_b,2*t_cfrp_min]);
    % spar simplified: no sandwich
    t_sp_fin_b_simple = max([(abs(s_t(1))*h_prof_fin^2/(4.8*Q_lam_c(1,1,1)))^(1/3),t_cfrp_min,t_sp_fin_s]);
    if(t_sp_fin_b_simple*rho_mat(1)<t_sp_fin*rho_mat(1)+t_spsw_fin*rho_mat(5))
        t_sp_fin=t_sp_fin_b_simple;
        t_spsw_fin=0;
    end
    % shell stress:
    s_t=T45*[0;0;T_fin/(A0*c_fin^2)];
    t_s_fin_s_crfp=abs(s_t(1))/(zul_c(1,1)/1.5);
    t_s_fin_s_grfp=abs(s_t(1))/(zul_c(3,1)/1.5);
    % shell buckling: (shell not simplified: inner shell thickness goes to zero anyways, core material can use full space)
    t_s_fin_b_crfp = (2*abs(T_fin/(A0_tail*c_fin^2)*0.65*c_fin)*0.65*c_fin*1.5*rho_mat(5)/(3*4.8*rho_mat(1)*Q_lam_c(1,1,1)))^(1/3);
    t_s_fin_b_grfp = (2*abs(T_fin/(A0_tail*c_fin^2)*0.65*c_fin)*0.65*c_fin*1.5*rho_mat(5)/(3*4.8*rho_mat(3)*Q_lam_c(1,1,3)))^(1/3);
    t_s_fin_crfp=max([t_s_wing_b_crfp;t_s_fin_s_crfp]); 
    t_s_fin_grfp=max([t_s_wing_b_grfp;t_s_fin_s_grfp]);
    t_sw_fin_crfp = abs(T/(A0_tail*c_fin^2)*0.65*c_fin)*0.65*c_fin*1.5/(3*4.8*Q_lam_c(1,1,1)*t_s_fin_crfp^2)-t_s_fin_crfp;
    t_sw_fin_grfp = abs(T/(A0_tail*c_fin^2)*0.65*c_fin)*0.65*c_fin*1.5/(3*4.8*Q_lam_c(1,1,3)*t_s_fin_grfp^2)-t_s_fin_grfp;
    t_sw_fin_crfp=min([A0_tail*c_fin^2/(s0_tail*c_fin);max([t_sw_fin_crfp, t_core_min])]);
    t_sw_fin_grfp=min([A0_tail*c_fin^2/(s0_tail*c_fin);max([t_sw_fin_grfp, t_core_min])]);
    t_s_fin_crfp=max([t_s_fin_crfp;t_s_fin_b_crfp/2+t_cfrp_min]); % resp. min theckness
    t_s_fin_grfp=max([t_s_fin_grfp;t_s_fin_b_grfp/2+t_gfrp_min]); % resp. min theckness
    if(t_s_fin_crfp*rho_mat(1)+t_sw_fin_crfp*rho_mat(5)<=...
            t_s_fin_grfp*rho_mat(3)+t_sw_fin_grfp*rho_mat(5)) % then select carbon
        t_s_fin=t_s_fin_crfp;
        t_sw_fin=t_sw_fin_crfp;
        mat_s_fin=1;
    else
        t_s_fin=t_s_fin_grfp;
        t_sw_fin=t_sw_fin_grfp;
        mat_s_fin=2;
    end

    % Fuselage: pull up and yaw (laminate orientation 0° helical angle)
    % -----------------------------------------------------------------
    delta_fus=L/N;
    y_fus=delta_fus/2:delta_fus:L-delta_fus/2;
    y_fus2=0:delta_fus:L;
    for j=1:size(y_fus,2)
        Y_fus(j,:)=y_fus-y_fus2(j); %distance matrix for collocation points
    end
    Y_fus=Y_fus.*(sign(Y_fus)+1)/2;
    p_fus=ones(1,N);
    Qy_fus=(Fy-psi_dd*y_fus.*(L-y_fus)/L*fus_m/N-(hor_m+fin_m)*g*L*psi_dd);
    Qz_fus=(Fz+(g-theta_dd)*y_fus.*(L-y_fus)/L*fus_m/N-(hor_m+fin_m)*(g-theta_dd*L));
    Mby_fus=(Fz*(L-y_fus)+(g-theta_dd)*y_fus.*(L-y_fus)/L*fus_m/N.*(L-y_fus)/2+(L-y_fus)*(hor_m+fin_m)*(g-theta_dd*L))+T_hor;
    Mbz_fus=(Fy*(L-y_fus)-psi_dd*y_fus.*(L-y_fus)/L*fus_m/N.*(L-y_fus)/2-(hor_m+fin_m)*g*L.*(L-y_fus)*psi_dd)+T_fin;
%     Mb_fus=sqrt((T_fin+L*Q_fin-n_fin/2*g*(fus_m)*L/2-n_fin*hor_m*g)^2 + ...
%         (Fz*L+(1-theta_dd*L/2/g)*g*(fus_m)*L/2+(1-theta_dd*L/g)*g*(fin_m+hor_m)*L)^2); %attention: second part equals steady moment (to be compensated w/ pld, e.g.
%     Mb_fus=max([Mb_fus,abs(2*T_hor)]); %L*g*m_tot sort of hackish - doesn't converge
%     Q_fus=sqrt(Q_fin^2+Q_hor^2); %long beam, thus not relevant
%     T_fus=Mbx_fin;
    Mb_fus=sqrt(Mby_fus.^2+Mbz_fus.^2);
    s_t = max(Mb_fus)/(pi*D^3/8)*D/2;

    % stress
    t_fus_s = 1.5*s_t/zul_c(1,1);

    % max. 2° deflection
    t_fus_d=sum(delta_fus*Mb_fus/(pi*D^3/8)/(2/180*pi)/(Q_lam_c(1,1,1)));

    % bending buckling, 
    % cf Designing composite structures: lay-up selection
    % P M Weaver, Department of Aerospace Engineering, University of Bristol,

    t_fus_b=sqrt(1.5*s_t/(0.3*Q_lam_c(1,1,1))*(D/2));

    t_fus=max([t_cfrp_min,t_fus_b,t_fus_d]);
    
    %% recalculate mass distribution
    damping=0.8;
    % wing mass distr
    wing_m_distr =damping*wing_m_distr + (1-damping)*mult*(1.2*ones(size(y,2),1)*delta*(...
        rho_mat(mat_s_wing)*t_s_wing*s0*c+rho_mat(5)*t_sw_wing*s0*c...
        +2*b_f*c*t_f_wing*rho_mat(1)...
        +h_prof*t_sp_wing*rho_mat(1) + h_prof*t_spsw_wing*rho_mat(5)));
    fin_m =damping*fin_m + (1-damping)*mult*(1.2*b_fin*(...
        rho_mat(mat_s_fin)*t_s_fin*s0_tail*c_fin+rho_mat(5)*t_sw_fin*s0_tail*c_fin...
        +2*b_f*c_fin*t_f_fin*rho_mat(1)...
        +h_prof_fin*t_sp_fin*rho_mat(1) + h_prof_fin*t_spsw_fin*rho_mat(5)));
    hor_m =damping*hor_m + (1-damping)*mult*(1.2*b_hor*(...
        rho_mat(mat_s_hor)*t_s_hor*s0_tail*c_hor+rho_mat(5)*t_sw_hor*s0_tail*c_hor...
        +2*b_f*c_hor*t_f_fin*rho_mat(1)...
        +h_prof_hor*t_sp_hor*rho_mat(1) + h_prof_hor*t_spsw_hor*rho_mat(5)));
    fus_m=damping*fus_m + (1-damping)*mult*(1.2*L*D*pi*t_fus*rho_mat(2));
    m_struct_new=fus_m+hor_m+fin_m+2*sum(wing_m_distr);
    m_tot_new=m_struct+m_distr+m_pld+m_propulsion;
    rel_delta_m=abs(1-m_struct_new/m_struct);
    m_tot=m_tot_new;
    %disp(m_struct);
    m_struct=m_struct_new;
end

% max. glide ratio
delta_v_gr=1;
c_L_gr=mean(c_L_gr_max);
v_gr = sqrt(m_tot*g/(0.5*rho*A*c_L_gr)); 
while abs(delta_v_gr)>0.01
    Re=min([rho*c*v_gr/mu,ReList(length(ReList))]);
    v_gr_old = v_gr;
    c_L_gr = interp1(log(ReList),c_L_gr_max,log(Re));
    v_gr = sqrt(m_tot*g/(0.5*rho*A*c_L_gr)); % logarithmic interp.
    delta_v_gr=(v_gr-v_gr_old)/v_gr;
end
c_D_gr = interp1(log(ReList),c_D_gr_max,log(Re));

% stall speed
delta_v_s=1;
c_L_s=mean(c_L_max);
v_s = sqrt(m_tot*g/(0.5*rho*A*c_L_s)); 
while abs(delta_v_s)>0.01
    Re_s=min([rho*c*v_s/mu,ReList(length(ReList))]);
    v_s_old = v_s;
    c_L_s = interp1(log(ReList),c_L_max,log(Re_s));
    v_s = sqrt(m_tot*g/(0.5*rho*A*c_L_s)); % logarithmic interp.
    delta_v_s=(v_s-v_s_old)/v_s;
end
c_L_s_min = interp1(log(ReList),c_L_min,log(Re_s));

%disp(m_struct);
thicknesses.t_f_wing  = t_f_wing;  % wing flange thickness
thicknesses.t_s_wing  = t_s_wing;  % wing shell thickness (overall)
thicknesses.t_sp_wing = t_sp_wing; % wing spar thickness (overall)
thicknesses.t_spsw_wing = t_spsw_wing; % wing spar sandwich thickness
thicknesses.mat_s_wing = mat_s_wing; % 1=f_cfrp, 3=f_gfrp
thicknesses.t_sw_wing = t_sw_wing; % wing sandwich thickness
thicknesses.t_fus     = t_fus; % fuselage shell thickness
thicknesses.t_f_hor  = t_f_hor;  % hor flange thickness
thicknesses.t_s_hor  = t_s_hor;  % hor shell thickness (overall)
thicknesses.t_sp_hor = t_sp_hor; % hor spar thickness (overall)
thicknesses.t_spsw_hor = t_spsw_hor; % hor spar sandwich thickness
thicknesses.mat_s_hor = mat_s_hor; % 1=f_cfrp, 3=f_gfrp
thicknesses.t_sw_hor = t_sw_hor; % hor sandwich thickness
thicknesses.t_f_fin  =t_f_fin;  % fin flange thickness
thicknesses.t_s_fin  = t_s_fin;  % fin shell thickness (overall)
thicknesses.t_sp_fin = t_sp_fin; % fin spar thickness (overall)
thicknesses.t_spsw_fin = t_spsw_fin; % fin spar sandwich thickness
thicknesses.mat_s_fin = mat_s_fin; % 1=f_cfrp, 3=f_gfrp
thicknesses.t_sw_fin = t_sw_fin; % fin sandwich thickness
masses.m_fuselage=fus_m;
masses.m_horstabilizer=hor_m;
masses.m_fin=fin_m;
masses.m_wing=2*sum(wing_m_distr);
masses.m_prop=m_propulsion;
velocities.v_m=v_m;
velocities.v_d=v_d;
velocities.v_gr=v_gr;
velocities.v_cn=v_cn;
velocities.v_s=v_s;
polar.ReList=ReList;
polar.c_L_cn = c_L_cn_max;
polar.c_D_cn = c_D_cn_max;
polar.c_L_gr = c_L_gr_max;
polar.c_D_gr = c_D_gr_max;
polar.c_L_max = c_L_max;
polar.c_L_min = c_L_min;
polar.c_L_array = c_L;  %The full CD/CL-polars
polar.c_D_array = c_D;


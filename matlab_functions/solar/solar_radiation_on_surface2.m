%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Moix Pierre-Olivier     %
%     January 2004         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFIED FOR THE INDEPENDENT SOLAR RADIATION CALCULATOR
%
% solar_radiation_on_surface.m
% modeling of the sun irradiance on a given surface at a given time
% The model is based on Duffie & Beckman: "Solar Engineering of Thermal Processes"
% the equations numbers is refered to this book
%
% Derived from a base script of N. Morel at start but fully revised with
% the equations of the book
%
%
% The function calculates an estimate of the instantaneous
% solar radiation during at a time ts0 on a surface
% of arbitrary slope and orientation.
%
% Parameters:
%     slope in degrees: horizontal=0, vertical=90
%     orien in degrees South=0, East=90, West=-90, North=180
%     ts0 is time in seconds
%     altitude in meters
%     albedo, from 0 to 1
%        typically albedo is 
%               fresh snow 0.8 to 0.9
%               green land with grass 0.12 to 0.25
%               sand   0.25 to 0.45
%               forest 0.05 to 0.2
%               sea and ocean   0.02 to 0.05 if sun height >30 degrees
%                               0.02 to 0.2  if sun height <10 degrees 
%
%
%
% The output of the function is an array with the
% following components (the components not calculated
% because of 'flagSolOutput' value are set to zero):
% r(1)=global radiation on the surface
% r(2)=direct radiation on the surface
% r(3)=diffuse radiation on the surface
% r(4)=global horizontal radiation
% r(5)=direct horizontal radiation
% r(6)=diffuse horizontal radiation
%
% r(7)=incidence_angle=incidence angle theta on surface (radians)
% r(8)=azimut_angle=direction of the sun compared to the direction of the south
% r(9)=zenith_angle;
%


%aa_interp=interp1(1:12,aa,day/365*12);
%    bb_interp=interp1(1:12,aa,day/365*12);
%    cc_interp=interp1(1:12,aa,day/365*12);

function r=solar_radiation_on_surface2(slope,orien,t,day,latitude,longitude, altitude,albedo)
day=day-1; % different convention: day one in input means 0 days over
ts0=(day)*86400+t;


%%%%%%%%%%%%%%%%%%%%%
% not modified
timezone=round(longitude*24/360); % [hours], West>0  
%We place ourself at Greenwich time in any case
%as basis, this will be much easier:
ts0=ts0-timezone*60*60;

%%%%%%%%%%%%%%%%
% Limits
max_incidence_angle=89.8*pi/180; %To choose with cell properties
cosmin=cos(max_incidence_angle);    % For the limit of direct irradiance over this there 
% is total reflectance or the sun is behind the cell
%fminDiff=0.2; % minimum diffuse fraction

% Changes for simplifications of the inputs
phi=latitude*pi/180;
cosphi=cos(phi);
sinphi=sin(phi);


%%%%%%%%%%%%%%%%%%%
% declination
del=23.45*(pi/180)*sin(2*pi*(284+day)/365); % [rad]
cosdel=cos(del);
sindel=sin(del);
da=2*pi*day/365;


%%%%%%%%%%%%%%%%%%%%%%%
% sunrise and sunset solar times [s] 
omega_s=acos(-sindel/cosdel*sinphi/cosphi); %hour angle at sunrise 1.6.10
tsunrise=43200*(1-omega_s/pi);
tsunset=43200*(1+omega_s/pi);


%%%%%%%%%%%%%%%%%%%%%
% et = equation of time, in seconds, 
% correction due to the speed variation of the earth around the sun
et=(0.0072*cos(da)-0.0528*cos(2*da)-0.0012*cos(3*da)...
    -0.1229*sin(da)-0.1565*sin(2*da)-0.0041*sin(3*da))*3600;


%%%%%%%%%%%%%%%%%%%%%%%%
% time difference and other constants
tdiff=3600*(timezone-longitude/15)-day*86400;%+et; do not consider et, user enters sideral time
tsSol0=ts0+tdiff; % solar times, one day range 



%%%%%%%%%%%%%%%%%%%%%%
% solar angles (Duffie & Beckman, pp 15 ff)
ha=pi*(1-tsSol0/43200); % this is the omega in the book, ha is for hour angle

sinha=sin(ha);
cosha=cos(ha);
beta=pi/180*slope;
sin_beta=sin(beta);
cos_beta=cos(beta);
gam=pi/180*orien;
singam=sin(gam);
cosgam=cos(gam);
cost=sindel*sinphi*cos_beta-sindel*cosphi*sin_beta*cosgam...
    +cosdel*cosphi*cos_beta*cosha+cosdel*sinphi*sin_beta*cosgam*cosha...
    +cosdel*sin_beta*singam*sinha;  %equation 1.6.2

if (abs(cost)>1+eps) 
    error(['*** solar: cos(t)=',num2str(cost)]); 
end

incidence_angle=acos(cost);  % [rad]

%the zenith angle: (is incidence angle for an horizontal surface)
cosz=cosphi*cosdel*cosha+sinphi*sindel;
if (abs(cosz)>1+eps)
    error(['*** solar: zenith angle=',num2str(cosz)]); 
end
zenith_angle=acos(cosz);

%h=pi/2-acos(cosz);  %sun height h (called alpha s in the book), complement of the zenith angle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solar angle for which the sun azimut is > pi/2 or < -pi/2, this happens for days longer than 12 hours
% In that case we must determine in which quadrant we are 
% equation 1.6.6g

omega_ew=acos(sindel/cosdel*cosphi/sinphi);


% %%%%%%%%%%%%%%%%%%%%%%%%%
% %With the notations of the book:
sinz=sqrt(1-cosz*cosz);
sina=cosdel*sinha/sinz; % this is eq 1.6.6b a is for gamma prime s
if (abs(ha)<=omega_ew), C1=1;
else C1=-1; end
if ((phi-del)>=0), C2=1;
else C2=-1; end
if ((ha)>=0), C3=1;
else C3=-1; end
gamma_s=C1*C2*asin(sina)+C3*(1-C1*C2)*pi/2;

%correction to put back between -pi and pi:
if (gamma_s<-pi), gamma_s=gamma_s+2*pi; end
if (gamma_s>pi), gamma_s=gamma_s-2*pi; end
% %it's OK, checked it gives the same value
% %%%%%%%%%%%%%%%%%%%

azimut_angle=gamma_s;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Radiation for a completely clear sky
% reference: T.Markus, E.Morris, "Building, Climate and Energy",
% Pitman, London (1980)

%corrections to apply for each month for radiation:
aa=[1233,1230,1214,1185,1135,1103,1088,1085,1107,1151,1192,1220,1233,1230];
bb=[142,142,144,156,180,196,205,207,201,177,160,149,142,142]*0.001;
cc=[57,58,60,71,97,121,134,136,122,92,73,57,57,57]*0.001;


if (tsSol0<tsunrise || tsSol0>tsunset), % then we are during night
    r(1:6)=[0,0,0,0,0,0];
else % day
    % linear interpolation
    i=floor((day+15)/365*12+1);
    r=((day+15-(i-1)*365/12)/30);
    aa_interp=aa(i)+(aa(i+1)-aa(i))*r;%interp1(0:12,aa,(day+15)/365*12);
	bb_interp=bb(i)+(bb(i+1)-bb(i))*r;%interp1(0:12,bb,(day+15)/365*12);
	cc_interp=cc(i)+(cc(i+1)-cc(i))*r;%interp1(0:12,cc,(day+15)/365*12);
    m=35/sqrt(1224*cosz*cosz+1);
    tau=(1-altitude/44308).^5.257;
    qndirCS=aa_interp*exp(-bb_interp*m*tau);
    if (cosz>=0),
        qhdirCS=qndirCS*cosz;
        qhdiffCS=qndirCS*cc_interp;
    else
        qhdirCS=0; qhdiffCS=0;
    end
    qhtotCS=qhdirCS+qhdiffCS;
    if (cosz<cosmin || cost<cosmin), % low sun or sun behind surface
        qdirCS=0;
    else % normal case with a direct component
        qdirCS=qhdirCS*cost/cosz;
    end
    qdiffCS=0.5*(qhdiffCS*(1+cos_beta)+qhtotCS*albedo*(1-cos_beta));
    qtotCS=qdirCS+qdiffCS;
    r(1:6)=[qtotCS,qdirCS,qdiffCS,qhtotCS,qhdirCS,qhdiffCS];
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % Radiation for a completely clear sky
% % reference: Duffie & Beckman: "Solar Engineering of Thermal Processes" 2.8
% % Hottel method 1976
% % 
% 
% %atmospheric parameters:
% a0_star=0.4237-0.00821*(6-altitude/1000)^2;
% a1_star=0.5055+0.00595*(6.5-altitude/100)^2;
% k_star=0.2711+0.01858*(2.5-altitude/1000)^2;
% 
% %corrections factors for climate types
% if abs(latitude)>65
%     climate='Subarctic summer';
%     r0=0.99;
%     r1=0.99;
%     rk=1.01;
% elseif abs(latitude)<20
%     climate='Tropical'
%     r0=0.95;
%     r1=0.98;
%     rk=1.02;
% else
%     if month>3&month<9
%         climate='Midlattitude summer';
%         r0=0.97;
%         r1=0.99;
%         rk=1.02;
%     else
%         climate='Midlattitude winter';
%         r0=1.03;
%         r1=1.01;
%         rk=1.00;
%     end
% end
% 
% a0=a0_star*r0;
% a1=a1_star*r1;
% k=k_star*rk;
% 
% % we have to put a limitation on the coz to avoid division by 0:
% if cosz<cosmin, % low sun or sun behind surface
%         limit_it=0;
%     else % normal case with a direct component
%         limit_it=1;
%     end
%     
% tau_b=(a0+a1*exp(-k/cosz))*limit_it;
% 
% %radiation received over the atmosphere in space G_on an horizontal surface
% G_on=1367*(1+0.033*cos(2*pi*day/365));
% 
% %G_cnb=G_on*tau_b;
% %The horizontal received is
% G_cb=G_on*tau_b*cosz;
% 
% 
% %The radiation on the sloped surface is not computed yet with this method



r(7)=incidence_angle;
r(8)=azimut_angle;

r(9)=zenith_angle;

%r(10)=gamma_s2;
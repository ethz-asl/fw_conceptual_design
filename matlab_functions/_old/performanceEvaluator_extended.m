% Evaluates the performance of a solar airplane
% Authors:
%  - Framework by S. Leutenegger (10/2009)
%  - Extension by P. Oettershagen (03/2014)
%       - assumes given P_Level from flight tests for airplane
%--------------------------------------------------------------------------
% Input Arguments:
% ================
% climb         1/0, allow climbing
% plane         airplane parameters  
% environment:  environmental conditions 
%   .lat        latitude [deg]
%   .T_ground   temperature on the ground [°K]
%   .day        day of the month
%   .month      month of the year
%   .clearness  clearness of the sky [-]
%   .albedo     earth surface albedo [-]
% MDflag:       Multiday-Flag (0: SingleDay-Sim only, 1: 1st day of MD-sim, 2: 2nd day of MD-sim)
% MDinputs:     Inputs provided by the performanceEvaluator_MD, only
%               necessary if MDflag==2

%the real function
function [t_excess, t_endurance, t_chargemargin, min_SoC, flightdata, MDresults] = performanceEvaluator_extended(...
    climb,plane,environment,MDflag,MDinputs)

AC_STATE=ACS.UNDEFINED;

%Setup data storage arrays
t_array=[];
h_array=[];
irr_array=[];
bat_array=[];
P_solar_array=[];
P_elec_tot_array=[];
Re_array=[];
CL_array=[];
CD_array=[];
v_array=[];
P_prop_array=[];

% set the time discretization:
Dt = 300; %[s]
DEBUG=1;
t_fullycharged=NaN;
t_eq2=NaN;
    
% some parameters:
g  = 9.806; % [m/s^2]

E_bat_max = plane.bat.e_density * plane.bat.m * 3600; % maximum battery energy state
% --------------------------------
% Now simulate the course of a day
% --------------------------------

% Calculate the "Performance" in terms of endurance/range(speed) or 
% excess time, if continuos flight is possible

% -------------------------------------------------------------------
% STEP 1: FIND THE POWER EQUILIBRIA t_eq and t_eq2, i.e. when
%         P_Solar=P_elec_tot in the morning and the evening.
% -------------------------------------------------------------------
t = 0;
h=environment.h_0;
P_solar = 0;
eq_reached_flag = 1; % needed for distinguishing two cases: 
                     % power eq. reached or not

%STEP 1a: Calc Re,Cl,Cd for flight at h_0:
[rho,mu] = findAir(environment.h_0,environment.T_ground);

%Step1b: Calc Total electric power for level flight
P_level = plane.P_prop_level*sqrt(plane.P_prop_level_rho/rho)*(1.0+environment.turbulence);
P_elec_tot = plane.P_avi+P_level;

%Step1c: Find t_sunrise, i.e. when P_solar=0.
while P_solar<=0
irr_vec=solar_radiation_on_surface2(0,0,mod(t,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,0, h,environment.albedo);
    irr = irr_vec(1);
    P_solar = environment.clearness*irr*plane.SolarModules.Surface * plane.SolarModules.eta;
    t = t+Dt;
end
t_sunrise=t;

%Step1d: Find t_eq, that means the time when P_Solar=P_elect and thus
%        the charge can start
while P_solar < P_elec_tot
    irr_vec=solar_radiation_on_surface2(0,0,mod(t,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,0, h,environment.albedo);
    irr = irr_vec(1);
    P_solar = environment.clearness*irr*plane.SolarModules.Surface * plane.SolarModules.eta;
    t = t+Dt;
    if t >= 3600*12 % /!\ 12 would not be compatible with longitude shift 
        eq_reached_flag = 0; % the equilibrium point was never reached
        break; % abort the search
    end
end
t_eq=t;
%Step1e: Find t_eq2, that means the time when P_Solar=P_electot again
%        and thus the discharge would have to start if flying at h_0
while P_solar >= P_elec_tot
    irr_vec=solar_radiation_on_surface2(0,0,mod(t,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,0, h,environment.albedo);
    irr = irr_vec(1);
    P_solar = environment.clearness*irr*plane.SolarModules.Surface * plane.SolarModules.eta;
    t = t+Dt;
end
t_eq2=t;

% -------------------------------------------------------------------
% STEP 2: Simulate the day
% -------------------------------------------------------------------
if eq_reached_flag == 1
    % -------------------------------------------------------------------
    % STEP 2a: Check if the battery can be fully charged
    % -------------------------------------------------------------------
    t=12*3600; 
    E_bat=0;
    while (E_bat < E_bat_max/2)
        irr_vec=solar_radiation_on_surface2(0,0,mod(t,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,0, h,environment.albedo);
        irr = irr_vec(1);
        P_solar = environment.clearness*irr*plane.SolarModules.Surface * plane.SolarModules.eta;
        delta_E = Dt*(P_solar-P_elec_tot)*plane.bat.eta_chrg;
        E_bat = E_bat + delta_E;
        t = t + Dt;
        if(t>(t_eq)) %ERROR HERE? t_eq2 instead of t_eq
            break;
        end
    end
    % -------------------------------------------------------------------
    % STEP 2b: Determine earliest start time when starting with full batt
    %         (step backward in time)
    % -------------------------------------------------------------------
    E_bat = max((E_bat_max-2*E_bat)/2,0);
    while E_bat < E_bat_max
        if t<0;
            irr = 0; % set zero light before midnight
        else
            irr_vec=solar_radiation_on_surface2(0,0,mod(t,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,0, h,environment.albedo);
            irr = irr_vec(1);
        end
        P_solar = environment.clearness*irr*plane.SolarModules.Surface * plane.SolarModules.eta;
        E_bat = E_bat + Dt*(P_elec_tot-P_solar)/plane.bat.eta_dchrg;
        t = t - Dt;
    end
    t_start_full = t;
    %--------------------------------------------
    % STEP 2c: MAIN SIMULATION LOOP
    %--------------------------------------------
    h=environment.h_0;
    h_before=h;

    %STARTING POINT SETTINGS (determined by caller-function)
    if(MDflag==0 ||MDflag==1)
        E_bat=0;    %Single-Day Sim, start with empty bats at t_eq
        t=t_eq;
    elseif(MDflag==2)
        E_bat = MDinputs(1); %Multi-day sim, 2nd day. Start with (partially) charged batteries.
        t = MDinputs(2);
    end

    %Begin main loop
    while E_bat >=0 && h>=environment.h_0 && t<3600*24+t_eq

        % Altitude has changed, recalc Re,CL,CD,v,Plevel,Pelectot
        if(h~=h_before) 
            h_before=h;
            [rho,mu] = findAir(h,environment.T_ground);

            %Step1b: Calc Total electric power for level flight
            P_level = plane.P_prop_level*sqrt(plane.P_prop_level_rho/rho)*(1.0+environment.turbulence);
            P_elec_tot = plane.P_avi+P_level;
        end
        irr_vec=solar_radiation_on_surface2(0,0,mod(t,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,0, h,environment.albedo);
            irr = irr_vec(1);
        P_solar = environment.clearness*irr*plane.SolarModules.Surface * plane.SolarModules.eta;

        % -------------------------------------------------------------------
        % CONTROL LAW
        % -------------------------------------------------------------------
        if(E_bat>=E_bat_max) %fully charged! ->dayflight=> level flight/climb
            E_bat=E_bat_max; 
            %Set charge completed flag once
            if(isnan(t_fullycharged)) t_fullycharged=t; end

            %Case 1:
            if(h>environment.h_0)
                if(h>environment.h_max) h=environment.h_max; end

                %CONTROL-LAW: Put P_prop=0 only when P_Solar<P_0
                if(P_solar>=plane.P_avi && h<environment.h_max)
                    P_prop=P_solar-plane.P_avi;
                elseif(P_solar>=plane.P_avi && h==environment.h_max)
                    if(P_solar>=P_elec_tot)
                        %Standard control: Only give level power
                        P_prop=P_level;
                    else
                        P_prop=P_solar-plane.P_avi; %Start sinking cause not enough solar energy for level flight
                    end
                   %P_prop=min(P_level/n_propulsion,P_solar-P_0);   %Either Plevel if we have enough, or if not enough solar energy then P_s-P_0 -> start sinking    
                else
                    P_prop=0;
                end

                delta_E=Dt*(P_solar-P_prop-plane.P_avi);

                if(h>h_before) AC_STATE=ACS.CLIMBING;
                elseif(h==h_before) AC_STATE=ACS.MAXALTREACHED;
                elseif(h<=h_before) AC_STATE=ACS.SINKING_PROPON;
                end
            %Case 2: at h_0: Begin climb or stay    
            else
                if(AC_STATE==ACS.CHARGING) AC_STATE=ACS.LEVELFLIGHT; end

                P_prop=P_level+climb*max(0,P_solar-P_level-plane.P_avi);
                delta_E=Dt*(P_solar-plane.P_avi-P_prop);                    
            end
        else
            %Case 3: Above std_Altitude: Sinking with propeller power=0->u=-Plevel...
            if h>environment.h_0 
                P_prop=0;

                if(P_solar-plane.P_avi)>=0 delta_E=0; 
                else delta_E=Dt*(-plane.P_avi); 
                end
            %Case 4:  At std-altitude: charge battery only OR do level flight at night!    
            else
                h=environment.h_0;
                P_prop=P_level;
                delta_E=Dt*(P_solar-P_prop-plane.P_avi);

                if (delta_E>=0)     AC_STATE=ACS.CHARGING;
                elseif (delta_E<0)  AC_STATE=ACS.DECHARGING;
                end

                if(MDflag==1 && AC_STATE==ACS.DECHARGING && t>12*3600) %return to multiday-simulation
                    display('debug');
                    t_excess=NaN;
                    t_endurance=NaN;
                    t_chargemargin=NaN;
                    flightdata=0; %only for now... is gonna be set later...
                    MDresults=[t h E_bat P_level+plane.P_avi t_eq t_sunrise]; % Return this current state to the multiday-sim
                    return;
                end

            end

        end
        % -------------------------------------------------------------------
        % STATE PROPAGATION
        % -------------------------------------------------------------------
        %SINK or CLIMB (or LEVEL)
        if(abs(P_prop-P_level)/P_level>10*eps)
            if(P_prop>P_level) etaclimb=plane.eta_climb;
            elseif(P_prop<P_level) etaclimb=1.0;
            end
            P_prop=min(P_prop,plane.P_prop_max); %Artificially limit climbing power
            h=h+Dt*etaclimb*climb*(P_prop-P_level)/(plane.m*g);
        end    
        if (h<environment.h_0) h=environment.h_0;
        elseif(h>environment.h_max) h=environment.h_max;
        end
        %charge or decharge?
        if delta_E>0
            delta_E=delta_E*plane.bat.eta_chrg;
        else
            delta_E=delta_E/plane.bat.eta_dchrg;
        end
        E_bat = min(E_bat+delta_E, E_bat_max);
        %time
        t = t + Dt;

        %Store flight data
        t_array=[t_array,t];
        h_array=[h_array,h];
        %Re_array=[Re_array,Re];
        %CL_array=[CL_array,CL];
        %CD_array=[CD_array,CD];
        %v_array=[v_array,v];
        bat_array=[bat_array,E_bat];
        irr_array=[irr_array,irr];
        P_solar_array=[P_solar_array,P_solar];
        P_elec_tot_array=[P_elec_tot_array,P_elec_tot];
        P_prop_array=[P_prop_array,P_prop];

    end
    %--------------------------------------------
    % MAIN SIMULATION LOOP END
    %--------------------------------------------
    t_end = t;
else % if the equilibrium is never reached:
    E_bat = E_bat_max/2; % symmetry assumption for the discharge
    t=12*3600; % symmetry assumption for the discharge
    while E_bat >= 0
        irr_vec=solar_radiation_on_surface2(0,0,mod(t,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,0, h,environment.albedo);
        irr = irr_vec(1);
        P_solar = environment.clearness*irr*plane.SolarModules.Surface * plane.SolarModules.eta;
        delta_E = Dt*(P_solar-P_elec_tot)/plane.bat.eta_dchrg;
        E_bat = E_bat + delta_E;
        t = t + Dt;
    end
    t_end = t;
    t_start_full = 24*3600-t;
end

%---------------------------------------------------------------------
% DONE! CALCULATE PEFORMANCE RESULTS!
%---------------------------------------------------------------------
if eq_reached_flag == 1
    if(t_end-24*3600<t_eq)
        t_excess=(t_end-(24*3600+t_eq))/3600;
    else
        t_excess = plane.bat.eta_dchrg*E_bat/P_elec_tot/3600;
    end
    if(isnan(t_fullycharged)~=1 && isnan(t_eq2)~=1)
        t_chargemargin=(t_eq2-t_fullycharged)/3600;
    else
        t_chargemargin=NaN;
    end
else
    t_excess = NaN;
    t_chargemargin=NaN;
end
t_endurance = max(24*3600-2*t_start_full,t_end-t_start_full)/3600;
if t_excess > 0
    t_endurance = NaN; % the term endurance non-existent
end
i_min_E_Bat=find(t_array==(86400+t_eq),1,'first');
if(isempty(i_min_E_Bat)) i_min_E_Bat=numel(t_array); end
min_E_Bat = bat_array(i_min_E_Bat);
min_SoC = min_E_Bat/E_bat_max;

%Save all the data
flightdata.t_array=t_array;
flightdata.h_array=h_array;
flightdata.Re_array=Re_array;
flightdata.CL_array=CL_array;
flightdata.CD_array=CD_array;
flightdata.v_array=v_array;
flightdata.bat_array=bat_array;
flightdata.irr_array=irr_array;
flightdata.P_solar_array=P_solar_array;
flightdata.P_elec_tot_array=P_elec_tot_array;
flightdata.P_prop_array=P_prop_array;
%Also set b,AR,m_bat for flightdata, but this has to be overwritten later!!
flightdata.b=NaN;
flightdata.AR=NaN;
flightdata.m_bat=NaN;

MDresults=0; %Set to zero cause not needed in single day simulation
%---------------------------------------------------------------------
% DEBUG OUTPUTS ONLY
%---------------------------------------------------------------------
if(DEBUG==1 && (environment.dayofyear==173 ))%|| environment.dayofyear==249 ||environment.dayofyear==252))
    str=strcat('24h flight overview: DayOfYear #',num2str(environment.dayofyear));
    figure('Name',str);
    ax(1)=subplot(3,1,1);
    dh_array=zeros(size(h_array));
    dh_array(1)=0;
    for i=2:size(h_array,2)
        dh_array(i)=h_array(i)-h_array(i-1);
    end
    plot(t_array/3600,h_array);
    legend('h');
    ylabel('Altitude above MSL [m]');

%    ax(end+1)=subplot(4,1,2);    
%         plot(t_array/3600,Re_array/100,...
%                 t_array/3600,CL_array*1000,...
%                 t_array/3600,CD_array*10000,...
%                 t_array/3600,dh_array*10);
%         legend('Re/100','Cl*1000','Cd*10000','dh*10');
    ax(end+1)=subplot(3,1,2);
    plotyy(t_array/3600,bat_array/1e6,t_array/3600,bat_array/E_bat_max*100.0);
    ylabel('Battery');
    legend('Battery Energy [J]', 'Battery SoC [%]');
    ax(end+1)=subplot(3,1,3);
    %plot(t_array/3600,irr_array,...
    plot(t_array/3600,P_solar_array,...
         t_array/3600,P_elec_tot_array,...
         t_array/3600,P_prop_array,...
         t_array/3600,plane.P_avi);
    ylabel('Power [W] or [W/m^2]');
    xlabel('Time of Day [h]');
    legend('P_{Solar}[W]','P_{electot}[W]','P_prop[W]','P_{avionics}[W]');
    linkaxes(ax,'x');
    %DEBUG=0;
end

if(DEBUG==1)
    str=strcat('d: ',num2str(environment.dayofyear),' t_exc: ', num2str(t_excess),' t_sr: ', num2str(t_sunrise), ' t_eq: ',num2str(t_eq), ' t_eq2: ', num2str(t_eq2));
    str=strcat(str, ' t_fc: ',num2str(t_fullycharged),'P_eltot: ',num2str(P_elec_tot));
    display(str);
end


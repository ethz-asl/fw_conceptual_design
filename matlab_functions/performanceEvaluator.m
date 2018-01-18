% *************************************************************************
%    performanceEvaluator.m: UAV Simulation & Performance Calculation 
% *************************************************************************
% Descr.: This is the core module that performs the forward-simulation
%   of the UAV in time by considering the airplane characteristics and the
%   the current environment. Also calculates the expected performance 
%   metrics. Simulates one-day or multi-day flight with different
%   initial conditions (see settings) and either only conceptually-designed 
%   aircraft (from AirplaneDesign.m) or fixed configurations (from 
%   AirplaneAnalysis.m) that can involve data input from flight tests.
% Authors: S. Leutenegger (2009), P. Oettershagen (2015)
% *************************************************************************
% Input Arguments:
% ****************
% params        technological parameters
% plane         airplane characteristics  
% environment:  environmental conditions 
%   .lat        latitude [deg]
%   .T_ground   temperature on the ground [°K]
%   .day        day of the month
%   .month      month of the year
%   .clearness  clearness of the sky [-]
%   .albedo     earth surface albedo [-]
% settings      simulation settings
% *************************************************************************

function [results, flightdata] = performanceEvaluator( params, plane, environment, settings)

% set the time discretization:
Dt = settings.dt; %[s]

%Setup data storage arrays
flightdata.b = plane.struct.b;
flightdata.AR = plane.struct.AR;
flightdata.m_bat = plane.bat.m;
flightdata.t_array = [];
flightdata.h_array = [];
flightdata.irr_array = [];
flightdata.eta_array = [];
flightdata.bat_array = [];
flightdata.P_solar_array = [];
flightdata.P_elec_tot_array = [];
flightdata.Re_array = [];
flightdata.CL_array = [];
flightdata.CD_array = [];
flightdata.v_array = [];
flightdata.P_prop_array = [];
flightdata.P_charge_array = [];

results.t_eq2 = NaN;
results.t_fullcharge = NaN;
results.t_chargemargin = NaN;
results.t_excess = NaN;
results.t_endurance = NaN;
results.min_SoC = NaN;
results.t_sunrise = NaN;
results.t_max = NaN;
results.t_sunset = NaN;

flags.sunrise_reached = 1;
flags.sunset_reached = 1;
flags.eq_reached = 1;                   % Flag: Power eq. reached or not
flags.eq_always_exceeded = 0;           % Floag: Power eq. always exceeded (happens e.g. far north in summer when sun does not set)?
flags.max_solarpower_reached = 1;
flags.full_charge_reached = 1;


% -------------------------------------------------------------------
% STEP 1: PREPARATORY CALCULATIONS
%         (e.g. finding sunrise and equality times, etc)
% -------------------------------------------------------------------
% Simple initializations
t = 0;
h = environment.h_0;
P_solar = 0;
irr_vec = zeros(1,9);

% Data cleaning
environment.dayofyear = round(environment.dayofyear);

%STEP 1a: Simple pre-calculations
E_bat_max = params.bat.e_density * plane.bat.m; % maximum battery energy state

if(~(exist('plane', 'var') && isfield(plane, 'ExpPerf') && isfield(plane.ExpPerf, 'solar') && isfield(plane.ExpPerf.solar, 'surface'))) 
    plane.ExpPerf.solar.surface = plane.struct.b^2/plane.struct.AR*params.solar.rWngCvrg;
end

%STEP 1b: Calculate air properties
[rho,mu] = findAir(environment.h_0,environment.T_ground);

%Step 1c: Calc Total electric power for level flight
if(exist('plane', 'var') && isfield(plane, 'ExpPerf') && isfield(plane.ExpPerf, 'P_prop_level'))
    P_elec_level = plane.ExpPerf.P_prop_level * (plane.m/plane.ExpPerf.m)^1.5 * sqrt(plane.ExpPerf.rho_P_prop_level/rho) * (1.0+environment.turbulence);
else
    P_elec_level = CalcPFromPolars(plane, params, rho, mu) * (1.0+environment.turbulence);
end
P_elec_level_tot_nom = P_elec_level + plane.avionics.power + plane.payload.power; %This is the nominal (not instantaneous) total electric level-flight power
results.P_elec_level_tot_nom = P_elec_level_tot_nom;

%Step 1d: Find t_sunrise, i.e. when P_solar=0.
while irr_vec(1) <= 0
    irr_vec = solar_radiation_on_surface2(0,0,mod(t+environment.add_solar_timeshift,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,environment.lon, h,environment.albedo);
    if (t==0 && irr_vec(1)>0) %Special case 1: Sun does not set during night
        flags.sunset_reached = 0;
        flags.sunrise_reached = 0;
        results.t_sunrise=0;
        break;
    elseif (t >= 3600*24) %Special case 2: Sun does not rise a all
        flags.sunrise_reached = 0;
        results.t_sunrise = 0;
        break;
    end
    results.t_sunrise = t;
    t = t+Dt;
end


%Step 1e: Find t_max, i.e. t(P_solar = P_solar_max)
if(~flags.sunset_reached)
    results.t_max = 12*3600; % TODO: This is in solar time at the moment, needs to be changed to local time however
else
    d_irr = 0;
    while d_irr >= 0
        irr_old = irr_vec(1);
        irr_vec = solar_radiation_on_surface2(0,0,mod(t+environment.add_solar_timeshift,86400),(environment.dayofyear+floor(t/86400)),...
                    environment.lat,environment.lon, h,environment.albedo);
        d_irr = irr_vec(1) - irr_old;
        if t >= 3600*15
            flags.max_solarpower_reached = 0; % the maximum solar power point was never reached
            results.t_max = NaN;
            break; % abort the search
        end
        results.t_max = t;
        t = t+Dt;
    end
end

%Step 1f: Find t_eq, that means the time when P_Solar=P_elect and thus the charge can start
t = results.t_sunrise;
while P_solar < P_elec_level_tot_nom
    [P_solar,irr_vec,~] = CalculateIncomingSolarPower(t,h,environment,plane,settings,params);
    if (t > results.t_max || ~flags.max_solarpower_reached)
        flags.eq_reached = 0; % the equilibrium point was never reached
        results.t_eq = NaN;
        break; % abort the search
    elseif (t==results.t_sunrise && P_solar > P_elec_level_tot_nom)
        flags.eq_always_exceeded = 1; % the equilibrium point is already exceeded at sunrise -> P_solar is never smaller P_elec_level_tot_nom!
        results.t_eq = NaN;
        break;
    end
    results.t_eq = t;
    t = t+Dt;
end

%Step 1g: Find t_eq2, that means the time when P_Solar=P_electot again
%        and thus the discharge would have to start if flying at h_0
if(flags.eq_reached && ~flags.eq_always_exceeded)
    t = results.t_max;
    while P_solar >= P_elec_level_tot_nom
        [P_solar,irr_vec,~] = CalculateIncomingSolarPower(t,h,environment,plane,settings,params);
        t = t+Dt;
    end
    results.t_eq2 = t;
else
    results.t_eq2 = NaN;
end

%Step 1h: Find sunset time
if(flags.sunset_reached)
    while irr_vec(1) > 0
        irr_vec = solar_radiation_on_surface2(0,0,mod(t+environment.add_solar_timeshift,86400),(environment.dayofyear+floor(t/86400)),...
                    environment.lat,environment.lon, h,environment.albedo);
        t = t+Dt;
    end
    results.t_sunset = t;
else
    results.t_sunset = NaN;
end

% -------------------------------------------------------------------
% STEP 2: Pre-simulate the day
% -------------------------------------------------------------------

if flags.eq_reached == 1
    E_bat = 0;
    if(~flags.eq_always_exceeded) t = results.t_eq;
    else t = 0;
    end
elseif (flags.eq_reached == 0)
    t = results.t_max;
    if(~flags.max_solarpower_reached) t = 12*3600; end
    E_bat = E_bat_max/2;
end

% Determine earliest start time when starting with full batt (step backward in time)
while E_bat < E_bat_max
    t = t - Dt;
    if t < results.t_sunrise;
        P_solar = 0; % set zero light before sunrise
    else
        [P_solar,irr_vec,~] = CalculateIncomingSolarPower(t,h,environment,plane,settings,params);
    end
    E_bat = E_bat + Dt * (P_elec_level_tot_nom-P_solar) / params.bat.eta_dchrg;
end
t_start_full = t;
    
%--------------------------------------------
% STEP 3: MAIN SIMULATION LOOP
%--------------------------------------------
h = environment.h_0;
h_before = h;

% Initial conditions (determined throgh settings-variable by caller-function)
if(settings.SimType == 0)
    if(flags.eq_reached == 1)
        E_bat = 0;
        if(~flags.eq_always_exceeded) t = results.t_eq;
        else t = 0;
        end
    else
        E_bat = E_bat_max;
        t = t_start_full;
    end
elseif(settings.SimType == 1)
    E_bat = settings.InitCond.SoC * plane.bat.m * params.bat.e_density;
    t = settings.InitCond.t;
end

t_sim_end = settings.SimTimeDays * 3600 * 24;
if(~isnan(results.t_eq)) t_sim_end = t_sim_end + results.t_eq; end
% Note: Important to have multiple of days + t_eq here! Sim needs to
% stop at X-days + t_eq, otherwise t_excess calculation could be wrong!

turb_before = environment.turbulence;
while E_bat >=0 && h>=environment.h_0 && t<=t_sim_end
    %--------------------------------------------
    % STEP 3a: Recalculate current flight conditions
    %--------------------------------------------

    tmod = mod(t,86400);
    
    % Pre-calculate potential level-power increase due to atmospheric turbulence, both in general as well as day-turbulence
    turb = environment.turbulence;
    t_dayturb_start = results.t_sunrise + 2.1*3600; %Heuristics only. Identified on AtlantikSolar 81h-flight
    t_dayturb_end = results.t_sunset - 1.5*3600;    %Heuristics only. Identified on AtlantikSolar 81h-flight
    t_dayturb_mid = 0.5*(t_dayturb_start+t_dayturb_end);
    if(tmod > t_dayturb_start && tmod < t_dayturb_end) %Implement linear symmetric turbulence distribution
        cur_day_turb = 2*environment.turbulence_day *(1-abs((tmod-t_dayturb_mid)/(t_dayturb_start-t_dayturb_mid)));
        turb = turb + cur_day_turb; 
    end
    
    %Step1b: Calc Total electric power for level flight
    if(turb ~= turb_before || h ~=h_before) 
        [rho,mu] = findAir(h,environment.T_ground); % Altitude has changed, recalc Re,CL,CD,v,Plevel,Pelectot
        h_before = h;
        
        if(exist('plane', 'var') && isfield(plane, 'ExpPerf') && isfield(plane.ExpPerf, 'P_prop_level'))
            P_elec_level = plane.ExpPerf.P_prop_level * (plane.m/plane.ExpPerf.m)^1.5 * sqrt(plane.ExpPerf.rho_P_prop_level/rho) * (1.0 + turb);
        else
            P_elec_level = CalcPFromPolars(plane, params, rho, mu) * (1.0 + turb);
        end
        P_elec_level_tot = P_elec_level + plane.avionics.power + plane.payload.power;
    end
    turb_before = turb;
    
    % Recalculate parameters for optimal/higher-speed cruise (not
    % implemented anymore)
    %if(params.optGRcruise==1 && AC_STATE==ACS.MAXALTREACHED)
    %    OptCruisePars=PreparePars_OptCruise(m_wo_bat+bat.m,rho,mu,A_wing,A_wing/b,g,polar);
    %end

    if(tmod<results.t_sunrise || tmod>results.t_sunset)
        P_solar = 0;
        P_solar_etas = [NaN NaN];
        irr_vec = NaN*ones(1,9);
    else
        [P_solar,irr_vec, P_solar_etas] = CalculateIncomingSolarPower(t,h,environment,plane,settings,params);
    end
    
    % -------------------------------------------------------------------
    % STEP 3b: CONTROL LAW. P_prop is the only control input.
    % -------------------------------------------------------------------
    % TODO for future: This can probably be simplified a lot.

    extraP = P_solar - P_elec_level - plane.avionics.power - plane.payload.power;
    if(E_bat >= E_bat_max)
        if(extraP >= 0 && h < environment.h_max) % Enough Power, Not Maximum Altitude
             P_prop = settings.climbAllowed * extraP + P_elec_level;
        elseif(extraP >= 0 && h >= environment.h_max) % Enough Power, Maximum Altitude
             P_prop = P_elec_level;
        elseif(extraP < 0 && h > environment.h_0) % Not Enough Power, Higher than Needed
             P_prop = max( 0 , P_elec_level + extraP );
        elseif(extraP < 0 && h <= environment.h_0) % Not Enough Power, Lower than Needed
             P_prop = P_elec_level;
         end
     else
         if(h>environment.h_0) % Not Fully Charged, Altitude Higher than Needed
             P_prop = 0;
         else % Not Fully Charged, Altitude Lower than Needed
             P_prop = P_elec_level;
          end
      end
    
    %  %Not implemented anymore (For : Enough Power, Maximum Altitude)                           
    %  if(params.optGRcruise==1) 
    %  [vtmp,Retmp,CLtmp,CDtmp]=CalcFlightPars_OptCruise((P_solar-P_0)*n_propulsion,P_level,m_wo_bat+bat.m,rho,mu,A_wing,A_wing/b,g,polar,OptCruisePars);
    %  if(isnan(Retmp)==0) %because CalcFlightPars can fail when P_prop->P_level due to inaccuracies in the interpolation 
    %  v=vtmp;Re=Retmp;CL=CLtmp;CD=CDtmp;
    %  P_prop=P_solar-P_0; %Max. Power Cruise at level flight
    %  end
    %  end

    % -------------------------------------------------------------------
    % STEP 3c: STATE PROPAGATION
    % -------------------------------------------------------------------
    
    % Step 3c1): Calculate the deltas of the CURRENT TIME STEP of the state variables
    % dh/dt : Sink, Climb, or level?
    if(abs(P_prop-P_elec_level)/P_elec_level > 10*eps)
        if(P_prop>P_elec_level) etaclimb = params.prop.eta_climb; % When climbing, efficieny may be different than when sinking
        elseif(P_prop<P_elec_level) etaclimb = 1.0;
        end
        P_prop = min(P_prop, plane.prop.P_prop_max); %Artificially limit climbing power
    end
    % dE/dt: Charge or discharge?
    P_elec_tot = P_prop + plane.avionics.power + plane.payload.power;
    delta_E  = Dt * (P_solar - P_elec_tot);
    if delta_E>0
        delta_E = ChargeLimiter(E_bat/E_bat_max, delta_E/Dt * params.bat.eta_chrg, params, plane)*Dt;
    else
        delta_E = delta_E/params.bat.eta_dchrg;
    end
    
    % Step 3c2): Store flight data of the CURRENT TIME STEP (i.e. including current deltas, but NOT the updated state variables yet)
    flightdata.t_array = [flightdata.t_array, t];
    flightdata.h_array = [flightdata.h_array, h];
    flightdata.bat_array = [flightdata.bat_array, E_bat];
    flightdata.irr_array = [flightdata.irr_array, [irr_vec(1) ; irr_vec(2) ; irr_vec(7)]];
    flightdata.eta_array = [flightdata.eta_array, [P_solar_etas(1) ; P_solar_etas(2)]];
    flightdata.P_solar_array = [flightdata.P_solar_array, P_solar];
    flightdata.P_elec_tot_array = [flightdata.P_elec_tot_array, P_elec_tot];
    flightdata.P_prop_array = [flightdata.P_prop_array, P_prop];
    flightdata.P_charge_array = [flightdata.P_charge_array, delta_E/Dt];
    %flightdata.Re_array = [flightdata.Re_array, Re];
    %flightdata.CL_array = [flightdata.CL_array, CL];
    %flightdata.CD_array = [flightdata.CD_array, CD];
    %flightdata.v_array = [flightdata.v_array, v];
    
    % Step3: Calculate state variables of the NEXT TIME STEP
    t = t + Dt;
    E_bat = min(E_bat+delta_E, E_bat_max);
    if(abs(P_prop-P_elec_level)/P_elec_level > 10*eps)
        h = h + Dt*etaclimb*settings.climbAllowed*(P_prop-P_elec_level)/((plane.m_no_bat+plane.bat.m)*params.physics.g);
        if (h<environment.h_0) h=environment.h_0;
        elseif(h>environment.h_max) h=environment.h_max;
        end
    end
end
%--------------------------------------------
% MAIN SIMULATION LOOP END
%--------------------------------------------
t_end = t;

%---------------------------------------------------------------------
% STEP 4: Calculate Performance Results
%---------------------------------------------------------------------
% Note that all results.XXX fields have been initialized to NaN before.
if flags.eq_reached == 1
    % Step 4.1: Excess Time & min-SoC
    if(flags.eq_always_exceeded)
        results.t_excess = Inf();
        results.min_SoC = 1;
    %elseif(t_end-24*3600 < results.t_eq)
    elseif(mod(t_end,86400) < results.t_eq || mod(t_end,86400) > results.t_eq2)
        % Case a) We have not made it through the night (negative t_exc !)
        if(mod(t_end,86400) < results.t_eq) results.t_excess = (mod(t_end,86400)-results.t_eq) / 3600;
        elseif(mod(t_end,86400) > results.t_eq2) results.t_excess = (mod(t_end,86400)-86400-results.t_eq) / 3600;
        end
    else
        % Case b) We have made it through the night (positive t_exc!). We do 
        % extensive searches here to find not any but the correct (of  the 
        % last day) t_exc, so we basically do backwards search.
        idx_lastteq = find(mod(flightdata.t_array-results.t_eq,86400) < Dt,1,'last');
        if(isempty(idx_lastteq)) 
            display('WARNING no t_eq for calculation of t_exc found'); 
        end
        results.t_excess = flightdata.bat_array(idx_lastteq) * params.bat.eta_dchrg / P_elec_level_tot_nom / 3600;
        results.min_SoC = flightdata.bat_array(idx_lastteq)/E_bat_max;
    end
    % Step 4.2: Charge margin
    if(flags.eq_always_exceeded)
        results.t_chargemargin = Inf();
    else    
        % We do extensive searches here to find not any but the correct (of 
        % the last day) charge margin, so we basically do backwards search.
        idx_current = numel(flightdata.t_array);
        idx_lastfullcharge_end = find((flightdata.bat_array == E_bat_max) & (flightdata.t_array > (settings.SimTimeDays-1)*86400),1,'last');
        %if(flightdata.t_array(idx_lastfullcharge_end) < (settings.SimTimeDays-1)*86400) idx_lastfullcharge_end=[]; end;
        idx_lastteqbeforelastfullcharge = find(mod(flightdata.t_array(1:idx_lastfullcharge_end)-results.t_eq,86400) < Dt,1,'last');
        if(isempty(idx_lastteqbeforelastfullcharge)) idx_lastteqbeforelastfullcharge=1; end %Handle case of not t_eq before
        idx_lastfullcharge_start = idx_lastteqbeforelastfullcharge - 1 + find (flightdata.bat_array(idx_lastteqbeforelastfullcharge:idx_lastfullcharge_end) == E_bat_max,1,'first');
        if(isempty(idx_lastfullcharge_end) || isempty(idx_lastfullcharge_start))
            results.t_fullcharge=NaN;
            results.t_fullcharge=NaN;
        else
            results.t_fullcharge = mod(flightdata.t_array(idx_lastfullcharge_start),86400);
            results.t_chargemargin = (flightdata.t_array(idx_lastfullcharge_end+1)-flightdata.t_array(idx_lastfullcharge_start))/3600.0;
        end
    end
end

% Step 4.3: Endurance
if (results.t_excess < 0 || isnan(results.t_excess))
    results.t_excess = 0;
    results.min_SoC = 0;
    if(settings.SimType==0) results.t_endurance = (t_end-t_start_full)/3600; %results.t_sunrise)/3600; %WRONG when defined as in [LeuteneggerJIRS]
    elseif(settings.SimType==1) results.t_endurance = (t_end-settings.InitCond.t)/3600;
    end
end

%---------------------------------------------------------------------
% DEBUG OUTPUTS ONLY
%---------------------------------------------------------------------
if(settings.DEBUG==1)
    %Plot_BasicSimulationTimePlot(flightdata,environment,params, plane);
    Plot_BasicSimulationTimePlot(flightdata,environment,params, plane);
end

if(settings.DEBUG==1)
    str=strcat('d: ',num2str(environment.dayofyear),' t_exc: ', num2str(results.t_excess),' t_endur: ', num2str(results.t_endurance),' t_sr: ', num2str(results.t_sunrise), ' t_eq: ',num2str(results.t_eq), ' t_eq2: ', num2str(results.t_eq2));
    str=strcat(str, ' t_fc: ',num2str(results.t_fullcharge),' P_eltot: ',num2str(P_elec_tot),' Min. SoC: ',num2str(results.min_SoC));
    display(str);
end


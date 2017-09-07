% *************************************************************************
%        AirplaneAnalysis.m: Solar-powered UAV Conceptual Analysis
% *************************************************************************
% Descr.: Use this file to perform detailed analysis of your solar- 
%   powered UAV, i.e. analyse a selected configuration w.r.t operation
%   at different days of the year, different latitudes, or different
%   atmospheric clearness (e.g. clouds) and turbulence (e.g. wind) 
%   values.
% Authors: P. Oettershagen, S. Leutenegger (2009-2015), based on A. Noth
% *************************************************************************

clear variables
close all
addpath(genpath('matlab_functions')) 

% -------------------------------------------------------------------------
% STEP 1: ANALYSIS SETUP
% -------------------------------------------------------------------------

% Analysis-variables have to be specified here! The choices are VAR.DAY_OF_YEAR,
% VAR.CLEARNESS, VAR.TURBULENCE, VAR.LATITUDE, VAR.POWER (all which are dfined
% in the file Var.m). If you only want to specify one or two variables, simply 
% provide a constant value for the remaining (i.e. the second/third) one
%P-Std
% vars(1) = VAR.DAY_OF_YEAR;
% vars(1).values = [5*365/12+21 5*365/12+30 6*365/12+15];%[5*30.4166+21 : 30.55/4 : 7*30.5+21;]
% vars(2) = VAR.LATITUDE;
% vars(2).values = 47.6;
% vars(3) = VAR.TURBULENCE;
% vars(3).values = 0;

%P2
vars(1) = VAR.DAY_OF_YEAR;
vars(1).values = [5*365/12+21:10:11*365/12+21];%[5*30.4166+21 : 30.55/4 : 7*30.5+21;]
vars(2) = VAR.LATITUDE;
vars(2).values = 0:30:90;
vars(3) = VAR.TURBULENCE;
vars(3).values = 0;

% Airplane general technological parameters first
initParameters;
params.bat.chrg_lim_type = 2;   % Enable charge limiting using experimental data

%This is the plane-specific data
plane.payload.mass = 0.0;
plane.bat.m = 2.918+0*0.293+0*0.04863;
plane.m_no_bat = 6.92-2.918 + plane.payload.mass;
plane.struct.b = 5.65;
plane.struct.AR = 18.5;
plane.m = plane.m_no_bat+plane.bat.m;
plane.ExpPerf.m = 6.92;    %
plane.ExpPerf.solar.surface = 88 * (0.125^2 - 4*70.36E-6);
plane.ExpPerf.P_prop_level = 35.8;          %
plane.ExpPerf.rho_P_prop_level = 1.10;     % Density at which power curve of aircraft was recorded
plane.avionics.power = 6.0;
plane.payload.power = 0;
plane.prop.P_prop_max = 180.0;
    
%These are the (default) environment parameters
environment.dayofyear = 5*30.5+30;          % Only used if you do not specify this as a VARIABLE above
environment.lat = 47.6;                     % Rafz
environment.lon = 8.53;
environment.h_0 = 416+120;                  % with 120m AGL flight altitude for enough safety
environment.h_max = 700;
environment.T_ground = 31.3+273.15;
environment.turbulence = 0;                 % Only used if you do not specify this as a VARIABLE above
environment.turbulence_day = 0.0;%0.307;    % Relative increase of power consumption during the day, e.g. due to thermals
environment.clearness = 1.0;                % Only used if you do not specify this as a VARIABLE above
environment.albedo = 0.12;
environment.add_solar_timeshift = -3600;    % [s], due to Daylight Saving Time (DST), actually used for solar income calculations
environment.plot_solar_timeshift = -1.533;  % [h], just used for plotting results (to plot them in solar time), does not affect anything else

%Evaluation settings
settings.DEBUG = 0;                         % Force DEBUG mode
settings.dt = 100;                          % Discretization time interval [s]
settings.climbAllowed = 0;
settings.SimType = 0;                       % 0 = Start on t_eq, 1 = start on specified Initial Conditions
settings.SimTimeDays = 2;                   % Simulation Time in days (e.g. 1 = std. 24h simulation)
settings.InitCond.SoC = 0.577;%0.46;               % State-of-charge [-]
settings.InitCond.t = (11+14.0/60.0)*3600.0;%9.0*3600 + 32*60;     % [s]launch time
settings.useAOI = 1;                        % 1 to enable the use of angle-of-incidence dependent solar module efficiency
settings.useDirDiffRad = 1;                 % 1 to enable the use of separate diffuse and direct radiation solar module efficiencies

% -------------------------------------------------------------------------
% STEP 2: Calculate analysis results
% -------------------------------------------------------------------------

% Number of configurations calculated
N = numel(vars(1).values) * numel(vars(2).values) * numel(vars(3).values);
disp(['Number of configurations to be calculated: ' num2str(N)]);
h=waitbar(0,'Progress');

str='';
ctr = 0;
for i = 1:numel(vars(3).values)
    for k = 1:numel(vars(2).values)
        for j = 1:numel(vars(1).values)
            
            ctr = ctr+1;
            varval(3)=vars(3).values(i);
            varval(2)=vars(2).values(k);
            varval(1)=vars(1).values(j);
            
            %Assign variables dynamically
            idx = find(vars == VAR.CLEARNESS,1,'first');
            if ~isempty(idx) ; environment.clearness = varval(idx); end
            idx = find(vars == VAR.TURBULENCE,1,'first');
            if ~isempty(idx) ; environment.turbulence = varval(idx); end
            idx = find(vars == VAR.DAY_OF_YEAR,1,'first');
            if ~isempty(idx) ; environment.dayofyear = varval(idx); end
            idx = find(vars == VAR.LATITUDE,1,'first');
            if ~isempty(idx) ; environment.lat = varval(idx); end
            idx = find(vars == VAR.POWER,1,'first');
            if ~isempty(idx) ; plane.ExpPerf.P_prop_level = varval(idx); end
            
            %Execute simulation
            [results(i,k,j), flightdata(i,k,j)] = performanceEvaluator(params ,plane, environment, settings);
            
            if(fabs(environment.plot_solar_timeshift) > 0.01)
                results(i,k,j).t_eq2 = results(i,k,j).t_eq2 + environment.plot_solar_timeshift * 3600;
                results(i,k,j).t_fullcharge = results(i,k,j).t_fullcharge + environment.plot_solar_timeshift * 3600;
                results(i,k,j).t_sunrise = results(i,k,j).t_sunrise + environment.plot_solar_timeshift * 3600;
                results(i,k,j).t_max = results(i,k,j).t_max + environment.plot_solar_timeshift * 3600;
                results(i,k,j).t_sunset = results(i,k,j).t_sunset + environment.plot_solar_timeshift * 3600;
                results(i,k,j).t_eq = results(i,k,j).t_eq + environment.plot_solar_timeshift * 3600;
            end 
            
            completedRatio = ((i-1)*numel(vars(2).values)*numel(vars(1).values) + (k-1)*numel(vars(1).values) + j)/N;
            waitbar(completedRatio,h,[num2str(completedRatio*100.0,'Progress: %.0f\n') '%']);
            
            str = [str sprintf('#%d| Set: DoY=%g,Lat=%g,P=%g,CLR=%g,Turb=%g   Res:Soc_min=%g%%,T_exc=%gh,T_cm=%gh,T_end=%gh   CharTimes:t_sr=%gh t_eq1=%gh t_fc=%gh t_fc90=NA t_eq2=%gh t_ss=%gh\n',ctr,...
                environment.dayofyear,environment.lat,plane.ExpPerf.P_prop_level,environment.clearness,environment.turbulence,...
                results(i,k,j).min_SoC*100,results(i,k,j).t_excess,results(i,k,j).t_chargemargin,results(i,k,j).t_endurance,...
                results(i,k,j).t_sunrise/3600,results(i,k,j).t_eq/3600,results(i,k,j).t_fullcharge/3600, results(i,k,j).t_eq2/3600, results(i,k,j).t_sunset/3600)];
        end
    end
end
close(h)

display('*** Performance solutions ***');
display(str);

% -------------------------------------------------------------------------
% STEP 3: PLOTTING
% -------------------------------------------------------------------------
% Note: Plotting scripts are located in the matlab_functions/PlotScripts
% folder. Please modify and call these scripts if you want to modify the 
% plots

Plot_AirplaneAnalysis_Standard(results, vars);
%Plot_AirplaneAnalysis_ASFinalPaper_PlotOrderChanged(results, [], environment, plane, params, flightdata, vars,1);

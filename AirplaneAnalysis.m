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
vars(1) = VAR.DAY_OF_YEAR;
vars(1).values = round(5*30.5+21);%:5:round(6*30.5+26);
vars(2) = VAR.CLEARNESS;
vars(2).values = 1.0;%0.5:0.25:1;
vars(3) = VAR.TURBULENCE;
vars(3).values = 0;

% Airplane general technological parameters first
initParameters;
params.bat.chrg_lim_type = 2;   % Enable charge limiting using experimental data

%This is the plane-specific data
plane.m_no_bat = 6.92-2.918;
plane.struct.b = 5.65;
plane.struct.AR = 18.5;
plane.bat.m = 2.918+0*0.293+0*0.04863;
plane.m = plane.m_no_bat+plane.bat.m;
plane.ExpPerf.m = plane.m_no_bat+plane.bat.m;    %
plane.ExpPerf.solar.surface = 88 * (0.125^2 - 4*70.36E-6);
plane.ExpPerf.P_prop_level = 37;          %
plane.ExpPerf.rho_P_prop_level = 1.095;     % Density at which power curve of aircraft was recorded
plane.avionics.power = 5.5;
plane.payload.power = 0;
plane.prop.P_prop_max = 180.0;
    
%These are the environment parameters
environment.dayofyear = round(5*30.5+21);
environment.lat = 47.6;                     % Rafz
environment.lon = 8.5;
environment.h_0 = 416+120;                  % with 120m AGL flight altitude for enough safety
environment.h_max = 700;                   %
environment.T_ground = 28+271.15;
environment.turbulence = 0;
environment.turbulence_day = 0.0;           % Relative increase of power consumption during the day, e.g. due to thermals
environment.clearness = 1.0;
environment.albedo = 0.12;
environment.add_solar_timeshift = -3600;    % [s], due to Daylight Saving Time (DST)

%Evaluation settings
settings.DEBUG = 0;                         % Force DEBUG mode
settings.dt = 200;                          % Discretization time interval [s]
settings.climbAllowed = 1;
settings.SimType = 0;                       % 0 = Start on t_eq, 1 = start on specified Initial Conditions
settings.SimTimeDays = 2;                   % Simulation Time in days (e.g. 1 = std. 24h simulation)
settings.InitCond.SoC = 0.46;               % State-of-charge [-]
settings.InitCond.t = 9.0*3600 + 32*60;     % [s]launch time
settings.useAOI = 1;                        % 1 to enable the use of angle-of-incidence dependent solar module efficiency
settings.useDirDiffRad = 1;                 % 1 to enable the use of separate diffuse and direct radiation solar module efficiencies

% -------------------------------------------------------------------------
% STEP 2: Calculate analysis results
% -------------------------------------------------------------------------

% Number of configurations calculated
N = numel(vars(1).values) * numel(vars(2).values) * numel(vars(3).values);
disp(['Number of configurations to be calculated: ' num2str(N)]);
h=waitbar(0,'Progress');

for i = 1:numel(vars(3).values)
    for k = 1:numel(vars(2).values)
        for j = 1:numel(vars(1).values)
            
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
            [results(i,k,j), ~] = performanceEvaluator(params ,plane, environment, settings);
            
            completedRatio = ((i-1)*numel(vars(2).values)*numel(vars(1).values) + (k-1)*numel(vars(1).values) + j)/N;
            waitbar(completedRatio,h,[num2str(completedRatio*100.0,'Progress: %.0f\n') '%']);
        end
    end
end
close(h)

% -------------------------------------------------------------------------
% STEP 3: PLOTTING
% -------------------------------------------------------------------------
% Note: Plotting scripts are located in the matlab_functions/PlotScripts
% folder. Please modify and call these scripts if you want to modify the 
% plots

Plot_AirplaneAnalysis_Standard(results, vars);
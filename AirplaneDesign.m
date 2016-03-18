% *************************************************************************
%          AirplaneDesign.m: Solar-powered UAV Conceptual Design
% *************************************************************************
% Descr.: Use this file to perform the conceptual design of your solar- 
%   powered UAV, i.e. analyse its performance (in the form of excess time,
%   charge margin, endurance and minimum battery state-of-charge) as a
%   function of its design variables (wing span b, aspect ratio AR, battery
%   mass m_bat). One can also design configurations considering different
%   atmospheric clearness (e.g. clouds) and turbulence (e.g. wind) values.
% Authors: P. Oettershagen, S. Leutenegger (2009-2015), based on A. Noth
% *************************************************************************

% Initialize
clear variables;
close all;
clc;
addpath(genpath('matlab_functions')) 

% -------------------------------------------------------------------------
% STEP 1: DESIGN SETUP
% -------------------------------------------------------------------------
% Set the three variables to choose as design variables here. Choices are the
% labels defined in the file VAR.m (i.e. VAR.WING_SPAN, VAR.BATTERY_MASS,
% VAR.ASPECT_RATIO, VAR.CLEARNESS and VAR.TURBULENCE, VAR.DAY_OF_YEAR, VAR.LATITUDE). 
%
% There is basically two design ways:
% 1. Specify wing span, battery mass and aspect ratio ranges to design your 
%    airplane first
% 2. Then (optionally) choose the optimal wing span and aspect ratio from the 
%    first step, and set 'VAR.BATTERY_MASS','VAR.CLEARNESS' and 'VAR.TURBULENCE'
%    to optimize the partially-fixed configuration in more detail over the 
%    remaining variables.
%
% Example:
% vars(1)= VAR.WING_SPAN;
% vars(1).values = 3:1:5; %Analyse over wing spans from 3 to 5m in 1m steps

%P1
vars(1) = VAR.WING_SPAN;
vars(1).values = 3.5:0.1:6.5;
vars(2) = VAR.BATTERY_MASS;
vars(2).values = 1.0:0.1:7.0;
vars(3) = VAR.ASPECT_RATIO;
vars(3).values = 18.5;

%P2
% vars(1) = VAR.DAY_OF_YEAR; %VAR.BATTERY_MASS;
% vars(1).values = floor(0*30.5):5:floor(11*30.5+29); %2.7:0.1:3.0;
% vars(2) = VAR.LATITUDE; %VAR.WING_SPAN;
% vars(2).values = 0:2.5:70; %5.4:0.1:5.8;
% vars(3) = VAR.ASPECT_RATIO;
% vars(3).values = 18.5;

% P3
% vars(1) = VAR.CLEARNESS; %VAR.BATTERY_MASS;
% vars(1).values = 0.4:0.025:1;
% vars(2) = VAR.TURBULENCE; %VAR.WING_SPAN;
% vars(2).values = 0.0:0.025:0.6;; %5.4:0.1:5.8;
% vars(3) = VAR.DAY_OF_YEAR;
% vars(3).values = [floor(3*30.5+21), floor(5*30.5+21)];

% Airplane general technological parameters first
initParameters;
params.structure.corr_fact = 1.21;    % Structural mass correction factor. Set to 
                                      % * 1.0 to use the original model without correction.
                                      % * 1.21 to correspond to AtlantikSolar initial structural mass calculation by D. Siebenmann

% This is the default configuration for our design variables! 
% (which is only used if we don't design over b, m_bat or AR)
plane.struct.b = 5.6;
plane.struct.AR = 18.5;
plane.bat.m = 2.9;

%This is the other plane-specific data. 
plane.avionics.power = 5.5;
plane.avionics.mass = 0.6;
plane.payload.power = 0;
plane.payload.mass = 0.0;
plane.prop.P_prop_max = 180.0;

% Set environment
environment.dayofyear = 5*30.5+21;
environment.lat = 47.6;                     % 1: Barcelona 2:Tuggen/CH
environment.lon = 8.5;
environment.h_0 = 416+120;                  % with 120m AGL flight altitude for enough safety
environment.h_max = 800;                   % Barcelona: 4000ft
environment.T_ground = 25+271.15;
environment.turbulence = 0;
environment.turbulence_day = 0.0;           % Relative increase of power consumption during the day, e.g. due to thermals
environment.clearness = 1.0;
environment.albedo = 0.12;
environment.add_solar_timeshift = -3600;    % [s], due to Daylight Saving Time (DST)

%Evaluation settings
settings.DEBUG = 0;                         % Force DEBUG mode
settings.dt = 200;                          % Discretization time interval [s]
settings.climbAllowed = 0;
settings.SimType = 0;                       % 0 = Start on t_eq, 1 = start on specified Initial Conditions
settings.SimTimeDays = 2;                   % Simulation Time in days (e.g. 1 = std. 24h simulation)
settings.InitCond.SoC = 0.46;               % State-of-charge [-]
settings.InitCond.t = 9.0*3600 + 32*60;     % [s]launch time
settings.evaluation.findalt = 0;            % if 1, it finds the maximum altitude for eternal flight
%settings.optGRcruise       =  0;           % 1 to allow cruise at optimal glide ratio & speed when max altitude reached 
settings.useAOI = 0;                        % 1 to enable the use of angle-of-incidence dependent solar module efficiency
settings.useDirDiffRad = 0;                 % 1 to enable the use of separate diffuse and direct radiation solar module efficiencies

% -------------------------------------------------------------------------
% STEP 2: Calculate performance results
% -------------------------------------------------------------------------

% Number of configurations calculated
N = numel(vars(1).values) * numel(vars(2).values) * numel(vars(3).values);
disp(['Number of configurations to be calculated: ' num2str(N)]);
h=waitbar(0,'Progress');

% Calculate performance results
for i = 1:numel(vars(3).values)
    for k = 1:numel(vars(2).values)
        for j = 1:numel(vars(1).values)
            
            varval(3)=vars(3).values(i);
            varval(2)=vars(2).values(k);
            varval(1)=vars(1).values(j);
            
            %Assign variables dynamically
            idx = find(vars == VAR.WING_SPAN,1,'first');
            if ~isempty(idx) ; plane.struct.b = varval(idx); end
            idx = find(vars == VAR.BATTERY_MASS,1,'first');
            if ~isempty(idx) ; plane.bat.m = varval(idx); end
            idx = find(vars == VAR.ASPECT_RATIO,1,'first');
            if ~isempty(idx) ; plane.struct.AR = varval(idx); end
            idx = find(vars == VAR.CLEARNESS,1,'first');
            if ~isempty(idx) ; environment.clearness = varval(idx); end
            idx = find(vars == VAR.TURBULENCE,1,'first');
            if ~isempty(idx) ; environment.turbulence = varval(idx); end
            idx = find(vars == VAR.DAY_OF_YEAR,1,'first');
            if ~isempty(idx) ; environment.dayofyear = varval(idx); end
            idx = find(vars == VAR.LATITUDE,1,'first');
            if ~isempty(idx) ; environment.lat = varval(idx); end
            
            [PerfResults(i,k,j),DesignResults(i,k,j),flightdata(i,k,j)] = ...
               evaluateSolution(plane,environment,params,settings);
           
            completedRatio = ((i-1)*numel(vars(2).values)*numel(vars(1).values) + (k-1)*numel(vars(1).values) + j)/N;
            waitbar(completedRatio,h,[num2str(completedRatio*100.0,'Progress: %.0f\n') '%']);
        end
    end
end
close(h)

% -------------------------------------------------------------------------
% STEP 3: Plotting
% -------------------------------------------------------------------------
% Note: Plotting scripts are located in the matlab_functions/PlotScripts
% folder. Please modify and call these scripts if you want to modify the 
% plots

Plot_AirplaneDesign_ASFinalPaper(PerfResults, DesignResults, environment, plane, params, flightdata, vars);
%Plot_AirplaneDesign_ASFinalPaper(PerfResults, DesignResults, environment, plane, params, flightdata, vars);

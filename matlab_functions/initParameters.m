%==========================================================================
%=== Global Design of solar-powered airplanes
%=== Initialization of technological parameters
%=== Stefan Leutenegger and Philipp Oettershagen
%=== 2009-2015
%==========================================================================

clear params;

%==========================================================================
%=== TECHNOLOGICAL PARAMETERS GROUP =======================================
%==========================================================================
%--- ADAPT THESE params -----------------------------------------------

%=== Propulsion Group =====================================================
params.prop.k_prop          = 0.0011;   % Mass/Power [kg/W]
params.prop.eta_ctrl        = 0.97;     % Eff. of motor controller [-]
params.prop.eta_mot         = 0.8;      % Efficiency of motor [-]
params.prop.eta_grb         = 0.97;     % Efficiency of gearbox [-]
params.prop.eta_plr         = 0.78;     % Efficiency of propeller [-]
params.prop.eta_climb       = 0.5;      % Additional loss factor when climbing [-]
params.prop.number          = 1;        % Number of propulsion units [-]

%=== Battery and Step-down Converter ======================================
params.bat.eta_chrg        = 0.95;     % Eff. of charge process [-]
params.bat.eta_dchrg       = 0.95;     % Eff. of discharge process [-]
params.bat.eta_bec         = 0.65;     % Eff. of bec (5V stepdown) [-]
params.bat.e_density       = 243*3600; % Energy density of battery including battery-cables sensors etc.
params.bat.distr           = 1;        % Distributed battery mass? 1/0
params.bat.chrg_lim_type   = 1;        % 0 = Don't use, 1=use experim. charge limit curve, 2 = use params below 
params.bat.chrg_lim_SoC1   = 0.9;      % SoC [-] above which charge limiting applies
params.bat.chrg_lim_Prel1  = 0.5;      % Maximum charge relative to battery capacity (i.e. in C's) when NOT limiting (i.e. below SoC1)
params.bat.chrg_lim_Prel2  = 0.02;     % Specific charge termination power, i.e. P_chrg(SoC=100%) = chrg_lim_Prel2 * e_bat

%=== Solar cells ==========================================================
params.solar.k_sc           = 0.4;      % Mass density of sc [Kg/m2]. ESTIMATED based on 200microns cells
params.solar.k_enc          = 0.2;      % Mass dens. of encaps. and cables! [Kg/m2]
params.solar.k_mppt         = 1/2368;   % Mass/Power of mppt [kg/W]
params.solar.rWngCvrg       = 0.9484;   % relative solar cell wing coverage area/without winglet area
params.solar.eta_sc         = 0.223;    % Efficiency of solar cells (more precisely: modules) [-]
params.solar.eta_cbr        = 0.95;     % Eff. of cambered conf. [-]
params.solar.eta_mppt       = 0.96;     % Efficiency of mppt [-]
params.solar.k_temp         = 0.003;    % Solar power (and thus eta) reduction in [1/K]

%==== Structure ===========================================================
params.structure.shell      = 0;        % 1 for shell wing, 0 for rib wing
params.structure.corr_fact  = 1.0;      % Correction factor that can be used to correct the structural mass estimation
                                        % results if an experimental data point (i.e. a finished structure) is available.

%--------------------------------------------------------------------------

%=== Physics ==============================================================
params.physics.g            = 9.81;      % Earth Acceleration [m/s^2]

% TODO check still needed? Or put (default) settings here?
%=== Performance evaluation
params.evaluation.clmb      = 1;         % 1 to allow altitude changes
params.evaluation.findalt   = 0;         % if 1, it finds the maximum
                                           % altitude for eternal flight
                                           
%==========================================================================
%=== AIRPLANE CHARACTERISTICS (Defaults only) =============================
%==========================================================================
%--- ADAPT THESE params -----------------------------------------------

%==== Avionics ============================================================
plane.avionics.mass        = 0.4;      % Mass of contr./electr/harness. [kg]
plane.avionics.power       = 5.5;      % Power required for control [W]

%==== Payload =============================================================
plane.payload.mass         = 0.0;      % Mass of payload [kg]
plane.payload.power        = 3.0;      % Power required for payload [W]

%==== Propulsion =============================================================
plane.prop.P_prop_max      = 180.0;    % Max. prop/climbing power




environment.usemars        = 0;

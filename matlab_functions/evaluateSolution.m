% *************************************************************************
%                         evaluateSolution.m 
% *************************************************************************
% Descr.: Performs the actual (structural) design of an airplane as a fct.
%   of its wing span, aspect ratio and battery mass. Then evaluates this 
%   configuration and returns the evaluated performance results
% Authors: S. Leutenegger (2009), P. Oettershagen (2015)
% *************************************************************************

function [PerfResults,DesignResults,flightdata] ...
 = evaluateSolution(plane,environment,params,settings)

% Input processing
% ================

% set the masses of the various components
% ========================================

% Solar module
DesignResults.m_solar = plane.struct.b^2/plane.struct.AR*params.solar.rWngCvrg*...
    (params.solar.k_sc+params.solar.k_enc);

% MPPT
I_max = 1000; %[W/m^2]
DesignResults.m_mppt = params.solar.k_mppt*I_max*plane.struct.b^2/plane.struct.AR ...
	*params.solar.eta_sc*params.solar.eta_cbr;

% total point mass
DesignResults.m_central = DesignResults.m_mppt + plane.payload.mass+plane.avionics.mass...
    +(1-params.bat.distr)*plane.bat.m;

% total distributed mass
DesignResults.m_distr = DesignResults.m_solar + plane.bat.m * params.bat.distr;

% Automatic Preliminary Design (S. Leutenegger / M. Jabas)
% ========================================================
if(params.structure.shell==1)
    [DesignResults.m_struct,masses,thicknesses,velocities,plane.polar]=StructureDesigner(plane.struct.b,...
        plane.struct.AR, DesignResults.m_central,DesignResults.m_distr,params.propulsion.number,environment.usemars,params.structure.corr_fact);
else
    [DesignResults.m_struct,masses,thicknesses,velocities,plane.polar]=StructureDesignerRibWing(plane.struct.b,...
        plane.struct.AR, DesignResults.m_central,DesignResults.m_distr,params.prop.number,environment.usemars,params.structure.corr_fact);
end
%Save structural masses to DesignResults (array) for later use
DesignResults.m_prop = masses.m_prop;
DesignResults.m_bat = plane.bat.m;
DesignResults.m_no_bat = DesignResults.m_distr + DesignResults.m_central + DesignResults.m_struct + DesignResults.m_prop - plane.bat.m;
plane.m_no_bat = DesignResults.m_no_bat;


% Now simulate the course of a day
% ================================

% Calculate the "Performance" in terms of endurance/range(speed) or 
% excess time, if continuous flight is possible

if (settings.evaluation.findalt~=1)
        [PerfResults,flightdata] = performanceEvaluator(params, plane, environment, settings);
elseif (settings.evaluation.findalt == 1)   
    PerfResults.t_excess=0;
    hL = 0; %min
    hU = 30000; % max
    environment.h_max = hU;
    h = hL:(hU-hL)/60.0:hU;
     
    for i = 1:numel(h)
        environment.h_0 = h(i);
        [PerfResults, flightdata] = performanceEvaluator(params, plane, environment, settings);
        if (isnan(PerfResults.t_excess) || PerfResults.t_excess<=0)
            break;
        else
            PerfResults.h_max = environment.h_0;
        end
    end
end
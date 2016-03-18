function [ PSolar, irr_vec, etas] = CalculateIncomingSolarPower(t,h,environment,plane,settings,params)
%CALCULATEINCOMINGSOLARPOWER Summary of this function goes here
%   Small wrapper function for calculating the absolute incoming solar power from the incoming irradiance and the plane and
%   environmental properties
    
    % Precalcs
    c_temp = 1.0 - params.solar.k_temp*(environment.T_ground - (273.15 + 25));
    eta_solar = params.solar.eta_sc * c_temp * params.solar.eta_mppt;

    % Irradiation
    irr_vec = solar_radiation_on_surface2(0,0,mod(t+environment.add_solar_timeshift,86400),(environment.dayofyear+floor((t+environment.add_solar_timeshift)/86400)),...
                environment.lat,environment.lon, h,environment.albedo);
    
    % AOI-handling
    epsilonAOI = 1.0;
    epsilonDiff = 1.0;
    if (settings.useAOI)
        AOI = irr_vec(7);
        epsilonAOI = interp1(params.solar.angle_AOI,params.solar.epsilon_AOI, AOI, 'linear', 0.0);
        epsilonDiff = params.solar.epsilon_diff;
    end
    
    % Irradiation levels
    % Could also integrate irradiation levels Isolar here
    
    %P-Solar calculation
    PSolar_global = environment.clearness * irr_vec(1) * plane.ExpPerf.solar.surface * eta_solar * params.solar.eta_cbr * epsilonAOI;
    if(settings.useDirDiffRad)
        PDiffuse = environment.clearness * irr_vec(3) * plane.ExpPerf.solar.surface * eta_solar * epsilonDiff;
        PDirect = environment.clearness * irr_vec(2) * plane.ExpPerf.solar.surface * eta_solar * params.solar.eta_cbr * epsilonAOI;
        PSolar = PDiffuse + PDirect;
    else
        PSolar = PSolar_global;
    end
    
    %Re-calculate eta's
    etas(1) = PSolar / (environment.clearness * irr_vec(1) * plane.ExpPerf.solar.surface);
    etas(2) = epsilonAOI;
    
end
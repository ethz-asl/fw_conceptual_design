function [ PSolar, irr_vec ] = CalculateIncomingSolarPower(t,h,environment,plane,eta_solar)
%CALCULATEINCOMINGSOLARPOWER Summary of this function goes here
%   Small wrapper function for calculating the absolute incoming solar power from the incoming irradiance and the plane and
%   environmental properties
    
    irr_vec = solar_radiation_on_surface2(0,0,mod(t+environment.add_solar_timeshift,86400),(environment.dayofyear+floor(t/86400)),...
                environment.lat,0, h,environment.albedo);
    PSolar = environment.clearness * irr_vec(1) * plane.ExpPerf.solar.surface * eta_solar;
end
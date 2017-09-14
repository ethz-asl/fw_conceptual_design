function [ P_elec_level ] = CalcPFromPolars( plane, params, rho, mu)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    %STEP 1a: Calc Re,Cl,Cd for flight at h_0:
    A_wing = plane.struct.b^2/plane.struct.AR;
    v = sqrt(2*(plane.m_no_bat+plane.bat.m)*params.physics.g/(rho*A_wing*mean(plane.polar.c_L_cn)));
    
    %refine iteratively
    delta_v = 1;
    while abs(delta_v) > 0.01
        Re = min([rho*A_wing/plane.struct.b*v/mu,plane.polar.ReList(length(plane.polar.ReList))]);
        if Re<plane.polar.ReList(1)
            P_elec_level = -1.0;
            display('Error calculating airspeed from polar. Reynolds number too low');
            return;
        end
        v_old = v;
        CL = interp1(log(plane.polar.ReList),plane.polar.c_L_cn,log(Re)); % logarithmic interp.
        v = sqrt((plane.m_no_bat+plane.bat.m)*params.physics.g/(0.5*rho*A_wing*CL)); 
        delta_v=(v-v_old)/v;
    end
    CD = interp1(log(plane.polar.ReList),plane.polar.c_D_cn,log(Re));

    %Step1b: Calc Total electric power for level flight
    P_level = 1/2*rho*A_wing*CD*v^3;
    P_elec_level = P_level / (params.prop.eta_ctrl * params.prop.eta_mot * params.prop.eta_grb * params.prop.eta_plr);
end


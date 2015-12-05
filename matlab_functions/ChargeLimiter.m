function [ P_bat_limited ] = ChargeLimiter(SoC_current, P_bat, parameters, plane)
%CHARGELIMITER This function represents a battery-charge-limiter function,
%which can be based either on theoretical (model-based) or experimental
%data.
   
    if(parameters.bat.chrg_lim_type == 1)
        % Model-based (theoretical) charge limiting curve
        
        % Max admissible charging power in [W] at any SoC
        P_bat_max_overall = parameters.bat.chrg_lim_Prel1 * parameters.bat.e_density * plane.bat.m /3600;
        
        % Use an exponentially-decreasing (exp(-c*x)) fct., where P_bat_limited=Prel2*P_max 
        % at SoC2(=1.0). Find the constant first        
        c =-log(parameters.bat.chrg_lim_Prel2);
        
        if(SoC_current < parameters.bat.chrg_lim_SoC1)
            P_bat_limited = P_bat;
            if(P_bat_limited > P_bat_max_overall) P_bat_limited = P_bat_max_overall; end
        elseif(SoC_current >= 1.0)
            P_bat_limited = 0;
        else
            P_bat_limited = exp (-c*(SoC_current-parameters.bat.chrg_lim_SoC1)/(1.0-parameters.bat.chrg_lim_SoC1))*P_bat_max_overall;
            if(P_bat < P_bat_limited)
                P_bat_limited = P_bat;
            end
        end    
    elseif(parameters.bat.chrg_lim_type == 2)
        % Experimentally-recorded charge limiting curve, measured on AS-P Test Flight Day #23
        SoC = [0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.9 0.89 0.86 0.82 0.78];
        MaxChargingCurrent = [0.6 0.8 1.0 1.5 2.0 3.0 3.5 4.0 5.0 6.0 7.0 8.0 10.0];
        Voltage = [25.36 25.33 25.28 25.24 25.2 25.17 25.15 25.13 25.1 25.08 25.03 25.0 24.9];

        MaxChargingPower = MaxChargingCurrent .* Voltage;

        if(SoC_current < min(SoC))
            P_bat_limited = P_bat;
        elseif(SoC_current >= 1.0)
            P_bat_limited = 0;
        elseif(SoC_current > max(SoC))
            P_bat_limited = MaxChargingPower(1);
        else
            P_bat_limited = interp1(SoC, MaxChargingPower, SoC_current);
            if(P_bat < P_bat_limited)
                P_bat_limited = P_bat;
            end
        end    
    else
        % Don't use charge limiting at all
        P_bat_limited = P_bat;
    end
end


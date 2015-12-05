function [rho,mu]=findAir(h,T_ground)

T_std=288.15;   %ICAO standard temperature at 0m
[rho,a,T,P,mu]=stdatmo(h,T_ground-T_std);

% g           =   9.806;
% h_b         =   0;           % Reference Altitude [m]
% rho_b       =   1.225;       % Air density [kg/m^3] at reference altitude
% T_b         =   T_ground;    % Temperature [°K] at reference altitude
% rho         =   rho_b*exp(-g*0.0289644*(h-h_b)... 
%                 /(8.31432*T_b)); % Air density [kg/m^3]
% if h<11000
%     T       = T_b*(1-22.558e-6*h);
% elseif h< 20000
%     T       = 0.7519*T_b;
% elseif h< 32000
%     T       = 0.7519*T_b + 0.001*(h-20000);
% else
%     T=200;
%     disp('Too high; temperature not modeled');
% end
% 
% % Dynamic air viscosity [N*S/m^2] Sutherland
% mu          =   18.27e-6*(T/T_b)^(1.5)*(T_b+110.4)/(T+110.4);

end
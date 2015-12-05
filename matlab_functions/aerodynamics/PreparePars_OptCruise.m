function [OptCruisePars]=PreparePars_OptCruise(m,rho,mu,A,c_wing,g,polar)
%Performs a precalculation of important cruise parameters. These are
%constant over rho, thus only have to be recalculated if height chance
%occurs -> computational advantages.

% For every Re, get the CL required to lift the A/C
v_array=polar.ReList.*mu./rho./c_wing;
CL_array=2*m*g/rho/A./v_array.^2;

%get respective CD at that CL from polar and calculate power necessary
for i=1:length(polar.ReList)
    %CLtest=2*m*g/rho/A/v(i)^2;
    if(CL_array(i)>polar.c_L_max(i) || CL_array(i)<polar.c_L_min(i)) 
        CD_array(i)=NaN;
    else
        CD_array(i)=interp1(polar.c_L_array,polar.c_D_array{1,i},CL_array(i),'linear');
    end
end
P_array=0.5*rho.*CD_array*A.*(v_array.^3);

% plot(polar.ReList,CL_array);
% hold on
% plot(polar.ReList,CD_array);
% hold on
% plot(polar.ReList,P_array);

%interpolate over P to get Re
min_idx=find(isnan(P_array)==0,1,'first');
max_idx=find(isnan(P_array)==0,1,'last');

OptCruisePars.v_array=v_array;
OptCruisePars.CL_array=CL_array;
OptCruisePars.CD_array=CD_array;
OptCruisePars.P_array=P_array;
OptCruisePars.min_idx=min_idx;
OptCruisePars.max_idx=max_idx;

end
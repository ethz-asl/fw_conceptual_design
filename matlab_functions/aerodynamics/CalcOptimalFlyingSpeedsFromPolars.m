function [v_rmax v_tmax] = CalcOptimalFlyingSpeedsFromPolars()
    % THIS FUNCTION IS UNTESTED AND HAS TO BE ADAPTED TO WORK
    
    %Just a calculation of optimum flying speed (lowest power and best glide
    %ratio)
    
    g = parameters.physics.g;
    if(environment.usemars==1)
        [rho,mu] = findMarsAir(environment.h,environment.T_ground);
    else
        [rho,mu] = findAir(environment.h,environment.T_ground);
    end
    v_tmax=sqrt(2*(m_distr+m_central+m_struct+masses.m_prop)*g/(rho*b^2/AR*mean(polar.c_L_cn)));
    %refine iteratively
    delta_v_tmax=1;
    while abs(delta_v_tmax)>0.01
        Re=min([rho*b/AR*v_tmax/mu,polar.ReList(length(polar.ReList))]);
        v_tmax_old = v_tmax;
        CL = interp1(log(polar.ReList),polar.c_L_cn,log(Re));
        v_tmax = sqrt((m_distr+m_central+m_struct+masses.m_prop)*g/(0.5*rho*b^2/AR*CL)); % logarithmic interp.
        delta_v_tmax=(v_tmax-v_tmax_old)/v_tmax;
    end
    CL = interp1(log(polar.ReList),polar.c_L_gr,log(Re));
    v_rmax=sqrt(2*(m_distr+m_central+m_struct+masses.m_prop)*g/(rho*b^2/AR*CL));
    %refine iteratively
    delta_v_rmax=1;
    while abs(delta_v_rmax)>0.01
        Re=min([rho*b/AR*v_rmax/mu,polar.ReList(length(polar.ReList))]);
        v_rmax_old = v_rmax;
        CL = interp1(log(polar.ReList),polar.c_L_gr,log(Re));
        v_rmax = sqrt((m_distr+m_central+m_struct+masses.m_prop)*g/(0.5*rho*b^2/AR*CL)); % logarithmic interp.
        delta_v_rmax=(v_rmax-v_rmax_old)/v_rmax;
    end
end
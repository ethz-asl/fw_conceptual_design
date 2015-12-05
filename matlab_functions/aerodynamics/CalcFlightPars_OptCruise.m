function [v,Re,CL,CD]=CalcFlightPars_OptCruise(P,Plevel,m,rho,mu,A,c_wing,g,polar,OptCruisePars)

P=(1-0.6*(P-Plevel)/P)*P; 

Re=interp1(OptCruisePars.P_array(OptCruisePars.min_idx:OptCruisePars.max_idx),polar.ReList(OptCruisePars.min_idx:OptCruisePars.max_idx),P,'linear'); %first guess

v=Re*mu/rho/c_wing;
CL=interp1(polar.ReList,OptCruisePars.CL_array,Re);
CD=interp1(polar.ReList(OptCruisePars.min_idx:OptCruisePars.max_idx),OptCruisePars.CD_array(OptCruisePars.min_idx:OptCruisePars.max_idx),Re);

%To take into account decreased efficiency at higher speeds
%Update this to represent real propulsion efficiency over P_prop

% %Calc CD/CL^1.5 (=C1) for optimal cruise
% C1=sqrt(rho*A/2/(m*g)^3)*P;
% 
% %At all reynolds numbers, if CD/CL^1.5==C1, then the A/C can be at level flight
% %consuming the power specified by P. 
% 
% for i=1:length(polar.ReList)
%     temp=polar.c_D_array{1,i};
%     C1vsRe(i,:)=real(temp./(polar.c_L_array.^1.5));
% end
% 
% %we need a positive CL
% min_idx=find(polar.c_L_array(1,:)>0,1,'first');
% for i=1:length(polar.ReList)
%     temp= numel(polar.c_L_array);
%     for j=temp:-1:min_idx
%         if(isnan(C1vsRe(i,j))==0) max_idx(i)=j; break; end
%     end
%     idx_Hit(i)=interp1(C1vsRe(i,min_idx:max_idx(i)),min_idx:max_idx(i),C1);
%     CLvsRe_Hit(i)=interp1(C1vsRe(i,min_idx:max_idx(i)),polar.c_L_array(min_idx:max_idx(i)),C1);
% end
%     
% %Now we need to find v (->Re) such that both FL=Fg and PProp=P!
% for i=1:length(polar.ReList)
%     CL(i)=interp1(min_idx:max_idx(i),polar.c_L_array(min_idx:max_idx(i)),idx_Hit(i),'linear');
%     CD(i)=interp1(min_idx:max_idx(i),polar.c_D_array{1,i}(min_idx:max_idx(i)),idx_Hit(i),'linear');
%     C1check(i)=CD(i)/CL(i)^1.5;
%     v=polar.ReList(i)/rho/c_wing*mu;
%     FL_array(i)= 0.5*rho*CL(i)*A*v^2;
%     P_array(i) = 0.5*rho*CD(i)*A*v^3;
% end
% if(find(isnan(FL_array)==1))
%     display('hups');
% end
% min_idx=find(FL_array>0,1,'first');
% max_idx=find(isnan(FL_array)==0,1,'last');
% 
% FL=m*g;
% if(m*g<min(FL_array(min_idx:max_idx))) FL=min(FL_array(min_idx:max_idx)); end;    %Kindof a bad hack to prevent that we cannot find the FL
% Re1=interp1(FL_array(min_idx:max_idx),polar.ReList(min_idx:max_idx),FL,'linear');
% if(P<min(P_array(min_idx:max_idx))) P=min(P_array(min_idx:max_idx)); end;    %Kindof a bad hack to prevent that we cannot find the P
% Re2=interp1(P_array(min_idx:max_idx),polar.ReList(min_idx:max_idx),P,'linear');
% Re=(Re1+Re2)/2.0;
% v=Re/rho/c_wing*mu;
% %v1=Re1/rho/c_wing*mu;
% %v2=Re2/rho/c_wing*mu;
% 
% %Test:
% CLCheck=interp1(polar.ReList,CL,Re,'linear');
% CDCheck=interp1(polar.ReList,CD,Re,'linear');
% %Check results
% Pcheck=0.5*rho*CDCheck*A*v^3;
% FLcheck=0.5*rho*CLCheck*A*v^2;
% CL=CLCheck;
% CD=CDCheck;

% % % % % % % % % % % % For every Re, get the CL required to lift the A/C
% % % % % % % % % % % v_array=polar.ReList.*mu./rho./c_wing;
% % % % % % % % % % % CL_array=2*m*g/rho/A./v_array.^2;
% % % % % % % % % % % 
% % % % % % % % % % % %get respective CD at that CL from polar and calculate power necessary
% % % % % % % % % % % for i=1:length(polar.ReList)
% % % % % % % % % % %     %CLtest=2*m*g/rho/A/v(i)^2;
% % % % % % % % % % %     if(CL_array(i)>polar.c_L_max(i) || CL_array(i)<polar.c_L_min(i)) 
% % % % % % % % % % %         CD_array(i)=NaN;
% % % % % % % % % % %     else
% % % % % % % % % % %         CD_array(i)=interp1(polar.c_L_array,polar.c_D_array{1,i},CL_array(i),'linear');
% % % % % % % % % % %     end
% % % % % % % % % % % end
% % % % % % % % % % % P_array=0.5*rho.*CD_array*A.*(v_array.^3);
% % % % % % % % % % % 
% % % % % % % % % % % % plot(polar.ReList,CL_array);
% % % % % % % % % % % % hold on
% % % % % % % % % % % % plot(polar.ReList,CD_array);
% % % % % % % % % % % % hold on
% % % % % % % % % % % % plot(polar.ReList,P_array);
% % % % % % % % % % % 
% % % % % % % % % % % %interpolate over P to get Re
% % % % % % % % % % % min_idx=find(isnan(P_array)==0,1,'first');
% % % % % % % % % % % max_idx=find(isnan(P_array)==0,1,'last');


% if(isnan(Re))
%     %Special case, all powers are too high. Generate an additional point
%     %, which is the minimum Re at which the required FL can still be
%     %generated. Use this for interpolation then.
%     
%     %we assume linearity of CLMax between the two points of interest. Then
%     %the minimum Re, for which FL can be reached, can be expressed through
%     %a third-order polynomial
%     idx2=find(isnan(P_array)==0,1,'first');
%     idx1=idx2-1;
%     CLMax1=polar.c_L_max(idx1);CLMax2=polar.c_L_max(idx2);
%     Re1=polar.ReList(idx1); Re2=polar.ReList(idx2);
%     Re_array=Re1:1000:Re2;
%     CLReq=m*g*2*rho/A*(c_wing./Re_array/mu).^2;
%     CLMax=CLMax1+(Re_array-Re1)./(Re2-Re1)*(CLMax2-CLMax1);
%     for i=1:length(Re_array)
%         dCL(i)=abs(CLReq(i)-CLMax(i));
%     end
%     Retmpidx=find(dCL==min(dCL));
%     Retmp=Re_array(Retmpidx);
%     vtmp=Retmp*mu/c_wing/rho;
%     %Now check
%     FLMax=0.5*rho*CLMax(Retmpidx)*A*vtmp^2;
%     if(FLMax<m*g) Retmpidx=Retmpidx+1; end;
%     %Now calc CD & P at that point
%     CD1=interp1(polar.c_L_array, polar.c_D_array{1,idx1},CLMax1);
%     CD2=interp1(polar.c_L_array, polar.c_D_array{1,idx2},CLMax2);
%     CDtmp=interp1([CLMax1 CLMax2],[CD1 CD2],CLMax(Retmpidx));
%     Ptmp=0.5*rho*CDtmp*A*vtmp^3;
%     
%     %ReGuess Re
%     Re=interp1([Ptmp P_array(min_idx)],[Retmp polar.ReList(min_idx)],P,'linear');
%     
% %     p(4)=-m*g;
% %     p(3)=0;
% %     p(2)=A/2/rho*mu^2/c_wing^2*(CLMax1-Re1/(Re2-Re1)*(CLMax2-CLMax1));
% %     p(1)=A/2/rho*mu^2/c_wing^2/(Re2-Re1);
% %     
% %     ReSolutions=roots(p);
% %     for i=1:length(ReSolutions)
% %         if(imag(ReSolutions(i))==0 && ReSolutions(i)>=0) ReMin=ReSolutions(i); end
% %     end
% %     %Check: 
% %     CLAvailAt_ReMin=interp1(polar.ReList,polar.c_L_max,ReMin);
% %     v=ReMin*mu/rho/c_wing;
% %     CLRequAt_ReMin=2*m*g/rho/A/v^2;
% end
    
%Check
% F_Check=0.5*rho*CL*A*v^2;
% PCheck=0.5*rho*CD*A*v^3;
% dF=F_Check-m*g;
% dP=PCheck-P;
end

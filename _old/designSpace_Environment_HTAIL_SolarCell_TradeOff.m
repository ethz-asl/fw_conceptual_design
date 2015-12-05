%==========================================================================
% Analyze the Design Space of the senseSoar Airplane, but with respect
% to the environmental conditions actually prevalent! This function thus
% has to be called after the optimization of the Airplane (which has been
% done through the designSpaceNew-function)
%
% Philipp Oettershagen
% 05/2012
%==========================================================================

% TODO!!!:
% - Set VMax at level flight if max. altitude reached (->use all power!)
% - Multiday flight simulation verification


% Initialize
clear variables;
close all;
clc;
% use double core- processing
 %matlabpool open 2
 %options = optimset('UseParallel','always');
initParameters;

% Editable Section:
% =================

% Set fixed aircraft parameters
%**************************************************************************
% MOD for different AC configs
%**************************************************************************
%config1: BASECONFIG. b=5.4m, 100 cells (16@tail)
b=5.6;
AR=18.0;
delta_mass=0.0;
nrcells=100;

%config2: b=5.4m, 84cells (0@tail)
% b=5.4;
% AR=17.7;
% delta_mass=-0.06; %cable & HTP weight, solar cells deducted automatically, wing weight constant
% nrcells=84;

%config3: b=5.656m, 88cells (0@tail)
% b=5.656;
% AR=18.54;
% delta_mass=-0.06; %cable & HTP weight, 
% nrcells=88;

%config4: b=5.656m, 104cells (16@tail)
% b=5.656;
% AR=18.54;
% delta_mass=0; 
% nrcells=104;

%calculate new AC-config parameters
parameters.avionics.mass      =  0.15+delta_mass; %Set delta mass here
parameters.solar.rWngCvrg     =  nrcells*0.125^2/b/(b/AR);
% ATTENTION: ADD CL^1.5/CD changes manually!!!

% Set environment and aircraft parameters to be variated
clearness_min=0.30;
clearness_max=1.0;
clearness_step=0.05;
turbulence_min=0.0;
turbulence_max=0.7;
turbulence_step=0.05;
m_bat_min  = 2.0;
m_bat_max  = 2.0;
m_bat_step = 0.2;

% Set environment, payload and airfoil
environment.month = 6;
environment.day = 21;
environment.dayofyear = 5*30.5+21;
environment.h = 100;
environment.hmax = 5500*0.3048;        %Maximum altitude not be exceeded [m]
environment.T_ground = 300.0; %288.15;
environment.lat = 38.0;
environment.albedo = 0.1;       %0.1 for atlantic, more if above ground
environment.clearness = NaN;    %Clearness or Cloud-Coverage-factor(CCF)
environment.turbulence = NaN;   %Turbulence, defined here as Plevel_turb/Plevel_noturb=(1+turbulence)
environment.usemars = 0;
payload.mass = 0.18;
payload.power = 3;

% Change some parameters
parameters.multidaysim        =  1;        % if 0, simulates a single day (-> pure excess time).
                                           % if 1, simulates multiday flight
parameters.propulsion.number  =  1;        % Number of propulsion units [-]
parameters.structure.shell    =  0;        % 1 for shell wing, 0 for rib wing
parameters.evaluation.clmb    =  1;        % 1 to allow altitude changes
parameters.optGRcruise        =  0;        % 1 to allow cruise at optimal glide ratio & speed when max altitude reached 
parameters.evaluation.findalt =  0;        % if 1, it finds the maximum
                                           % altitude for eternal flight
parameters.dt                 =  100;      % Discretization time interval [s]
                                        
%--------------------------------------------------------------

% Evaluation
% ==========

% Number of configurations calculated:
N = ((clearness_max-clearness_min)/clearness_step+1)...
    *((turbulence_max-turbulence_min)/turbulence_step+1)...
    *((m_bat_max-m_bat_min)/m_bat_step+1);
%disp(['Number of configurations to be calculated: ' num2str(N)])
%disp(['Expected processing time: ' num2str(N*10/3600) ' h'])

% Start the evaluation
h=waitbar(0,'Progress');
i=0;
tic
l=0;

%data storage matrices (3D)
clearness_array=clearness_min:clearness_step:clearness_max;
turbulence_array=turbulence_min:turbulence_step:turbulence_max;
m_bat_array=m_bat_min:m_bat_step:m_bat_max;
t_excess=zeros(length(turbulence_array),length(clearness_array),length(m_bat_array));
t_endurance=zeros(length(turbulence_array),length(clearness_array),length(m_bat_array));
t_chargemargin=zeros(length(turbulence_array),length(clearness_array),length(m_bat_array));
%continuousflight=zeros(length(turbulence_array),length(clearness_array),length(m_bat_array));
v_tmax=zeros(length(turbulence_array),length(clearness_array),length(m_bat_array));
m_struct=zeros(length(turbulence_array),length(clearness_array),length(m_bat_array));

for m_bat=m_bat_array
    l=l+1;
    m=0;
    for turbulence=turbulence_array
        environment.turbulence=turbulence; %Just a hack cause matlab doesnt allow to put environment.turbulence @for loop
        m=m+1;
        n=0;
        for clearness=clearness_array
            environment.clearness=clearness; 
            n=n+1;
            
            [performance,tmpflightdata,polar,masses] = ...
               evaluateSolution(b,AR,m_bat,payload,environment,parameters);
            % store
            t_excess(m,n,l)=performance.t_excess;
            t_endurance(m,n,l)=performance.t_endurance;
            t_chargemargin(m,n,l)=performance.t_chargemargin;
            %continuousflight(m,n,l)=performance.continuousflight;
            v_tmax(m,n,l)=performance.v_tmax;
            m_struct(m,n,l)=masses.m_struct;
            if(isempty(tmpflightdata)==0) 
                flightdata_array(m,n,l)=tmpflightdata; 
            end
            flightdata_array(m,n,l).AR=AR;
            flightdata_array(m,n,l).b=b;
            flightdata_array(m,n,l).m_bat=m_bat;
            
            % waiting times
            i=i+1;
            elapsed = toc;
            remaining = N*elapsed/i-elapsed;
            rhours=floor(remaining/3600);
            rminutes=floor((remaining-rhours*3600)/60);
            rseconds=floor((remaining-rhours*3600-rminutes*60));
            %disp(['Time remaining: ' num2str(rhours,'%02.0f') ':' ...
            %    num2str(rminutes,'%02.0f') ':' num2str(rseconds,'%02.0f')]);
            waitbar(i/N,h,['Time remaining: ' num2str(rhours,'%02.0f') ':' ...
                num2str(rminutes,'%02.0f') ':' num2str(rseconds,'%02.0f')]);
        end
    end
end
close(h)

%Plotting
l=1;
i=0;
for m_bat=m_bat_min:m_bat_step:m_bat_max
    if(length(clearness_array)<2 ||length(turbulence_array)<2) continue; end
    
    i=i+1;
    %h12(l)=figure(l+1); %TOCHANGE
    
    %subplot(2,2,1)
    [c2,hc2]=contourf(clearness_min:clearness_step:clearness_max,...
        turbulence_min:turbulence_step:turbulence_max, t_excess(:,:,l),...
        200,'Linestyle','none','ShowText','off');
    xlabel('CCF [-]')
    ylabel('Turbulence [-]');
    title(['Excess Time [h] (m_bat=' num2str(m_bat) 'kg)']);
    %caxis([-15,15])
    caxis([0,max(max(max( t_excess(:,:,:))))])
    colorbar
    
    nc = get(hc2,'Children');
    temp = 100;
    select_i=zeros(length(nc),1);
    for i = 1:length(nc)
       ud1 = get(nc(i),'UserData');   
       if (abs(ud1) < temp)
           temp = abs(ud1);
       end
    end
    filler=1;
    for i = 1:length(nc)
       ud1 = get(nc(i),'UserData');   
       if (abs(ud1) == temp)
           select_i(filler)=i;
           filler=filler+1;
       end
    end
    j=1;
    if i>0
        while(select_i(j)~=0)
            set(nc(select_i(j)),'Linestyle','-');
            set(nc(select_i(j)),'LineWidth',2);
            j=j+1;
        end
    end

    
%     subplot(2,2,2)
    figure
    contourf(clearness_min:clearness_step:clearness_max,...
        turbulence_min:turbulence_step:turbulence_max, t_endurance(:,:,l),...
        200,'linestyle','none');
    xlabel('ccf [-]')
    ylabel('turbulence [-]');
    title(['endurance [h] (m_bat=' num2str(m_bat) 'kg)']);
    caxis([0,48])
    colorbar
    
%     subplot(2,2,3)
%     contourf(b_min:b_step:b_max,...
%         m_bat_min:m_bat_step:m_bat_max, v_tmax(:,:,l)',...
%         200,'Linestyle','none');
%     xlabel('Wingspan (b) [m]')
%     ylabel('Battery mass (m\_bat) [kg]');    
%     title(['Nominal Speed [m/s] (AR=' num2str(AR) ')']);
%     caxis([min(min( v_tmax(:,:,l))),max(max( v_tmax(:,:,l)))])
%     colorbar

%     subplot(2,2,3)
    figure    
    contourf(clearness_min:clearness_step:clearness_max,...
        turbulence_min:turbulence_step:turbulence_max, t_chargemargin(:,:,l),...
        200,'Linestyle','none');
    xlabel('CCF [-]')
    ylabel('Turbulence [-]');   
    title(['Charge margin [h] (m_bat=' num2str(m_bat) 'kg)']);
    caxis([0,max(max(max( t_chargemargin(:,:,:))))])
    colorbar
    
%     subplot(2,2,4)
%     contourf(clearness_min:clearness_step:clearness_max,...
%         turbulence_min:turbulence_step:turbulence_max, t_endurance(:,:,l),...
%         200,'Linestyle','none');
%     xlabel('CCF [-]')
%     ylabel('Turbulence [-]');
%     title(['Structural Mass [kg] (m_bat=' num2str(m_bat) 'kg)']);
%     caxis([min(min( m_struct(:,:,l))),max(max(m_struct(:,:,l)))])
%     colorbar
    
    l=l+1;
end

%DetailedAnalysisForm({AR_array,b_array,m_bat_array,flightdata_array,t_excess,t_endurance,t_chargemargin});
%end use double core- processing
 %matlabpool close
 %print(h12(1),'test2.emf','-dmeta','-painters')
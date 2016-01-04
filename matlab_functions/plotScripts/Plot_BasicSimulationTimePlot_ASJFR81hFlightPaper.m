function [] = Plot_BasicSimulationTimePlot_ASJFR81hFlightPaper(flightdata,environment,params, plane)
% Basic plot of flight data vs. time
%   Plot_BasicSimulationTimePlot(...) plots the flightdata (altitude, batters
%   state, all power input and output components) over time.

    time = flightdata.t_array ./3600;

    %str=strcat('24h flight overview',' b=',num2str(flightdata.b),' AR=',num2str(flightdata.AR),' mbat=',num2str(flightdata.m_bat));%: DayOfYear=',num2str(environment.dayofyear),' CLR=',num2str(environment.clearness),' Turb=',num2str(environment.turbulence),' b=',num2str(plane.struct.b),' AR=',num2str(plane.struct.AR),' mbat=',num2str(plane.bat.m));
    %figure('Name',str);
    figure
    
    [ax,hLine1,hLine2]=plotyy([time',time',time'],...
            [flightdata.P_solar_array',flightdata.P_elec_tot_array',flightdata.P_charge_array'],...
           time',...
            flightdata.bat_array'/3600);
    l1=legend('P_{solar} [W]','P_{out} [W]','P_{bat} [W]','E_{bat} [Wh]');
    xlabel('Time [h]');
    ylabel(ax(1),'Power [W]');
    ylabel(ax(2),'Energy [Wh]');
    ylim(ax(1),[-60.0 max(flightdata.P_solar_array)*1.5]);
    ylim(ax(2),[-104.0 max(flightdata.bat_array/3600)*1.05]);
    
    hLine2.LineWidth=2;
    hLine2.Color=[1 0 0];
    ax(2).YColor=[1 0 0];
    hLine1(1).Color=[0 0 1];
    hLine1(2).Color=[0 1 0];
    set(ax,'FontSize',24);
        
    %set(gcf, 'Position', [100 100 150 150]);
end


function Plot_AnalysisSpace_DifferentPlaces()
    %**************************************************************************
    %*** PLOTTING
    %**************************************************************************

    today=now;
    doy = today - datenum(2015,1,1) + 1

    % PLOT1: T_EXCESS
    figure('Name','tExc')
    for k=1:numel(powers)
        subplot(numel(powers),1,k);
        for i=1:numel(lat)
            plot(days, squeeze(res.t_excess(i,k,:)),'LineWidth',2);
            hold all
            str{i}=strcat(place{i}, '@(h_0=', num2str(h_0(i)), 'm/m_{bat}=',num2str(plane.bat.m),'kg), P_{prop, level}=',num2str(powers(k)),'W');
        end
        legend(str{1:end});

        ylim([0 max(max(res.t_excess(:,k,:)))+0.5]);
        t_exc_lim=0.15*plane.bat.m*plane.bat.e_density/(plane.P_prop_level+plane.P_avi);
        line([min(days) max(days)],[t_exc_lim t_exc_lim],'Color','red');
        text((min(days)+max(days))/2,t_exc_lim-0.03,'Critical excess time','FontSize',12,'Color','red')
        xlabel('day of year [-]');
        ylabel('Modified excess time [h]');

        line([173 173],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(173,0.5,'  June 21st','FontSize',12)
        line([204 204],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(204,0.5,'  July 21st','FontSize',12)
        line([234 234],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(234,0.5,'  August 21st','FontSize',12)
        line([265 265],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(265,0.5,'  Sept. 21st','FontSize',12)
        %line([doy doy],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color',[1 0.5 0.5]);
        %text(doy,0.5,'  TODAY','FontSize',12)
    end

    % PLOT2: MINIMUM STATE OF CHARGE
    figure('Name','SoC')
    for k=1:numel(powers)
        subplot(numel(powers),1,k);
        for i=1:numel(lat)
            plot(days, squeeze(minSoC(i,k,:)),'LineWidth',2);
            hold all
            str{i}=strcat(place{i}, '@(h_0=', num2str(h_0(i)), 'm/m_{bat}=',num2str(plane.bat.m),'kg), P_{prop, level}=',num2str(powers(k)),'W');
        end
        legend(str{1:end});

        ylim([0 max(max(minSoC(:,k,:)))+0.1]);
        line([min(days) max(days)],[0.10 0.10],'Color','red');
        text((min(days)+max(days))/2,0.10-0.02,'Critical SoC','FontSize',12,'Color','red')
        xlabel('day of year [-]');
        ylabel('Minimum State of Charge [%]');

        line([173 173],[0 max(max(minSoC(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(173,0.3,'  June 21st','FontSize',12)
        line([204 204],[0 max(max(minSoC(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(204,0.3,'  July 21st','FontSize',12)
        line([234 234],[0 max(max(minSoC(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(234,0.3,'  August 21st','FontSize',12)
        line([265 265],[0 max(max(minSoC(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(265,0.3,'  Sept. 21st','FontSize',12)
        %line([doy doy],[0 max(max(minSoC(:,k,:)))+0.5],'LineStyle','--','Color',[1 0.5 0.5]);
        %text(doy,0.3,'  TODAY','FontSize',12)
    end
end
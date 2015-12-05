function Plot_AnalysisSpace_FSRPaper()
    %**************************************************************************
    %*** PLOTTING
    %**************************************************************************

    min_end = 8.0*ones(numel(days),1);
    month = days/30.5;
    cd = round(numel(days)/2); %center day, just middle of the plot

    %PLOT1: T_ENDURANCE
    figure

    for k = 1:numel(powers)
        daylength = (res.t_sunset(i,k,:)-res.t_sunrise(i,k,:))/3600;

        subplot(numel(powers),1,k);

        %plot(month, squeeze(daylength),'LineWidth',2,'Color','red');
        % Green background
        clf
        area(month, squeeze(res.t_endurance(numel(clearness),k,:)) ,'LineWidth',1,'FaceColor',[0.4 1.0 0.4],'EdgeColor',[0.4 1.0 0.4]);
        hold all
        % DayLength
        area(month, squeeze(daylength) ,'LineWidth',2,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
        plot(month, squeeze(daylength) ,'LineWidth',2,'Color','r');

        %h = text(month(cd)-0.5,squeeze(daylength(cd))+0.2,' Length of Daytime'); 
        %set(h,'BackgroundColor','w');
        dayshift = 95;    
        %Plot MinEndurance
        plot(month, min_end,'LineWidth',2,'Color','k');
        h = text(month(cd+dayshift),min_end(1),[' CLR=' num2str(0.0) ' (min. endurance)']); 
        set(h,'BackgroundColor',[1 0.7 0.7]);

        %Plot first clearness values
        plot(month, squeeze(res.t_endurance(1,k,:)),'--','Color',[0.4 0.4 0.4]);
        h = text(month(cd+dayshift),res.t_endurance(1,k,cd+dayshift)-0.8,[' CLR=' num2str(clearness(1))]); 
        set(h,'BackgroundColor',[1 0.7 0.7]);

        %legend('Flight Endurance Min/Max [h]','Flight Endurance [h]','Length of Daytime [h]');

        for i = 2:numel(clearness)-1
            plot(month, squeeze(res.t_endurance(i,k,:)),'--','Color',[0.4 0.4 0.4]);
            h = text(month(cd+dayshift),res.t_endurance(i,k,cd+dayshift)-1.0,[' CLR = ' num2str(clearness(i))]);
            set(h,'BackgroundColor',[0.4 1.0 0.4]);
        end
        plot(month, squeeze(res.t_endurance(numel(clearness),k,:)),'LineWidth',2,'Color','k');
        h = text(month(cd+dayshift),res.t_endurance(numel(clearness),k,cd+dayshift)+0.3,[' CLR=' num2str(clearness(end)) ' (max. endurance)']); 
        set(h,'BackgroundColor','w');

        ylim([7.6, 24]);
        xlim([0,12]);
        xlabel('Month of year [-]');
        ylabel('Flight Endurance [h] vs. atmosphere clearness (CLR)');
        months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
        set(gca,'XTick',1:12);
        set(gca,'XTickLabel',months(1:12));

        line([173/30.5 173/30.5],[0 max(max(res.t_endurance(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(173/30.5-0.4,23.0,'  June 21st')
        line([356/30.5 356/30.5],[0 16],'LineStyle','--','Color','black');
        text(356/30.5-0.8,17.0,'  Dec. 21st')

        set(gcf,'Color','w');
        set(findall(gcf,'-property','FontSize'),'FontSize',12)

        %Fake the legend (some bug here)
        figure
        plot(1,1,'LineWidth',2,'Color','k');
        hold all
        plot(1,2,'--','Color',[0.4 0.4 0.4]); 
        plot(1,2,'LineWidth',2,'Color','r');
        legend('Flight endurance [h] (max/min)','Flight endurance [h]','Length of daylight [h]');
        set(findall(gcf,'-property','FontSize'),'FontSize',12)

    end

    %PLOT1: T_ENDURANCE
    figure(1)
    for k = 1:numel(powers)
        subplot(numel(powers),1,k);
        plot(days, squeeze(res.t_excess(1,k,:)),'LineWidth',2);
        %hold all
        %plot(days, squeeze(res.t_excess(2,k,:)),'LineWidth',2);
        ylim([0 max(max(res.t_excess(:,k,:)))+0.5]);
        t_exc_lim = 0.15*plane.bat.m*plane.bat.e_density/(plane.P_prop_level+plane.P_avi);
        line([min(days) max(days)],[t_exc_lim t_exc_lim],'Color','red');
        text((min(days)+max(days))/2,t_exc_lim-0.03,'Critical excess time','FontSize',12,'Color','red')
        xlabel('day of year [-]');
        ylabel('Modified excess time [h]');

        %str1 = strcat('Barcelona, h_0 = ', num2str(h_0(2)), 'm , P_{prop, level}=',num2str(powers(k)),'W');
        str2 = strcat('Tuggen/CH, h_0=', num2str(h_0(1)), 'm , P_{prop, level}=',num2str(powers(k)),'W');
        legend(str2);

        line([173 173],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(173,0.5,'  June 21st','FontSize',12)
        line([204 204],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(204,0.5,'  July 21st','FontSize',12)
        line([215 215],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(215,0.5,'  August 2nd','FontSize',12)
    end

    %PLOT1: T_DAYLIGHT
    figure(1)
    for k = 1:numel(powers)
        subplot(numel(powers),1,k);
        plot(days, squeeze(res.t_excess(1,k,:)),'LineWidth',2);
        %hold all
        %plot(days, squeeze(res.t_excess(2,k,:)),'LineWidth',2);
        ylim([0 max(max(res.t_excess(:,k,:)))+0.5]);
        t_exc_lim = 0.15*plane.bat.m*plane.bat.e_density/(plane.P_prop_level+plane.P_avi);
        line([min(days) max(days)],[t_exc_lim t_exc_lim],'Color','red');
        text((min(days)+max(days))/2,t_exc_lim-0.03,'Critical excess time','FontSize',12,'Color','red')
        xlabel('day of year [-]');
        ylabel('Modified excess time [h]');

        %str1 = strcat('Barcelona, h_0=', num2str(h_0(2)), 'm , P_{prop, level}=',num2str(powers(k)),'W');
        str2 = strcat('Tuggen/CH, h_0=', num2str(h_0(1)), 'm , P_{prop, level}=',num2str(powers(k)),'W');
        legend(str2);

        line([173 173],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(173,0.5,'  June 21st','FontSize',12)
        line([204 204],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(204,0.5,'  July 21st','FontSize',12)
        line([215 215],[0 max(max(res.t_excess(:,k,:)))+0.5],'LineStyle','--','Color','black');
        text(215,0.5,'  August 2nd','FontSize',12)
    end
end
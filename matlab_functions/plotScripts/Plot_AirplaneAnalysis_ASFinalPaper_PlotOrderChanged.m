function Plot_AirplaneDesign_ASFinalPaper_PlotOrderChanged(PerfResults, DesignResults, environment, plane, params, flightdata, vars, plotnr)
    
    close all;
    
    %--------------------------------------------------------------------------
    %--- PLOT CONFIGURATION
    %--------------------------------------------------------------------------    
    % IMPORTANT: Set type of plot here (or get it via parameterlist of function)
    % 1 = standard plot axes, 2 = Pin/Pout plot axes
    choose_manual_labels = false;
    
    % Plot setup
    FontSize = 16;
    vmargin = 0.20;
    hmargin = 0.045;
    datevector = zeros(1,numel(vars(1).values));
    ylabel_delta=0.0;
    if(plotnr==1)
        datevector = datenum(2015,1*ones(1,numel(vars(1).values)),1*ones(1,numel(vars(1).values)))-2;
        vars(1).values = floor(vars(1).values);
    elseif(plotnr==2) 
        ylabel_delta=1.0;
        hmargin=0.056;
    elseif(plotnr==3)
        datevector = datenum(2015,1*ones(1,numel(vars(3).values)),1*ones(1,numel(vars(3).values)))-2;
        vars(3).values = floor(vars(3).values);
    end
    
    %set(0,'defaulttextinterpreter','latex');
    %--------------------------------------------------------------------------
    %--- OVERVIEW PLOTS
    %--------------------------------------------------------------------------
    if(plotnr<3)
        if(numel(vars(2).values)<2 ||numel(vars(1).values)<2) 
            display('ERROR: Not enough design variables specified, cannot plot results! Please specific at least the first two design variables.');
        else
            for i = 1:numel(vars(3).values)
                %Configure figure
                str = strcat(vars(3).shortname, '=', num2str(vars(3).values(i)));
                figure('Name',str,'OuterPosition',[120 0 1200 900]);
                colormap(flipud(colormap(hot))); %colormap(jet);
                set(gcf, 'Color', 'w');

                %Convert data first
                for k = 1:numel(vars(2).values)
                    for j = 1:numel(vars(1).values)
                        %display([i,' ',k,' ',j]);
                        temp_t_exc(k,j) = PerfResults(i,k,j).t_excess;
                        temp_min_SoC(k,j) = PerfResults(i,k,j).min_SoC;
                        temp_chargemargin(k,j) = PerfResults(i,k,j).t_chargemargin;
                        temp_endurance(k,j) = PerfResults(i,k,j).t_endurance;
                        temp_m_total(k,j) = 0.0;
                        if(~isempty(DesignResults)) DesignResults(i,k,j).m_no_bat + DesignResults(i,k,j).m_bat; end;
                        if(isnan(temp_chargemargin(k,j))) temp_chargemargin(k,j) = 0.0; end %Remove NaNs to have smooth contour plotting below
                    end
                end

                %NaN-check to avoid plotting errors below
                if(sum(~isnan(temp_endurance))==0) temp_endurance(:,:) = 0.0; end
                %if(sum(~isnan(temp_chargemargin))==0) temp_chargemargin(:,:) = 0.0; end

                %PLOT1 : Excess time. Draws contour surface, and contour lines
                ax(1) = subplot_tight(1,3,1,[vmargin hmargin]);
                grid on; set(ax(end),'GridLineStyle',':');
                hold on
                [c2,hc2] = contourf(vars(1).values+datevector,vars(2).values+ylabel_delta,temp_t_exc,[0.0:0.05:min(max(max(max(temp_t_exc))),24.0)],'Linestyle','none');
                %[c2,hc2] = contour(vars(1).values+datevector,vars(2).values+ylabel_delta,temp_t_exc,[0.05 0.05],'LineColor',[0 0 0]);
                [c2,hc2] = contour(vars(1).values+datevector,vars(2).values+ylabel_delta,temp_t_exc,[0:1:min(max(max(max(temp_t_exc))),24.0)],'LineColor',[0 0 0]);
                if(choose_manual_labels) clabel(c2,hc2,'manual'); 
                else clabel(c2,hc2,[0:1:min(max(max(max(temp_t_exc))),24.0)],'LabelSpacing',144,'FontSize',FontSize);
                end
                title('Excess Time [h]','FontSize',FontSize);
                if(plotnr == 2) ylabel(ax(i),'P_{out}/P_{out,nom} [-]','FontSize',FontSize);
                else h = ylabel(vars(2).name,'FontSize',FontSize);
                end
                h = xlabel(vars(1).name,'FontSize',FontSize);
                cbar_handle(1) = colorbar('FontSize',FontSize-2,'Location','northoutside');
                caxis([0,min(max(max(max(temp_t_exc))),24.0)]);

                % PLOT2 : Charge margin
                ax(end+1) = subplot_tight(1,3,2,[vmargin hmargin]);
                grid on; set(ax(end),'GridLineStyle',':');
                hold on
                [c2,hc2] = contourf(vars(1).values+datevector,vars(2).values+ylabel_delta,temp_chargemargin, [0.0:0.05:min(max(max(max(temp_chargemargin))),24.0)],'Linestyle','none');
                [c2,hc2] = contour(vars(1).values+datevector,vars(2).values+ylabel_delta,temp_chargemargin,'LevelStep',1.0,'LineColor',[0 0 0]);
                if(choose_manual_labels) clabel(c2,hc2,'manual'); 
                else clabel(c2,hc2,[0:1:min(max(max(max(temp_chargemargin))),24.0)],'LabelSpacing',144,'FontSize',FontSize);
                end
                title('Charge margin [h]','FontSize',FontSize);
                %ylabel(vars(2).name);
                xlabel(vars(1).name,'FontSize',FontSize);
                cbar_handle(end+1) = colorbar('FontSize',FontSize-2,'Location','northoutside');
                caxis([0,min(max(max(max(temp_chargemargin))),24.0)]);            

                % PLOT3 : Endurance
                ax(end+1) = subplot_tight(1,3,3,[vmargin hmargin]);
                grid on; set(ax(end),'GridLineStyle',':');
                hold on
                [c2,hc2] = contourf(vars(1).values+datevector,vars(2).values+ylabel_delta,temp_endurance, [0.0:0.1:max(max(max(temp_endurance)))],'Linestyle','none');
                [c2,hc2] = contour(vars(1).values+datevector,vars(2).values+ylabel_delta,temp_endurance,'LevelStep',2.0,'LineColor',[0 0 0]);
                if(choose_manual_labels) clabel(c2,hc2,'manual'); 
                else clabel(c2,hc2,[0:2:max(max(max(temp_endurance)))],'LabelSpacing',144,'FontSize',FontSize);
                end
                title('Endurance [h]','FontSize',FontSize);
                %ylabel(vars(2).name);
                xlabel(vars(1).name,'FontSize',FontSize);
                cbar_handle(end+1) = colorbar('FontSize',FontSize-2,'Location','northoutside');
                caxis([0,max(max(max(temp_endurance(:,:))))]);

                % PLOT4: Optional plots can be chosen here

                % FINAL CONFIG
                for i = 1:numel(ax)
                    set(ax(i),'FontSize',FontSize-2);
                    if(plotnr<2)
                        PlotEveryXthTick = 10;
                        set(ax(i),'XTick',floor(vars(1).values(1:PlotEveryXthTick:end)+datevector(1:PlotEveryXthTick:end)));
                        set(ax(i),'XTickLabel',datestr(floor(vars(1).values(1:PlotEveryXthTick:end)+datevector(1:PlotEveryXthTick:end)),'mmmdd'));
                        xlim([vars(1).values(1)+datevector(1) vars(1).values(end)+datevector(end)]);
                        ax(i).XTickLabelRotation=45;
                    elseif(plotnr==2)
                        xlabel(ax(i),'P_{in}/P_{in,nom} [-]','FontSize',FontSize);
                    end
                    set(ax(i),'layer','top')
                end
            end
        end
    elseif(plotnr==3)
        %Configure figure
        figure('Name','Payload plot','OuterPosition',[120 0 1200 680]);
        set(gcf, 'Color', 'w');
        FontSize = 14;
        hmargin = 0.055;

        %Convert data first
        for i = 1:numel(vars(3).values)
            for k = 1:numel(vars(2).values)
                for j = 1:numel(vars(1).values)
                    %display([i,' ',k,' ',j]);
                    temp_t_exc(i,k,j) = PerfResults(i,k,j).t_excess;
                    temp_min_SoC(i,k,j) = PerfResults(i,k,j).min_SoC;
                    temp_chargemargin(i,k,j) = PerfResults(i,k,j).t_chargemargin;
                    temp_endurance(i,k,j) = PerfResults(i,k,j).t_endurance;
                    if(temp_t_exc(i,k,j)<0.001)
                        temp_t_exc(i,k,j) = NaN;
                        temp_min_SoC(i,k,j) = NaN;
                        temp_chargemargin(i,k,j) = NaN;
                    end
                    temp_m_total(i,k,j) = 0.0;
                    if(~isempty(DesignResults)) DesignResults(i,k,j).m_no_bat + DesignResults(i,k,j).m_bat; end;
                    %if(isnan(temp_chargemargin(i,k,j))) temp_chargemargin(i,k,j) = 0.0; end %Remove NaNs to have smooth contour plotting below
                end
            end
        end

        StyleArr = {{'-b'};{'-g'};{':g'};{'-c'};{':c'};{'-r'};{':r'};{'--m'}};

        %PLOT1 : Excess time.
        ax(1) = subplot_tight(1,2,1,[vmargin hmargin]);
        grid on; set(ax(end),'GridLineStyle',':');
        hold on
        for k=1:numel(vars(2).values)
                plot(vars(3).values+datevector, temp_t_exc(:,k,k),StyleArr{k}{1,1},'LineWidth',2); %,CLRArr(i),'LineStyle',StyleArr(i)
                hold all
                str{k}=strcat('m_{pld}=',num2str(vars(1).values(k)),'kg, P_{pld}=',num2str(vars(2).values(k)),'W'); %(k-1)*numel(vars(1).values)+j
        end
        str{numel(vars(2).values)} = strcat(str{k},' (m_{bat}=3.8kg)');

        %legend(str);
        t_exc_lim=0.10*plane.bat.m*params.bat.e_density/3600/(plane.ExpPerf.P_prop_level+plane.avionics.power+plane.payload.power);
        line([min(vars(3).values+datevector) max(vars(3).values+datevector)],[t_exc_lim t_exc_lim],'Color','red');
        text((min(vars(3).values+datevector)+max(vars(3).values+datevector))/2.0,t_exc_lim+0.32,'Critical excess time (SoC=10%)','FontSize',12,'Color','red','BackgroundColor','w')
        xlabel('Day of year [-]','FontSize',FontSize);
        ylabel('Excess time T_{exc} [h]','FontSize',FontSize);

        %PLOT2 : Charge margin
        ax(end+1) = subplot_tight(1,2,2,[vmargin hmargin]);
        grid on; set(ax(end),'GridLineStyle',':');
        hold on
        for k=1:numel(vars(2).values)
                plot(vars(3).values+datevector, temp_chargemargin(:,k,k),StyleArr{k}{1,1},'LineWidth',2);
                hold all
                str{k}=strcat('m_{pld}=',num2str(vars(1).values(k)),'kg, P_{pld}=',num2str(vars(2).values(k)),'W'); %(k-1)*numel(vars(1).values)+j
        end
        str{numel(vars(2).values)} = strcat(str{k},' (m_{bat}=3.8kg)');

        xlabel('Day of year [-]','FontSize',FontSize);
        ylabel('Charge margin T_{cm} [h]','FontSize',FontSize);

        for i = 1:numel(ax)
            set(ax(i),'FontSize',FontSize);
            PlotEveryXthTick = 10;
            set(ax(i),'XTick',floor(vars(3).values(1:PlotEveryXthTick:end)+datevector(1:PlotEveryXthTick:end)));
            set(ax(i),'XTickLabel',datestr(floor(vars(3).values(1:PlotEveryXthTick:end)+datevector(1:PlotEveryXthTick:end)),'mmmdd'));
            xlim(ax(i),[vars(3).values(1)+datevector(1) vars(3).values(end)+datevector(end)]);
            ax(i).XTickLabelRotation=45;
        end
        legend(str,'FontSize',10);
    elseif(plotnr==4)
        %Configure figure
        figure('Name','Sensitivity plot','OuterPosition',[120 0 1200 680]);
        set(gcf, 'Color', 'w');
        FontSize = 14;
        hmargin = 0.115;

        %Convert data first
        k=1;
        for i = 1:numel(vars)
            for j = 1:numel(vars(1).values)
                %display([i,' ',k,' ',j]);
                temp_t_exc(i,k,j) = PerfResults(i,k,j).t_excess;
                temp_min_SoC(i,k,j) = PerfResults(i,k,j).min_SoC;
                temp_chargemargin(i,k,j) = PerfResults(i,k,j).t_chargemargin;
                temp_endurance(i,k,j) = PerfResults(i,k,j).t_endurance;
                if(temp_t_exc(i,k,j)<0.001)
                    temp_t_exc(i,k,j) = NaN;
                    temp_min_SoC(i,k,j) = NaN;
                    temp_chargemargin(i,k,j) = NaN;
                end
                temp_m_total(i,k,j) = 0.0;
                if(~isempty(DesignResults)) DesignResults(i,k,j).m_no_bat + DesignResults(i,k,j).m_bat; end;
                %if(isnan(temp_chargemargin(i,k,j))) temp_chargemargin(i,k,j) = 0.0; end %Remove NaNs to have smooth contour plotting below
            end
        end

        StyleArr = {{'b'};{'g'};{'c'};{'r'};};

        %PLOT1 : Excess time.
        ax(1) = subplot_tight(1,2,1,[vmargin hmargin]);
        ax(end).XGrid='on'; ax(end).GridLineStyle=':';
        hold on
        for k=1:numel(vars)
                plot(vars(1).values, squeeze(temp_t_exc(k,1,:)),StyleArr{k}{1,1},'LineWidth',2); %,CLRArr(i),'LineStyle',StyleArr(i)
                hold all
        end
        title('Excess time T_{exc} [h]','FontSize',FontSize);
        xlabel('Relative change of technological parameter','FontSize',FontSize);
        ylabel('Absolute performance metric','FontSize',FontSize);
        curticks_percent=[cellstr(num2str(get(ax(end),'xtick')'*100-100,'%+g'))]; 
        new_ticks = [char(curticks_percent),char(ones(size(curticks_percent,1),1)*'%')];
        set(ax(end),'xticklabel',new_ticks); 
        %legend('e_{bat}/e_{bat}^{nom}', '\eta_{sm}/\eta_{sm}^{nom}', '\eta_{prop}/\eta_{prop}^{nom}', 'm_{dry}/m_{dry}^{nom}','Location','Northwest');
        
        ax(end+1) = axes('Position',ax(end).Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
        nomval = temp_t_exc(1,1,find(vars(1).values == 1.0));
        ax(end).YLim(1)=ax(end-1).YLim(1) / nomval;
        ax(end).YLim(2)=ax(end-1).YLim(2) / nomval;
        ax(end).YGrid='on'; ax(end).GridLineStyle=':';
        ax(end).XTick=[];
        ylabel(ax(end),'Relative performance metric change','FontSize',FontSize);
        curticks_percent=[cellstr(num2str(get(ax(end),'ytick')'*100-100,'%+g'))]; 
        new_ticks = [char(curticks_percent),char(ones(size(curticks_percent,1),1)*'%')];
        set(ax(end),'yticklabel',new_ticks);

        %PLOT2 : Charge margin
        plothandles = gobjects(numel(vars),1);
        ax(end+1) = subplot_tight(1,2,2,[vmargin hmargin]);
        ax(end).XGrid='on'; ax(end).GridLineStyle=':';
        hold on
        for k=1:numel(vars)
                plothandles(k) = plot(vars(1).values, squeeze(temp_chargemargin(k,1,:)),StyleArr{k}{1,1},'LineWidth',2); %,CLRArr(i),'LineStyle',StyleArr(i)
                hold all
        end
        title('Charge margin T_{cm} [h]','FontSize',FontSize);
        xlabel('Relative change of technological parameter','FontSize',FontSize);
        ylabel('Absolute performance metric','FontSize',FontSize);
        curticks_percent=[cellstr(num2str(get(ax(end),'xtick')'*100-100,'%+g'))]; 
        new_ticks = [char(curticks_percent),char(ones(size(curticks_percent,1),1)*'%')];
        set(ax(end),'xticklabel',new_ticks); 
        
        ax(end+1) = axes('Position',ax(end).Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
        nomval = temp_chargemargin(1,1,find(vars(1).values == 1.0));
        ax(end).YLim(1)=ax(end-1).YLim(1) / nomval;
        ax(end).YLim(2)=ax(end-1).YLim(2) / nomval;
        ax(end).YGrid='on'; ax(end).GridLineStyle=':';
        ax(end).XTick=[];
        ylabel(ax(end),'Relative performance metric change','FontSize',FontSize);
        curticks_percent=[cellstr(num2str(get(ax(end),'ytick')'*100-100,'%+g'))]; 
        new_ticks = [char(curticks_percent),char(ones(size(curticks_percent,1),1)*'%')];
        set(ax(end),'yticklabel',new_ticks); 
        
        legend(plothandles,'e_{bat}/e_{bat}^{nom}', '\eta_{sm}/\eta_{sm}^{nom}', '\eta_{prop}/\eta_{prop}^{nom}', 'm_{dry}/m_{dry}^{nom}','Location','Southeast');
        
        for i = 1:numel(ax)
            set(ax(i),'XTickMode','manual','XLimMode','manual');
            set(ax(i),'FontSize',FontSize-2);
        end
    end
    %--------------------------------------------------------------------------
    %--- DETAILED TIME PLOTS
    %--------------------------------------------------------------------------
    %DetailedAnalysisForm({designVars.AR_array,designVars.b_array,designVars.m_bat_array,results.flightdata_array,results.t_excess,results.t_endurance,results.t_chargemargin});

    %--------------------------------------------------------------------------
    %--- SIMPLE TIME PLOTS
    %--------------------------------------------------------------------------
    % Allows to choose one configuration and plot detailed data.
%     for i = 1:numel(vars(3).values)
%         for k = 1:numel(vars(2).values)
%             for j = 1:numel(vars(1).values)
%                 if(flightdata(i,k,j).b == 4.5 && flightdata(i,k,j).m_bat == 7.0)
%                     Plot_BasicSimulationTimePlot(flightdata(i,k,j), environment, params, plane);
%                 end
%             end
%         end
%     end
    
end     
function Plot_AirplaneDesign_Standard(PerfResults, DesignResults, environment, plane, parameters, flightdata, vars)
    %--------------------------------------------------------------------------
    %--- OVERVIEW PLOTS
    %--------------------------------------------------------------------------
    for i = 1:numel(vars(3).values)
        if(numel(vars(2).values)<2 ||numel(vars(1).values)<2) 
            display('ERROR: Not enough design variables specified, cannot plot results! Please specific at least the first two design variables.');
            return;
        end

        %Set figure title
        str = strcat(vars(3).shortname, '=', num2str(vars(3).values(i)));
        figure('Name',str);
        
        %Convert data first
        for k = 1:numel(vars(2).values)
            for j = 1:numel(vars(1).values)
                display([i,' ',k,' ',j]);
                temp_t_exc(k,j) = PerfResults(i,k,j).t_excess;
                temp_min_SoC(k,j) = PerfResults(i,k,j).min_SoC;
                temp_chargemargin(k,j) = PerfResults(i,k,j).t_chargemargin;
                temp_endurance(k,j) = PerfResults(i,k,j).t_endurance;
                temp_m_total(k,j) = DesignResults(i,k,j).m_no_bat + DesignResults(i,k,j).m_bat;
            end
        end
        
        %NaN-check to avoid plotting errors below
        if(sum(~isnan(temp_endurance))==0) temp_endurance(:,:) = 0.0; end
        if(sum(~isnan(temp_chargemargin))==0) temp_chargemargin(:,:) = 0.0; end
        
        colormap(jet);
        
        %PLOT1 : Excess time. Draws contour surface, and contour lines
        subplot(2,2,1)
        [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_t_exc,[0.01:0.05:max(max(max(temp_t_exc)))],'Linestyle','none');
        hold on
        [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,'LevelStep',2.0,'LineColor',[0 0 0],'ShowText','on');
        title(['Excess Time [h] (' vars(3).shortname '=' num2str(vars(3).values(i)) ')']);
        ylabel(vars(2).name);
        xlabel(vars(1).name);
        colorbar

        % PLOT2 : Endurance
        subplot(2,2,2)
        [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_endurance, [0.1:0.1:max(max(max(temp_endurance)))],'Linestyle','none');
        hold on
        [c2,hc2]=contour(vars(1).values,vars(2).values,temp_endurance,'LevelStep',2.0,'LineColor',[0 0 0],'ShowText','on');
        title(['Endurance [h] (' vars(3).shortname '=' num2str(vars(3).values(i)) ')']);
        ylabel(vars(2).name);
        xlabel(vars(1).name);
        colorbar
        
        % PLOT3 : Charge margin
        subplot(2,2,3)
        [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_chargemargin, [0.01:0.01:max(max(max(temp_chargemargin)))],'Linestyle','none');
        hold on
        [c2,hc2]=contour(vars(1).values,vars(2).values,temp_chargemargin,'LevelStep',2.0,'LineColor',[0 0 0],'ShowText','on');
        title(['Charge margin [h] (' vars(3).shortname '=' num2str(vars(3).values(i)) ')']);
        ylabel(vars(2).name);
        xlabel(vars(1).name);
        colorbar
        
        %PLOT4: Optional plots can be chosen here
        
    %     PLOT4.1: Structural Mass
    %     subplot(2,2,4)
    %     contourf(designVars.b_array,...
    %         designVars.m_bat_array, m_struct(:,:,l)',200,...
    %         'Linestyle','none');
    %     xlabel('Wingspan (b) [m]')
    %     ylabel('Battery mass (m\_bat) [kg]');
    %     title(['Structural Mass [kg] (AR=' num2str(AR) ')']);
    %     caxis([min(min( m_struct(:,:,l))),max(max(m_struct(:,:,l)))])
    %     colorbar

%         % PLOT4.1: Combined Performace (T_exc+T_ChargeMargin)
%         subplot(2,2,4)
%         a=0.6;
%         b=1.0-a;
%         comb_perf=a*results.t_chargemargin+b*results.t_excess;
%         for i=1:size(results.flightdata_array,1)
%             for j=1:size(results.flightdata_array,2)
%                 for k=1:size(results.flightdata_array,3)
%                     if results.flightdata_array(i,j,k).minSoC<0.1
%                         comb_perf(i,j,k)=NaN;
%                     end
%                 end
%             end
%         end
% 
%         contourf(designVars.b_array,...
%             designVars.m_bat_array, comb_perf(:,:,l)',200,...
%             'Linestyle','none');
%         xlabel('Wingspan (b) [m]')
%         ylabel('Battery mass (m\_bat) [kg]');
%         title(['a*{t_cm}+b*t_{exc} [kg] (AR=' num2str(AR) ')']);
%         caxis([min(min( comb_perf(:,:,l))),max(max(comb_perf(:,:,l)))])
%         colorbar

        % PLOT4.3: Total airplane mass
        subplot(2,2,4)
        [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_m_total, [0.01:0.01:max(max(max(temp_m_total)))],'Linestyle','none');
        hold on
        [c2,hc2]=contour(vars(1).values,vars(2).values,temp_m_total,'LevelStep',1.0,'LineColor',[0 0 0],'ShowText','on');
        title(['Total mass [kg] (' vars(3).shortname '=' num2str(vars(3).values(i)) ')']);
        ylabel(vars(2).name);
        xlabel(vars(1).name);
        colorbar
    end

    %--------------------------------------------------------------------------
    %--- DETAILED TIME PLOTS
    %--------------------------------------------------------------------------
    %DetailedAnalysisForm({designVars.AR_array,designVars.b_array,designVars.m_bat_array,results.flightdata_array,results.t_excess,results.t_endurance,results.t_chargemargin});

    %--------------------------------------------------------------------------
    %--- SIMPLE TIME PLOTS
    %--------------------------------------------------------------------------
    Plot_BasicSimulationTimePlot(flightdata(end,end,end), environment, parameters, plane);
    
end     
function Plot_AirplaneAnalysis_Standard(results, vars)
    %**************************************************************************
    %*** PLOTTING
    %**************************************************************************
    % Note: This plot code will automatically plot a) multiple contour plots if
    %       3 variables are specified by the user, b) one contour plot for two
    %       variables, and c) simple 1D-plots for only one specified variable.
    today = now;
    doy = today - datenum(2015,1,1) + 1;

    for i = 1:numel(vars(3).values)

        %Convert first
        for k = 1:numel(vars(2).values)
            for j = 1:numel(vars(1).values)
                temp_t_exc(k,j) = results(i,k,j).t_excess;
                temp_min_SoC(k,j) = results(i,k,j).min_SoC;
                temp_chargemargin(k,j) = results(i,k,j).t_chargemargin;
                temp_endurance(k,j) = results(i,k,j).t_endurance;
            end
        end

        str = strcat(vars(3).shortname, '=', num2str(vars(3).values(i)));
        figure('Name',str);

        % PLOT1: T_EXCESS
        subplot(2,2,1);
        if(numel(vars(2).values)>1)
            contourf(vars(1).values,vars(2).values,temp_t_exc,[0.01:0.05:min(max(max(max(temp_t_exc))),24.0)],'Linestyle','none');
            hold on
            contour(vars(1).values,vars(2).values,temp_t_exc,'LevelStep',1.0,'LineColor',[0 0 0],'ShowText','on');
            ylabel(vars(2).name);
            colormap('jet');
            %colormap(flipud(colormap(hot))); %Set an inverted version of the "hot" colormap
            caxis([0.5,min(max(max(max(temp_t_exc))),24.0)])
            colorbar
        else
            plot(vars(1).values,temp_t_exc);
        end
        title('t_{excess}');
        xlabel(vars(1).name);

        % PLOT1: MIN_SOC
        subplot(2,2,2);
        if(numel(vars(2).values)>1)
            contourf(vars(1).values,vars(2).values,temp_min_SoC,[0.01:0.01:max(max(max(temp_min_SoC)))],'Linestyle','none');
            hold on
            contour(vars(1).values,vars(2).values,temp_min_SoC,'LevelStep',0.1,'LineColor',[0 0 0],'ShowText','on');
            ylabel(vars(2).name);
            colormap('jet');
            colorbar
        else
            plot(vars(1).values,temp_min_SoC);
        end
        title('SoC_{min}');
        xlabel(vars(1).name);

        % PLOT1: CHARGE MARGIN
        subplot(2,2,3);
        if(numel(vars(2).values)>1)
            contourf(vars(1).values,vars(2).values,temp_chargemargin,[0.01:0.05:min(max(max(max(temp_chargemargin))),24.0)],'Linestyle','none');
            hold on
            contour(vars(1).values,vars(2).values,temp_chargemargin,'LevelStep',1.0,'LineColor',[0 0 0],'ShowText','on');
            ylabel(vars(2).name);
            colormap('jet');
            colorbar
        else
            plot(vars(1).values,temp_chargemargin);
        end
        title('t_{chargemargin}');
        xlabel(vars(1).name);

        % PLOT1: ENDURANCE
        subplot(2,2,4);
        if(numel(vars(2).values)>1)
            contourf(vars(1).values,vars(2).values,temp_endurance,[0.01:0.5:max(0.51,max(max(max(temp_endurance))))],'Linestyle','none');
            hold on
            contour(vars(1).values,vars(2).values,temp_endurance,'LevelStep',1.0,'LineColor',[0 0 0],'ShowText','on');
            ylabel(vars(2).name);
            colormap('jet');
            colorbar
        else
            plot(vars(1).values,temp_endurance);
        end
        title('t_{endurance}');
        xlabel(vars(1).name);
    end
end


function quality=calculate_polar(airfoil_name,alpha_start,alpha_end,step,re)

%===============================================================
%=== DESCRIPTION Version 1.00
%===
%=== This function uses the calculate_polar_angle function to 
%=== compute the Lift and Drag coefficient of a profile at a
%=== Reynolds number re and for an angle interval. For all angles,
%=== the computation is tried first with few iterations, in order
%=== to save time. If no solution was found, the number of
%=== iterations is increased. A *.dat file containing the values
%=== in 3 columns (angle CL CD) is created in the same folder.
%=== Moreover, a *.png picture file that includes plots of polars
%=== CL=f(CD) and CL=f(angle) is also created. All those files
%=== takes the name of the airfoil concatenated with the reynolds
%=== number. Example: polar of airfoil "goe795.dat" calculated at
%=== re = 60000 -> goe795-60000.dat goe795-60000.png
%=== A quality factor is defined as the ratio of successful angles
%=== and the tried angles.
%===
%=== Note: in some cases, with extreme airfoils, crashes of Xfoil
%===       may happen. A timer was added in order to kill the
%===       Xfoil task (here XfoilP4.exe) with force (taskkill /f)
%===
%=== Input:  airfoil_name  file_name, without *.dat extension
%===         alpha_start   start angle
%===         alpha_end     end angle (note: start angle < end angle
%===         step          angle step (In general, 0.5 sufficient)
%===         re            reynolds number
%===
%=== Output: quality       Float between 0 and 1
%===
%=== Files created: *.dat  containing (angle CL CD) values
%===                *.png  polar pictures CL=f(CD) & CL=f(angle)
%===
%=== Dependencies: needs "calculate_polar_angle" function
%===
%=== Author: André Noth
%===         andre.noth@epfl.ch
%===         http://aero.epfl.ch
%===
%=== More info:  Xfoil program
%===             http://raphael.mit.edu/xfoil/
%===============================================================

%===============================================================
%=== Initialization
%===============================================================

polar_name = airfoil_name;                                                  % Polar will have the same name than the airfoil
polar_filename=strcat(strcat(polar_name,strcat('-',num2str(re))),'.dat');   % String defining the filename of the polar (adding Reynolds number and .dat)
fid2=fopen(polar_filename,'wt');                                            % Opens the new file that will be filled in -> File descriptor
array_index     = 1;                                                        % 
point_number    = 0;                                                        %
CL_array        = [];                                                       % Array that will contain all the values of Lift Coefficient
CD_array        = [];                                                       % Array that will contain all the values of Drag Coefficient
CM_array        = [];
alpha_array     = [];                                                       % Array that will contain all the values of angle


%===============================================================
%=== Calculation and storage of the values in *.dat
%===============================================================

for i=alpha_start:step:alpha_end
    disp(strcat('Airfoil: ',strcat(polar_name,strcat(' Angle :',num2str(i))))); % Displays the name of the airfoil actually computed
    crash_timer = timer('TimerFcn','!taskkill /im xfoilP4.exe /f', 'StartDelay', 5.0); % Prepare a timer in case of crash...
    start(crash_timer);                                                         % and start it
    [success,CL,CD,CM]=calculate_polar_angle(polar_name,i,re,2);                   % Try to calculate with 2 iterations...
    stop(crash_timer);                                                          % Stop the timer
    if success == 0                                                             % If computation with 2 iterations not succesful...
        crash_timer = timer('TimerFcn','!taskkill /im xfoilP4.exe /f', 'StartDelay', 5.0);
        start(crash_timer);
        [success,CL,CD,CM]=calculate_polar_angle(polar_name,i,re,10);              % Try to calculate with 10 iterations...
        stop(crash_timer);
    elseif success == 0                                                         % If computation with 10 iterations not succesful...
        crash_timer = timer('TimerFcn','!taskkill /im xfoilP4.exe /f', 'StartDelay', 10.0);
        start(crash_timer);
        [success,CL,CD,CM]=calculate_polar_angle(polar_name,i,re,50);              % Try to calculate with 50 iterations...
        stop(crash_timer);
    elseif success == 0                                                         % If computation with 50 iterations not succesful...
        crash_timer = timer('TimerFcn','!taskkill /im xfoilP4.exe /f', 'StartDelay', 15.0);
        start(crash_timer);
        [success,CL,CD,CM]=calculate_polar_angle(polar_name,i,re,200);             % Try to calculate with 200 iterations ;-)
        stop(crash_timer);
    end
    
    if (success == 1)                                                           % If the computation gave a solution for this angle...
        fprintf(fid2,strcat(num2str(i,'%.3f'),'\t'));                           % we write the angle...
        fprintf(fid2,strcat(num2str(CL,'%.3f'),'\t'));                          % the lift coefficient and...
        fprintf(fid2,strcat(num2str(CD,'%.3f'),'\t'));                          % the drag coefficient
        fprintf(fid2,strcat(num2str(CM,'%.3f'),'\n'));                          % the moment coefficient
        
        alpha_array(array_index)=i;                                             % For the plot later, we store the angle...
        CL_array(array_index)=CL;                                               % the lift coefficient and...
        CD_array(array_index)=CD;                                               % the drag coefficient
        CM_array(array_index)=CM;
        array_index=array_index+1;
        point_number = point_number + 1;                                        % We count the number of successful angles
        success = 0;                                                            % Success is re-initialized
    end
end

fclose(fid2);
quality = point_number / ((alpha_end-alpha_start)/step +1);                     % A quality is computed based on the ideal~/obtained number of points

%===============================================================
%=== Plot of polars and creation of *.png image files
%===============================================================

if (quality > 0)                                                                % If there is at least one point (quality>0)
    
    figure(1);
    
    subplot(1,2,1);
    plot([-10,20],[0,0],'-k','LineWidth',2);
    hold on;
    plot([0,0],[-0.5,2],'-k','LineWidth',3);
    plot(CD_array,CL_array,'-bs','LineWidth',1);
    axis([0 0.04 -0.5 2]);
    Title(strcat('Polar curves of airfoil ',strcat(':',polar_name)));
    set(gca,'XGrid','on')
    set(gca,'XTick',0:0.002:0.04)
    set(gca,'XTickLabel',{'0';'';'';'';'';'0.01';'';'';'';'';'0.02';'';'';'';'';'0.03';'';'';'';'';'0.04';'';'';'';'';})
    xlabel('CD');
    set(gca,'YGrid','on')
    set(gca,'YTick',-0.5:0.1:2)
    set(gca,'YTickLabel',{'-0.5';'';'';'';'';'0.0';'';'';'';'';'0.5';'';'';'';'';'1';'';'';'';'';'1.5';'';'';'';'';'2'})
    ylabel('CL');
    
    subplot(1,2,2);
    plot([-10,20],[0,0],'-k','LineWidth',2);
    hold on;
    plot([0,0],[-0.5,2],'-k','LineWidth',2);
    plot(alpha_array,CL_array,'-bs','LineWidth',1);
    Title(strcat('Re number ',strcat(':',num2str(re))));
    axis([-10 20 -0.5 2]);
    set(gca,'XGrid','on')
    set(gca,'XTick',-10:2:20)
    set(gca,'XTickLabel',{'-10';'';'';'';'';'0';'';'';'';'';'10';'';'';'';'';'20'})
    xlabel('Alpha');
    set(gca,'YGrid','on')
    set(gca,'YTick',-0.5:0.1:2)
    set(gca,'YTickLabel',{'-0.5';'';'';'';'';'0.0';'';'';'';'';'0.5';'';'';'';'';'1';'';'';'';'';'1.5';'';'';'';'';'2'})
    ylabel('CL');
    
    saveas(gcf,strcat(strcat(polar_name,strcat('-',num2str(re))),'.png'));
    close(gcf);
    
end
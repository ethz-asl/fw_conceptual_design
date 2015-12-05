function [success,CL,CD,CM] = calculate_polar_angle(airfoil_name,alpha,re,iteration)

%===============================================================
%=== DESCRIPTION Version 1.00
%===
%=== This function uses Xfoil program to calculate the Lift and
%=== Drag coefficient for an airfoil at an incidence angle and
%=== with a reynold number re. The various command that will be
%=== used to interact with Xfoil are stored in a file. The airfoil
%=== used for the calculation should be in folder /Airfoils/
%=== The Xfoil program should be in folder /Xfoil
%===
%=== Input:  airfoil_name  file_name, without *.dat extension
%===         alpha         incidence angle in degree
%===         re            reynolds number
%===
%=== Output: success       1 or 0 if calculation was possible or not
%===         CL            Lift coefficient
%===         CD            Drag coefficient
%===
%=== Files created: none
%===
%=== Dependencies: none (except Xfoil)
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

CL = 0;                                                                     % Initializes CL to 0
CD = 0;                                                                     % Initializes CD to 0
success = 0;
cd Xfoil;                                                                   % Goes to the Xfoil directory
if exist('polars.key')                                                      % If a command file already exists, we delete it
    !del polars.key
end

%===============================================================
%=== Interface to the Xfoil Program
%===============================================================

file=strcat('polars','.key');                                               % The consecutive keys for xfoil are stored in a .key file
fid=fopen(file,'wt');                                                       % Creates the file
fprintf(fid,'plop \n');                                                     %
fprintf(fid,'w \n');                                                        %
fprintf(fid,'0.2 \n');                                                      % Reduces the size of the window where plots appear
fprintf(fid,'\n');                                                          %
fprintf(fid,strcat('load ../Airfoils/',strcat(airfoil_name,'.dat \n')));    % Loads the airfoil
fprintf(fid,'oper \n');                                                     % Enters the operation mode
fprintf(fid,'visc \n');                                                     % Takes care of the viscosity of the fluid
fprintf(fid,strcat(num2str(re),'\n'));                                      % Sets the Reynolds number
fprintf(fid,'pacc \n');                                                     % Start recording the values...
fprintf(fid,'polar_1\n');                                                   % ... in a polar list named polar_1
fprintf(fid,'\n');                                                          %

for i=0:sign(alpha):alpha-sign(alpha)*mod(alpha,1)-sign(alpha)              % For high angles (>10°), crash of Xfoil can happen. It...
    fprintf(fid,strcat(strcat('a ',num2str(i)),'\n'));                      % can be avoided by going step by step to the value...
end                                                                         % calculating intermediate values

for j=0:1:iteration                                                         % The calculation of CL and CD for one angle needs many...
    fprintf(fid,strcat(strcat('a ',num2str(alpha)),'\n'));                  % iterations until the solution converges
end

fprintf(fid,'pwrt\n');                                                      % Save the polar point in a file...
fprintf(fid,strcat('polar_transition_file.dat','\n'));                      % named 'polar_transistion_file.dat'
fprintf(fid,' \n');                                                         %
fprintf(fid,'quit\n');                                                      % Xfoil is quit
fclose(fid);

!xfoilP4.exe<polars.key %>a

%===============================================================
%=== Treatment of the polar file and extraction of CD and CL
%===============================================================

if exist('polar_transition_file.dat')
    fid=fopen('polar_transition_file.dat','a+');                            % Open a transition file
    fseek(fid, 1, 'bof');                                                   % Go the beginning
    tline       = fgetl(fid);                                               % Get the line
    tline_last	= fgetl(fid);                                               % Get the last line

    while tline ~= -1                                                       % Go through all the file until the last line where the...
        tline_last = tline;                                                 % values corresponding to our angle is stored.
        tline = fgetl(fid);
    end

    if findstr(tline_last,'------')                                         % It means that no data are stored in the file...
        success = 0;                                                        % which means that Xfoil had didn't succeed in the calculation
    else
        [alpha_read, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr]=strread(tline_last,'%f %f %f %f %f %f %f');
        success = (alpha_read == alpha);                                    % Success is 1 if Xfoil calculated CD CL for the desired angle
    end
    fclose(fid);                                                            % We close the transition file
else                                                                        % If the transition file doesn't exist...
    success = 0;                                                            % it means that Xfoil had problems and stopped
    CL = 0;
    CD = 0;
    CM = 0;
end

%===============================================================
%=== Delete the files no more useful
%===============================================================

if exist('polar_transition_file.dat')
    !del polar_transition_file.dat
end
if exist('polar_1')
    !del polar_1
end
if exist('polars.key')
    !del polars.key
end
if exist('a')
    !del a
end
cd('..')

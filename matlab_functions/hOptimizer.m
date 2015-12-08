% find h-optimizing solution

options = optimset('Display','iter','TolFun',2,'TolX',1e-2);
initParameters;

global environment;

% Set environment, payload and airfoil
environment.month = 6;
environment.day = 21;
environment.h = 700;
environment.T_ground = 300;
%%environment.lat = 47;
environment.albedo = 0.2;
environment.clearness = 1;
payload.mass = 0.6;%0.2;%0.722;
payload.power = 4;%1.7;%18.4;

% Change some parameters
params.propulsion.number  =  2;        % Number of propulsion units [-]
params.structure.shell    =  0;        % 1 for shell wing, 0 for rib wing
params.evaluation.clmb    =  0;        % 1 to allow altitude changes
params.evaluation.findalt =  1;        % if 1, it finds the maximum
                                           % altitude for eternal flight
                                           
x0=[13.964754649863533  21.309585271646178  21.238851828246155];



lat_array=0:1:90;
b=zeros(length(lat_array),1);
AR=zeros(length(lat_array),1);
m_bat=zeros(length(lat_array),1);
v_tmax=zeros(length(lat_array),1);
h_max=zeros(length(lat_array),1);
m_struct=zeros(length(lat_array),1);

for (i=1:length(lat_array))
    disp(lat_array(i));
    environment.lat=lat_array(i);
    [x,fval,exitflag,output] = ...
        fminsearch(@evaluateSolutionSearch,x0,options);
    [performance,polar,masses] ...
        = evaluateSolution(x(1),x(2),x(3),payload,environment,params);
    b(i)=x(1);
    AR(i)=x(2);
    m_bat(i)=x(3);
    v_tmax(i)=performance.v_tmax;
    h_max(i)=performance.h;
    m_struct(i)=masses.m_struct;
    x0=x;
    store = ([lat_array;b;AR;m_bat;v_tmax;h_max;m_struct]);
    save store;
end

load store3
subplot(5,1,1)
plot(store3(1:68,1),store3(1:68,6)/1000)
ylabel('Maximum altitude [km]')
subplot(5,1,2)
plot(store3(1:68,1),store3(1:68,2))
ylabel('Wingspan b [m]')
subplot(5,1,3)
plot(store3(1:68,1),store3(1:68,3))
ylabel('Aspect ratio \Lambda [-]')
subplot(5,1,4)
plot(store3(1:68,1),store3(1:68,5))
ylabel('True airspeed [m/s]')
subplot(5,1,5)
plot(store3(1:68,1),store3(1:68,4))
ylabel('Battery mass [kg]')
xlabel('Latitude [°N]')
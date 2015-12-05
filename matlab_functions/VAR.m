% Definition of an enumeration of possible design/analysis variables
classdef VAR < handle
    properties
        %type;
        values;
    end
    enumeration
        UNDEFINED
        WING_SPAN
        BATTERY_MASS
        ASPECT_RATIO
        DAY_OF_YEAR
        LATITUDE
        CLEARNESS
        TURBULENCE
        POWER
    end
    methods
       function r = name(obj)
         if(strcmp(char(obj),'UNDEFINED'))
             r = 'UNDEFINED';
         elseif(strcmp(char(obj),'WING_SPAN'))
             r = 'Wing span b [m]';
         elseif(strcmp(char(obj),'BATTERY_MASS'))
             r = 'Battery mass m_{bat} [kg]';
         elseif(strcmp(char(obj),'ASPECT_RATIO'))
             r = 'Aspect ratio \lambda [-]';
         elseif(strcmp(char(obj),'DAY_OF_YEAR'))
             r = 'Day of year [-]';
         elseif(strcmp(char(obj),'LATITUDE'))
             r = 'Latitude \phi [°]';
         elseif(strcmp(char(obj),'CLEARNESS'))
             r = 'Atmospheric clearness [-]';
         elseif(strcmp(char(obj),'TURBULENCE'))
             r = 'Atmospheric turbulence [-]';
         elseif(strcmp(char(obj),'POWER'))
             r = 'Level flight el. prop. power [W]';
         end
       end
       
       function r = shortname(obj)
         if(strcmp(char(obj),'UNDEFINED'))
             r = 'UNDEFINED';
         elseif(strcmp(char(obj),'WING_SPAN'))
             r = 'Wing span';
         elseif(strcmp(char(obj),'BATTERY_MASS'))
             r = 'Battery mass';
         elseif(strcmp(char(obj),'ASPECT_RATIO'))
             r = 'Aspect ratio';
         elseif(strcmp(char(obj),'DAY_OF_YEAR'))
             r = 'Day of year';
         elseif(strcmp(char(obj),'LATITUDE'))
             r = 'Latitude';
         elseif(strcmp(char(obj),'CLEARNESS'))
             r = 'Clearness';
         elseif(strcmp(char(obj),'TURBULENCE'))
             r = 'Turbulence';
         elseif(strcmp(char(obj),'POWER'))
             r = 'Level flight el. prop. power';
         end
       end
    end
end
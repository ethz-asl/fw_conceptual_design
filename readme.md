# Conceptual Design and Analysis Tool for solar-powered UAVs #

## Based on ##
 - Philipp Oettershagen, Amir Melzer, Thomas Mantel, Konrad Rudin, Thomas Stastny, Bartosz Wawrzacz, Timo Hinzmann, Stefan Leutenegger, Kostas Alexis, Roland Siegwart, _Design of small hand‐launched solar‐powered UAVs: From concept study to a multi‐day world endurance record flight_, Journal of Field Robotics 34 (7), 1352-1377 . http://www.atlantiksolar.ethz.ch/wp-content/downloads/publications/JFR_81hFlight_paper_final.pdf
 - Philipp Oettershagen, Amir Melzer, Thomas Mantel, Konrad Rudin, Rainer Lotz, Dieter Siebenmann, Stefan Leutenegger, Kostas Alexis and Roland Siegwart, _A Solar-Powered Hand-Launchable UAV for Low-Altitude Multi-Day Continuous Flight_, International Conference on Robotics and Automation (ICRA) 2015. http://www.atlantiksolar.ethz.ch/wp-content/downloads/publications/AtlantikSolar_ICRA_2015_vFinal.pdf
 - Stefan Leutenegger, Mathieu Jabas, Roland Y. Siegwart, _Solar Airplane Conceptual Design and Performance Estimation_, Journal of Intelligent & Robotic Systems (2011), Volume 61, Issue 1, pp 545-561.
 - A. Noth, _Design of solar powered airplanes for continuous flight_, PhD thesis, ETH Zurich, 2008. http://www.asl.ethz.ch/research/asl/skysailor/Design_of_Solar_Powered_Airplanes_for_Continuous_Flight

## Usage ##
Run
 - `AirplaneDesign.m` to design your airplane (design variables are wing span b, aspect ratio AR, and battery mass m_bat) and analyse its performance (in the form of excess time t_exc, charge margin t_cm, endurance t_endurance, and minimum state-of-charge SoC).
 - `AirplaneAnalysis.m` to analyse your designed airplane with respect to a) other days of the year or other latitudes or b) meteorological disturbances (clouds or winds) in the form of the clearness and turbulence values.
    

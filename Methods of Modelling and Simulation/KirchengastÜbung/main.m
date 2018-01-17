% Uebung Kirchengast
% Felix Sonnleitner, 01430166
clc;
clear all;
close all;

Re     = 6371;           %km
z      = 650;            %km
theta0 = 0;              %rad
phi0   = [0 30 60 89.9] .* pi/180; %rad
N0     = 315;            %N-units
Hn     = 7;              %km
r0     = Re;             %km

begin = r0;
stop  = Re + z;
h = [0.1 1];

polarPlot            = true;
cartesianPlotFull    = true;
cartesianPlotReduced = true;
errorPlot            = true;

diary('Output.txt')
diary on
disp('-----------------------------START--------------------------')
for n = 1:numel(h)
    step  = h(n);
    
    % Calculation Euler and RK
    for i = 1:numel(phi0)
        [theta_euler(i,:), r_euler(i,:)] = f_euler(phi0(i), Hn, N0, theta0, begin, step, stop);
        [theta_rk(i,:), r_rk(i,:)]       = f_rungeKutta(phi0(i), Hn, N0, theta0, begin, step, stop);
        
        euler_curve(i,n) = theta_euler(i,end) * Re;
        rk_curve(i,n)    = theta_rk(i,end) * Re;
        disp(['Angle phi = ', num2str(phi0(i) * (180/pi)])
        disp(['horizontal distance d for euler beams = ', num2str(euler_curve(i,n)),' km for h = ', num2str(step)])
        disp(['horizontal distance d for runge-kutta beams = ', num2str(rk_curve(i,n)),' km for h = ', num2str(step)])
    end

    diff = abs((theta_rk - theta_euler)./theta_rk);
    
    % polar plot
    if polarPlot
        figure
        circle = linspace(0,2*pi,500);
        % Earthsurface plot
        polarplot(circle,Re.*ones(size(circle))-Re, 'o')
        title(['Polarplot: h = ',num2str(step)])
        ax = gca;
        ax.ThetaAxisUnits = 'degrees';
        ax.ThetaLim = [-30 30];
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'top';
        hold on
        % polarplot for each phi
        for line = 1:numel(phi0)
            polarplot(theta_euler(line,:), r_euler(line,:)-Re, '--')
            polarplot(theta_rk(line,:), r_rk(line,:)-Re)
        end
        % Satelliteorbit plot
        polarplot(circle,z.*ones(size(circle)))
        rlim([r_euler(1)-300-Re  r_euler(end)-Re+50]);
        rticks([r_euler(1)-300-Re 0 r_euler(end)-Re])
        thetaticks(-30:10:30)
        rticklabels({['r = ',num2str(r_euler(1)-300), ' km'], ['r = ', num2str(Re), ' km'],['r = ',num2str(r_euler(end)), ' km']})
        hold off
        legendstring = {'R_E = r_0', ...
                        '\phi_{0,E} = 0','\phi_{0,RK} = 0',...
                        '\phi_{0,E} = 30','\phi_{0,RK} = 30',...
                        '\phi_{0,E} = 60','\phi_{0,RK} = 60',...
                        '\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9',...
                        'Satellitenorbit'};
        legend(legendstring)
        print(['Polarplot_', num2str(n)],'-dpdf')
    end
    
   
   % cartesian plot
   if cartesianPlotFull
        figure
        hold on
        title(['Kartesischer Plot: h = ',num2str(step)])
        for line = 1:numel(phi0)
            d_euler = theta_euler(line,:) .* Re;
            d_rk    = theta_rk(line,:) .* Re;
            plot(r_euler(line,:), d_euler, '--')
            plot(r_rk(line,:), d_rk)
        end
        hold off
        xlabel('Radius [km]')
        ylabel('Abweichung [km]')
        legendstring2 = {
                        '\phi_{0,E} = 0','\phi_{0,RK} = 0',...
                        '\phi_{0,E} = 30','\phi_{0,RK} = 30',...
                        '\phi_{0,E} = 60','\phi_{0,RK} = 60',...
                        '\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9'
                        };
        legend(legendstring2, 'location', 'northwest')
        print(['Kartesisch_', num2str(n)],'-dpdf')
   end
   
   % cartesian plot limited
   if cartesianPlotReduced
        figure
        hold on
        title(['Kartesischer Plot fuer die ersten 30 km: h = ',num2str(step)])
        for line = 1:numel(phi0)
            d_euler = theta_euler(line,1:30/step) .* Re;
            d_rk    = theta_rk(line,1:30/step) .* Re;
            plot(r_euler(line,1:30/step), d_euler, '--')
            plot(r_rk(line,1:30/step), d_rk)
        end
        xlabel('Radius [km]')
        ylabel('Abweichung [km]')
        hold off
        legendstring2 = {
                        '\phi_{0,E} = 0','\phi_{0,RK} = 0',...
                        '\phi_{0,E} = 30','\phi_{0,RK} = 30',...
                        '\phi_{0,E} = 60','\phi_{0,RK} = 60',...
                        '\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9'
                        };
        legend(legendstring2)
        print(['Kartesisch_reduziert_', num2str(n)],'-dpdf')
   end
    
    % Error
    if errorPlot
        figure
        hold on
        title(['Errorplot (RK-Euler)/(RK) fuer die ersten 30 km: h = ',num2str(step)])
        for line = 1:numel(phi0)
            plot(r_euler(line,1:30/step), diff(line,1:30/step))
        end
        xlabel('Radius [km]')
        ylabel('Relativer Fehler [km/km]')
        hold off
        legendstring3 = {
                        '\phi_{0} = 0',...
                        '\phi_{0} = 30',...
                        '\phi_{0} = 60',...
                        '\phi_{0} = 89.9'
                        };
        legend(legendstring3)
        print(['Error_', num2str(n)],'-dpdf')
    end
    
    if numel(h) > 1
        % empty matrices
        theta_euler = [];
        r_euler     = [];
        theta_rk    = [];
        r_rk        = [];
    end
end
disp('----------------------------STOP--------------------------')
disp('')
diary off

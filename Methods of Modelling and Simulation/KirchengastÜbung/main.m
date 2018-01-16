% Ãœbung Kirchengast
% Felix Sonnleitner, 01430166
clc;
clear all;
close all;

Re     = 6371;           %km
z      = 650;            %km
theta0 = 0;              %Grad
phi0   = [0 30 60 89.9]; %Grad
N0     = 315;            %N-units
Hn     = 7;              %km
r0     = Re;             %km

begin = r0;
stop  = Re + z;
h = [0.01 0.1 1];

polarPlot            = false;
cartesianPlotFull    = false;
cartesianPlotReduced = false;
errorPlot            = true;

for n = 1:numel(h)
    step  = h(n);
    
    % Calculation Euler and RK
    for i = 1:numel(phi0)
        [theta_euler(i,:), r_euler(i,:)] = f_euler(phi0(i), Hn, N0, theta0, begin, step, stop);
        [theta_rk(i,:), r_rk(i,:)]       = f_rungeKutta(phi0(i), Hn, N0, theta0, begin, step, stop);
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
        ax.ThetaLim = [-30 30];
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'top';
        hold on
        % polarplot for each phi
        for line = 1:numel(phi0)
            polarplot(theta_euler(line,:), r_euler(line,:).'-Re, '--')
            polarplot(theta_rk(line,:), r_rk(line,:).'-Re)
        end
        % Satelliteorbit plot
        polarplot(circle,z.*ones(size(circle)))
        rlim([r_euler(1)-300-Re  r_euler(end)-Re]);
        hold off
        legendstring = {'R_E = r_0', ...
                        '\phi_{0,E} = 0','\phi_{0,RK} = 0',...
                        '\phi_{0,E} = 30','\phi_{0,RK} = 30',...
                        '\phi_{0,E} = 60','\phi_{0,RK} = 60',...
                        '\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9',...
                        'Satellitenorbit'};
        legend(legendstring)
    end
    
   
   % cartesian plot
   if cartesianPlotFull
        figure
        hold on
        title(['Kartesischer Plot: h = ',num2str(step)])
        for line = 1:numel(phi0)
            plot(r_euler(line,:), theta_euler(line,:), '--')
            plot( r_rk(line,:), theta_rk(line,:))
        end
        hold off
        legendstring2 = {
                        '\phi_{0,E} = 0','\phi_{0,RK} = 0',...
                        '\phi_{0,E} = 30','\phi_{0,RK} = 30',...
                        '\phi_{0,E} = 60','\phi_{0,RK} = 60',...
                        '\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9'
                        };
        legend(legendstring2)
   end
   
   % cartesian plot limited
   if cartesianPlotReduced
        figure
        hold on
        title(['Kartesischer Plot für die ersten 30 km: h = ',num2str(step)])
        for line = 1:numel(phi0)
            plot(r_euler(line,1:30/step),theta_euler(line,1:30/step), '--')
            plot(r_rk(line,1:30/step), theta_rk(line,1:30/step))
        end
        hold off
        legendstring2 = {
                        '\phi_{0,E} = 0','\phi_{0,RK} = 0',...
                        '\phi_{0,E} = 30','\phi_{0,RK} = 30',...
                        '\phi_{0,E} = 60','\phi_{0,RK} = 60',...
                        '\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9'
                        };
        legend(legendstring2)
   end
    
    % Error
    if errorPlot
        figure
        hold on
        title(['Errorplot (RK-Euler)/(RK) für die ersten 30 km: h = ',num2str(step)])
        for line = 1:numel(phi0)
            plot(r_euler(line,1:30/step), diff(line,1:30/step))
        end
        hold off
        legendstring3 = {
                        '\phi_{0} = 0',...
                        '\phi_{0} = 30',...
                        '\phi_{0} = 60',...
                        '\phi_{0} = 89.9'
                        };
        legend(legendstring3)
    end
    
    if numel(h) > 1
        % empty matrices
        theta_euler = [];
        r_euler     = [];
        theta_rk    = [];
        r_rk        = [];
    end
end
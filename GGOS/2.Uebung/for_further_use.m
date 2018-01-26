% FOR FURTHER USE: PLOT
%                   
% xp           = (R/omega_N) .* omega_calibrated(1,:);
% yp           = (R/omega_N) .* omega_calibrated(2,:);
% 
% xp_reference = (R/omega_N) .* reference(1,1:length(omega_calibrated)).*3600;
% yp_reference = (R/omega_N) .* reference(2,1:length(omega_calibrated)).*3600;
% 
% figure(2)
% hold on
% plot(xp,yp)
% plot(xp_reference,yp_reference)
% hold off
% title('Polar Motion at Earth Surface')
% legend('omega corrected', 'reference data')
% ylabel('y[m]')
% xlabel('x[m]')
% axis equal
% savefig('CalibratedPolarPlot.fig')
              
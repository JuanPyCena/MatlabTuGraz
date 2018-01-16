%ray tracing
clear all
close all
%format long
RE = 6371; %Erdradius [km]
hi = 1; % Schrittweite [km]
z = 650; %[km]
colourstring = 'rbmk';
phi0_vec = [0,30,60,89.9]./180.*pi;

legendstring = {'\phi_{0} = 0','\phi_{0} = 30','\phi_{0} = 60','\phi_{0} = 89.9'};

r_span = RE:hi:(RE+z);
theta_span_euler = nan(numel(r_span),numel(phi0_vec));
theta_span_euler(1,:) = 0; %Anfangswert bei ri = RE
theta_span_runge = theta_span_euler;

n_r = indexOfRefrection(r_span);

figure
hold on
title('n(r)')
plot(r_span,n_r)
xlabel('r0 + r [km]')
ylabel('n(r)')
hold off


for jj = 1:numel(phi0_vec)
    phi0 = phi0_vec(jj);
    disp(['Angle phi0 = ',num2str(phi0/pi*180)])
    psi = asin(RE*sin(pi-phi0)/(RE+z));
    theta_straight_beam = pi-(pi-phi0)-psi;
    d_straight = theta_straight_beam*RE;
    disp(['horizontal distance d for straight beams = ',num2str(d_straight),' km'])
    
    
    
    
    
    C = indexOfRefrection(RE).*RE.*sin(phi0);
    phi_r = @(n_r,r) asin(C./(n_r.*r));
    dtheta_dr = @(n_r,r) C./(n_r.*(r.^2).*cos(phi_r(n_r,r)));



    for ii = 2:numel(r_span)
        ri = r_span(ii-1);
        %Integration: Euler methode
        %y(x + h) = y(x) + h*f(x,y(x))
        % x = ri
        % y = theta_span
        n_ri = indexOfRefrection(ri);
        k1 = hi.*dtheta_dr(n_ri,ri);
        y_x= theta_span_euler(ii-1,jj);
        theta_span_euler(ii,jj) = y_x + k1;    
        %Integration: Runge Kutta 4th order
        % k1 = h*f(x,y(x)) Gleich wie bei Euler
        % k2 = h*f(x+h/2,y(x)+k1/2)
        % k3 = h*f(x+h/2,y(x)+k2/2)
        % k4 = h*f(x+h,y(x)+k3)
        % y(x+h) = y(x) + k1/6 + k2/3 + k3/3 + k4/6
        k2 = hi.*dtheta_dr(indexOfRefrection(ri+hi/2),ri+hi/2);
        k3 = k2;
        k4 = hi.*dtheta_dr(indexOfRefrection(ri+hi),ri+hi);
        theta_span_runge(ii,jj) = y_x + k1./6 + k2./3 + k3./3 + k4./6;

    end  
    d_euler_curve = (theta_span_euler(end,jj)).*RE;
    disp(['horizontal distance d for euler beams = ',num2str(d_euler_curve),' km'])
    d_runge_curve = (theta_span_runge(end,jj)).*RE;
    disp(['horizontal distance d for runge-kutta beams = ',num2str(d_runge_curve),' km'])
end




figure
circle = linspace(0,2*pi,100);
polarplot(circle,RE.*ones(size(circle))-RE)
title(['h = ',num2str(hi)])
ax = gca;
ax.ThetaLim = [-30 30];
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
hold on
for jj = 1:numel(phi0_vec)
polarplot(theta_span_euler(:,jj),r_span.'-RE)
polarplot(theta_span_runge(:,jj),r_span.'-RE)
%ax = gca; % current axes
%ax.RLim = [r_span(1)-10 r_span(end)+10];

end
polarplot(circle,z.*ones(size(circle)))
rlim([r_span(1)-500-RE  r_span(end)-RE]);
legendstring2 = {'r_{0}','\phi_{0,E} = 0','\phi_{0,RK} = 0','\phi_{0,E} = 30','\phi_{0,RK} = 30','\phi_{0,E} = 60','\phi_{0,RK} = 60','\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9','Satellitenorbit'};
legend(legendstring2)
%legend({'r_{0}','\phi_{0} = 0','\phi_{0} = 30','\phi_{0} = 60','\phi_{0} = 89.9','Satellitenorbit'})
hold off
legendstring2 = {'\phi_{0,E} = 0','\phi_{0,RK} = 0','\phi_{0,E} = 30','\phi_{0,RK} = 30','\phi_{0,E} = 60','\phi_{0,RK} = 60','\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9'};
figure
hold on
title(['h = ',num2str(hi)])
for jj = 1:numel(phi0_vec)
    d_euler = (theta_span_euler(:,jj)).*RE;
    d_runge = (theta_span_runge(:,jj)).*RE;
    plot(d_euler,r_span,'-')
    plot(d_runge,r_span,'--')
end
xlabel('horizontale Distanz d [km]')
ylabel('r(theta)-r0 [km]')
legend(legendstring2)
hold off

figure
hold on
title(['nur r0+30km mit h = ',num2str(hi)])
for jj = 1:numel(phi0_vec)
d_euler = (theta_span_euler(:,jj)).*RE;
d_runge = (theta_span_runge(:,jj)).*RE;

plot(d_euler((r_span-RE)<=30),r_span((r_span-RE)<=30),'-')
plot(d_runge((r_span-RE)<=30),r_span((r_span-RE)<=30),'--')
end
legend(legendstring2)
xlabel('horizontale Distanz d [km]')
ylabel('r(theta)-r0 [km]')
hold off

figure
hold on
title(['Differenz zwischen Euler und Runge-Kutta bei h = ',num2str(hi)])
for jj = 1:numel(phi0_vec)
    d_euler = (theta_span_euler(:,jj)).*RE;
    d_runge = (theta_span_runge(:,jj)).*RE;
    d_alter = abs(theta_span_euler(:,jj)-theta_span_runge(:,jj)).*180./pi.*RE;
    %Delta_theta = abs(theta_span_euler(:,jj)-theta_span_runge(:,jj));
    plot(r_span-RE,d_alter,'-')
    %plot(r_span-RE,Delta_theta,'*')
    xlabel('r [km]')
    ylabel('Differenz: Bezogen auf die horizontale Distanz d [km]')
    ylim([0 0.1])
end
legend(legendstring)
hold off

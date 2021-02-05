%%
c = 299792458.0; % Define speed of light (m/s)
aE = 6378137; % Define equitorial radius of Earth (m)
aP = 6356752.3142; % Define polar radius of Earth (m)
e = sqrt(1-aP^2/aE^2); % Eccentricity of Earth ellipsoid

%%
% Get data from files
fidNav = fopen('NavigationData.txt');
form = '%f %f %f %f %f %f %f';
rawNavData = textscan(fidNav, form,'headerlines', 2);
fclose(fidNav);
navData = cell2mat(rawNavData(1:6));
navData(:,2:4) = navData(:,2:4)*1000; % Convert satellite position from km to m
fidObs = fopen('ObservationData.txt');
form = '%f %f %f';
rawObsData = textscan(fidObs, form,'headerlines', 2);
fclose(fidObs);
obsData = cell2mat(rawObsData);

%%
f1 = 1575.42; % Define frequency of L1 satellite (MHz)
f2 = 1227.60; % Define frequency of L2 satellite (MHz)
gamma = (f1^2)/(f2^2);
rho_nought = (obsData(:,2)-gamma*obsData(:,3))/(1-gamma); % Calculate pseudo range measurement

%%
posR = [-5543600, -2054400, 2387500]; % Define initial guess of receiver position (m)
del_Cr = 0; % Define initial guess of receiver clock bias (s)
r_star = [posR del_Cr]'; % Define vector to be iterated
del_r_iter = 1; % Iteration counter
while true
    d = sqrt(sum((r_star(1:3)'-navData(:,2:4)).^2,2)); % Create initial guess of true range
    delT_i = navData(:,5); % Define satellite clock error
    rho_c = d+r_star(4)-c*delT_i; % Create expression of pseudorange
    % Create matrices for iteration:
    H = [(posR-navData(:,2:4))./d ones(length(obsData),1)];
    y = rho_nought-rho_c;
    
    % Least squares
    
    %del_r = (H'*H)\H'*y; % Perform least squares optimization
    
    %
    cvx_begin
        variable del_r(size(H,2))
        minimize( norm(H*del_r-y) )
    cvx_end
    
    r_star = r_star+del_r; % Improve guess
    if norm(del_r) < 10^-8 % Convergence condition
        break
    end
    if del_r_iter > 50 % Stop if position does not converge
        disp('Algorithm did not converge')
        return
    end
    del_r_iter = del_r_iter+1;
end
r_ECEF = r_star(1:3); % Define Earth center Earth fixed position of receiver (m)
rE = norm(r_ECEF); % Determine distance of receiver from center of Earth (m)
del_Tr = r_star(4)/c; % Calculate clock offset of receiver (s)
lambda = atan2d(r_ECEF(2),r_ECEF(1)); % Calculate longitude of observer (degrees)
% Calculate geocentric latitude of observer (degrees):
phi_geoc = atan2d(r_ECEF(3),sqrt(r_ECEF(1)^2+r_ECEF(2)^2));

% Calculate geodetic latitude of observer (degrees):
phi_old = phi_geoc; % Initial guess
phi_iter = 1; % Iteration counter
while true
    rN = aE/sqrt(1-e^2*sind(phi_old)^2); % Ellipsoid radius at position [m]
    % Ellipsoid height:
    if phi_old ~= 90 || phi_old ~= -90
        h = sqrt(r_ECEF(1)^2+r_ECEF(2)^2)/cosd(phi_old)-rN; % Ellipsoid height [m]
    else
        h = r_ECEF(3)/sind(phi_old)-rN+e^2*rN; 
    end
    % Geodetic latitude:
    phi = atand(r_ECEF(3)/sqrt(r_ECEF(1)^2+r_ECEF(2)^2)*(1-e^2*rN/(rN+h))^-1);
    if abs((phi-phi_old)/phi_old) < 10^-9 % Convergence condition
        break
    end
    if phi_iter > 50 % Stop if geodetic latitude does not converge
        disp('Algorithm did not converge')
        return
    end
    phi_old = phi;
    phi_iter = phi_iter+1;
end

% Display desired parameters:
disp(['Height = ', num2str(h),' m'])
disp(['Latitude = ', num2str(phi),' degrees'])
disp(['Longitude = ', num2str(lambda),' degrees'])
% disp(['Receiver Clock Offset=', num2str(del_Tr),' seconds'])

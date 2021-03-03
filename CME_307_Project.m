file_name = 'gnss_data.xlsx'; % File to read from
data_PR = read_GNSS_Data(file_name); % Read data and compute pseudoranges

t = unique(data_PR.ElapsedRealtimeMillis); % Times at which data is taken
n_old = 0; % Index in data
% Blank vectors to hold user position
userX = zeros(size(t));
userY = zeros(size(t));
userZ = zeros(size(t));
userb_u = zeros(size(t));
t_num = length(t);
locations = {'Chase Center','Oracle Park','Yerba Buena Gardens'};

for k = 1:t_num % Index in time vector
    
    y = [0;0;0;0]; % Arbitrary intial position and time
    t_indices = find(data_PR.ElapsedRealtimeMillis == t(k)); % Data points taken at t(k)
    n = max(t_indices); % Last index for current time
    [~,t_indices]=unique(data_PR.Svid(t_indices)); % Remove duplicates
    t_indices = t_indices+n_old;
    n_old = n;
    G = zeros(length(t_indices),4);
    
    while true
        
        r0 = sqrt((data_PR.X(t_indices)-y(1)).^2+(data_PR.Y(t_indices)-y(2)).^2 ...
            +(data_PR.Z(t_indices)-y(3)).^2); % Theoretical true range
        rho0 = r0-data_PR.B(t_indices)+y(4); % Theoretical pseudorange
        delta_rho = data_PR.rho_m(t_indices)-rho0;
        % Geometry matrix
        G(:,1) = -(data_PR.X(t_indices)-y(1))./r0;
        G(:,2) = -(data_PR.Y(t_indices)-y(2))./r0;
        G(:,3) = -(data_PR.Z(t_indices)-y(3))./r0;
        G(:,4) = ones(length(t_indices),1);
        
        cvx_begin quiet
            variable delta_y(size(G,2))
            minimize( norm(G*delta_y-delta_rho) )
        cvx_end
        
        if norm(delta_y(1:3)) < 10^-1
            break
        end
        
        y = y+delta_y;
        
    end
    
    userX(k) = y(1);
    userY(k) = y(2);
    userZ(k) = y(3);
    userb_u(k) = y(4);
    
    
    rE = 6378100;
    plotbound = 1.2*max(sqrt(data_PR.X(t_indices).^2 ...
        +data_PR.Y(t_indices).^2+data_PR.Z(t_indices).^2));
    
    h = figure;
    [x_Earth,y_Earth,z_Earth]=sphere(250);
    x_Earth = rE*x_Earth; y_Earth = rE*y_Earth; z_Earth = rE*z_Earth;
    s = surf(x_Earth,y_Earth,z_Earth,'edgecolor','none');
    alpha(0.5)
    pbaspect([1,1,1])
    xlim([-plotbound,plotbound])
    ylim([-plotbound,plotbound])
    zlim([-plotbound,plotbound])
    view_vec=[1,1,0.5];
    view(view_vec)
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'ZTick',[])
    title([locations(k),' Measurement'])
    hold on
    
    filename_gif = ['GPS_Sat_',num2str(k),'.gif'];
    filename_img = ['GPS_Sat_',num2str(k),'_Plot.png'];
    delete(filename_gif)
    for theta = 0:359
        
        R3=[cosd(theta), sind(theta), 0; -sind(theta), cosd(theta), 0; 0, 0, 1];
        userP = R3*[userX(k);userY(k);userZ(k)];
        satPos = R3*[data_PR.X(t_indices)';data_PR.Y(t_indices)';data_PR.Z(t_indices)'];
        
        userP_scatter = scatter3(userP(1),userP(2),userP(3),'r');
        satPos_scatter = scatter3(satPos(1,:),satPos(2,:),satPos(3,:),'k');
        
        drawnow
        % Capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if theta == 0
            imwrite(imind,cm,filename_gif,'gif','DelayTime',0.02,'Loopcount',inf);
        else
            imwrite(imind,cm,filename_gif,'gif','DelayTime',0.02,'WriteMode','append');
        end
        
        if theta == 240
            saveas(gcf,filename_img);
        end
        
        delete(userP_scatter);
        delete(satPos_scatter);
        
    end
    
end

wgs84 = wgs84Ellipsoid('meter');
[lat,lon,h] = ecef2geodetic(wgs84,userX,userY,userZ);

% Chase Center 37.767268, -122.385551
% Oracle Park 37.778034, -122.391543
% Yerba Buena Gardens 37.784897, -122.402524
d1 = distance(lat(1), lon(1), 37.767268, -122.385551, wgs84);
d2 = distance(lat(2), lon(2), 37.778034, -122.391543, wgs84);
d3 = distance(lat(3), lon(3), 37.784897, -122.402524, wgs84);


function data_PR = read_GNSS_Data(file_name)
    % Read GNSS data from specified file and compute the corresponding
    % pseudoranges in units of m
    
    c = 299792458; % Speed of light [m/s]
    sec_Wk = 604800; % Seconds in one week
    
    data_PR = readtable(file_name);
    data_PR.Nw = floor(-1e-9*data_PR.FullBiasNanos/sec_Wk); % Week number
    data_PR.tRxhw = data_PR.TimeNanos+data_PR.TimeOffsetNanos; % Corrected receiver time
    data_PR.bhw = data_PR.FullBiasNanos+data_PR.BiasNanos; % Adjusted bias in hardware clock
    data_PR.tRxGPS = data_PR.tRxhw-data_PR.bhw; % Receiver time in GPS frame
    data_PR.tRxw = data_PR.tRxGPS-data_PR.Nw*(sec_Wk*1e9); % Receiver time in current week frame
    data_PR.rho_ns = data_PR.tRxw-data_PR.ReceivedSvTimeNanos; % Pseudorange [ns]
    data_PR.rho_m = data_PR.rho_ns*(c*1e-9); % Pseudorange [m]
    
end
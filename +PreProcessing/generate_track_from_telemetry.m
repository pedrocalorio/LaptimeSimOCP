function Track = generate_track_from_telemetry(data,npoints,track_width)

% data is a matrix which was information from racing data
% npoints is the number of points the track will have

lateral_g = data(:,4)*9.806; %m/s²
speed = data(:,3)/3.6; % m/s

inv_radius = lateral_g./speed.^2;

old_sLap = linspace(data(1,2),data(end,2),length(data(:,2)));

% track distance with new points
Track.sLap = linspace(data(1,2),data(end,2),npoints);
% track curvature
Track.C = interp1(old_sLap,inv_radius,Track.sLap);
% vehicle speed for that lap
Track.Vxx = interp1(old_sLap,speed,Track.sLap);
% elapsed time 
Track.time = interp1(old_sLap,data(:,1),Track.sLap);
% track yaw angle with error
Track.psi = cumtrapz(Track.time,Track.C.*Track.Vxx);
% x and y coordinates 
Track.x = cumtrapz(Track.time,cos(Track.psi).*Track.Vxx);
Track.y = cumtrapz(Track.time,sin(Track.psi).*Track.Vxx);

%% Streching algorithm 
% used to make sure start and end point are the same
% Based on "On Minimum Time Vehicle Manoeuvring: The Theoretical Optimal
% Lap" D.Casanova

e_psi = abs(Track.psi(end)) - 2*pi;

delta_psi = -e_psi/Track.sLap(end);
Track.psi = cumtrapz(Track.time,(Track.C-delta_psi).*Track.Vxx);
Track.x = cumtrapz(Track.time,cos(Track.psi).*Track.Vxx);
Track.y = cumtrapz(Track.time,sin(Track.psi).*Track.Vxx);

delta_x = (Track.x(1)-Track.x(end))/Track.sLap(end);
delta_y = (Track.y(1)-Track.y(end))/Track.sLap(end);

%% Renaming some varaibles to output

Track.aYaw = Track.psi;
Track.curv = Track.C;
Track.d_sLap = [diff(Track.sLap), Track.sLap(end) - Track.sLap(end-1)];
Track.xCar = cumtrapz(Track.time,(cos(Track.psi)+delta_x).*Track.Vxx);
Track.yCar = cumtrapz(Track.time,(sin(Track.psi)+delta_y).*Track.Vxx);
Track.TrackWidth = track_width;

Track.xCarLeft = Track.xCar + track_width/2*gradient(Track.yCar)./(gradient(Track.xCar).^2 + gradient(Track.yCar).^2).^0.5;
Track.yCarLeft = Track.yCar - track_width/2*gradient(Track.xCar)./(gradient(Track.xCar).^2 + gradient(Track.yCar).^2).^0.5;
Track.xCarRight = Track.xCar - track_width/2*gradient(Track.yCar)./(gradient(Track.xCar).^2 + gradient(Track.yCar).^2).^0.5;
Track.yCarRight = Track.yCar + track_width/2*gradient(Track.xCar)./(gradient(Track.xCar).^2 + gradient(Track.yCar).^2).^0.5;


end
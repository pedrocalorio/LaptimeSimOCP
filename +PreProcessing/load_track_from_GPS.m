function Track = load_track_from_GPS(track_string,nGrid)

% reads the data from the sheet 
data = readmatrix(track_string);

% higher no. of segments causes trajectory to follow the reference line
nseg = length(data(:,1));

% x, y and track width data 
x   = data(:,1);
y   = data(:,2);
twr = data(:,3);
twl = data(:,4);

pathXY     = [x y];
boundaries = [twr twl];

stepLengths   = sqrt(sum(diff(pathXY,[],1).^2,2));
stepLengths   = [0; stepLengths]; % add the starting point
cumulativeLen = cumsum(stepLengths);
finalStepLocs = linspace(0,cumulativeLen(end), nseg);

finalPathXY          = interp1(cumulativeLen, pathXY, finalStepLocs);
finalTrackBoundaries = interp1(cumulativeLen, boundaries, finalStepLocs);

xt   = finalPathXY(:,1);
yt   = finalPathXY(:,2);
twrt = finalTrackBoundaries(:,1);
twlt = finalTrackBoundaries(:,2);

ind_step_preview_psi  = round(1 /  (mean(stepLengths)));
ind_step_review_psi   = round(1 / (mean(stepLengths)));
ind_step_preview_curv = round(2 / (mean(stepLengths)));
ind_step_review_curv  = round(2 / (mean(stepLengths)));

ind_step_preview_psi  = max(ind_step_preview_psi, 1);
ind_step_review_psi   = max(ind_step_review_psi, 1);
ind_step_preview_curv = max(ind_step_preview_curv, 1);
ind_step_review_curv  = max(ind_step_review_curv, 1);

steps_tot_psi  = ind_step_preview_psi + ind_step_review_psi;
steps_tot_curv = ind_step_preview_curv + ind_step_review_curv;

path_temp = [pathXY(end-ind_step_review_psi+1,:); pathXY; pathXY(1:ind_step_preview_psi,:)];

x_temp = path_temp(:,1);
y_temp = path_temp(:,2);

tangvecs = [ path_temp(1+steps_tot_psi:end, 1)' - path_temp(1:length(x_temp)-steps_tot_psi,1)';
                     path_temp(1+steps_tot_psi:end, 2)' - path_temp(1:length(x_temp)-steps_tot_psi,2)' ]';
% tangvecs = [ x_temp(1+steps_tot_psi) - x_temp(end-steps_tot_psi-1),
%                      path_temp(1+steps_tot_psi, 2) - path_temp(end-steps_tot_psi-1,2) ];

% # calculate psi of tangent vectors (pi/2 must be substracted due to our convention that psi = 0 is north)
psi = atan2(tangvecs(:, 2), tangvecs(:, 1)) - pi/ 2 ;
psi = psi - psi(1);
psi = +PreProcessing.normalize_psi(psi);

%% starts calculation of curvature

psi_temp = [psi; 0; psi(end-ind_step_review_psi+1,1)];
psi_temp(end) =  psi_temp(1:ind_step_preview_psi);

delta_psi = zeros(length(x),1);

for i=1:length(x)
    delta_psi(i) = +PreProcessing.normalize_psi(psi_temp(i + steps_tot_curv) - psi_temp(i));
end

s_points_cl = [cumulativeLen; 0; 0.0];
s_points = s_points_cl(1:length(s_points_cl)-1,1);

s_points_cl_reverse = flipud(-cumsum(flipud(stepLengths)));

s_points_temp = [s_points; s_points_cl_reverse(end-ind_step_review_curv,1) ];
s_points_temp(end) = s_points_cl(end) + s_points(1:ind_step_preview_curv) ;

kappa = delta_psi ./ (s_points_temp(steps_tot_curv+1:end) - s_points_temp(1:length(s_points_temp)-steps_tot_psi,1));
kappa = smooth(cumulativeLen,kappa,16,'lowess');

% normal direction for each vertex
dx = gradient(xt);
dy = gradient(yt);
dL = hypot(dx,dy);

ddx = gradient(dx);
ddy = gradient(dy);

% offset curve - anonymous function
xoff = @(a) -a*dy./dL + xt;
yoff = @(a)  a*dx./dL + yt;

% offset data
offset = [-twlt twrt];
for i = 1:numel(xt)
    xin = xoff(offset(i,1));      % get inner offset curve
    yin = yoff(offset(i,1));
    
    xout  = xoff(offset(i,2));      % get outer offset curve
    yout  = yoff(offset(i,2));
end

%% Alternative way of calculating the curvature

% C_num = (dx.*ddy - ddx.*dy);
% C_den = (sqrt(dx.^2+dy.^2)).^3;
% 
% kappa_v2 = C_num./C_den;


Track = struct();
Track.distance = cumulativeLen;
Track.sLap = linspace(0,cumulativeLen(end),nGrid);
Track.aYaw = unwrap_signal(psi);
Track.curv = kappa;
Track.d_sLap = [diff(Track.sLap), Track.sLap(end) - Track.sLap(end-1)];
Track.xCar = xt;
Track.yCar = yt;
Track.TrackWidth = 1.00*mean(twr+twl);
Track.xCarLeft = xin;
Track.yCarLeft = yin;
Track.xCarRight = xout;
Track.yCarRight = yout;
Track.left_offset = offset(:,1);
Track.right_offset = offset(:,2);


% Track.aYaw = interp1(cumulativeLen,Track.aYaw,Track.sLap,'spline');
% Track.curv = interp1(cumulativeLen,Track.curv,Track.sLap,'spline');
% Track.d_sLap = interp1(cumulativeLen,Track.d_sLap,Track.sLap,'spline');
% Track.xCar = interp1(cumulativeLen,xt,Track.sLap,'spline');
% Track.yCar = interp1(cumulativeLen,yt,Track.sLap,'spline');
% Track.xCarLeft  = interp1(cumulativeLen,xin,Track.sLap,'spline');
% Track.yCarLeft  = interp1(cumulativeLen,yin,Track.sLap,'spline');
% Track.xCarRight = interp1(cumulativeLen,xout,Track.sLap,'spline');
% Track.yCarRight = interp1(cumulativeLen,yout,Track.sLap,'spline');
% Track.left_offset = interp1(cumulativeLen,offset(:,1),Track.sLap,'spline');
% Track.right_offset = interp1(cumulativeLen,offset(:,2),Track.sLap,'spline');

% % plot starting line
% plot([xin(1) xout(1)], [yin(1) yout(1)],'color','b','linew',2)
% % plot([xin(2) xout(2)], [yin(2) yout(2)],'color','k','linew',2)
% 
% % plot reference line
% plot(xt,yt,'--')
% hold on
% 
% % plot inner track
% plot(xin,yin,'color','k')
% 
% % plot outer track
% plot(xout,yout,'color','k')
% hold off
% axis equal


end

function signal = unwrap_signal(signal)


for i=2:length(signal)
    difference = signal(i)-signal(i-1);
    if difference > pi
        signal(i:end) = signal(i:end) - 2*pi;
    elseif difference < -pi
        signal(i:end) = signal(i:end) + 2*pi; 
    end
end

end
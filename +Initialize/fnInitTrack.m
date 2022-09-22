function [problem] = fnInitTrack(problem, bPlotTrackDef, generateTrack)

if generateTrack

    problem.dsSystem.generateTrack = generateTrack;
    
    %% Generate the track definition

    switch problem.options.method
        case 'trapezoid'
            nGrid = problem.options.trapezoid.nGrid;
        case 'hermiteSimpson'
            nGrid = 2*problem.options.hermiteSimpson.nSegment + 1;
        otherwise
            error('Invalid method.');
    end

%     td = struct();
%     
%     TrackWidth = 8;
%     td.TrackWidth = TrackWidth;
%     td.sLap = linspace(0,400,nGrid);
%     td.curv = zeros(1,nGrid);
%     td.curv( (td.sLap > median(td.sLap)) & (td.sLap < (median(td.sLap) + pi*15/1) )) = 1/20;
%     td.curv_nosmooth = td.curv; % Store
%     td.curv = smooth(td.curv, 0.3, 'lowess')';
% 
%     td.d_sLap   = [diff(td.sLap), td.sLap(end) - td.sLap(end-1)];
%     td.aYaw  	= cumsum(td.curv.*td.d_sLap);
%     td.xCar 	= cumsum(td.d_sLap.*cos(-td.aYaw));
%     td.yCar 	= cumsum(td.d_sLap.*sin(-td.aYaw));
%     
%     td.xCarLeft = td.xCar + TrackWidth/2*gradient(td.yCar)./(gradient(td.xCar).^2 + gradient(td.yCar).^2).^0.5;
%     td.yCarLeft = td.yCar - TrackWidth/2*gradient(td.xCar)./(gradient(td.xCar).^2 + gradient(td.yCar).^2).^0.5;
%     td.xCarRight = td.xCar - TrackWidth/2*gradient(td.yCar)./(gradient(td.xCar).^2 + gradient(td.yCar).^2).^0.5;
%     td.yCarRight = td.yCar + TrackWidth/2*gradient(td.xCar)./(gradient(td.xCar).^2 + gradient(td.yCar).^2).^0.5;
    
    td = struct();

    td.Distance = 0;
    lStraight = 200;

    RCorner = 50;    
    RCorner2 = 50;

    lCorner  = pi*RCorner*2;
    lCorner2 = pi*RCorner2*2;
    TrackWidth = 10;
    
    %%
    td.Curvature(1:(lStraight/2)) = 0;
    td.Curvature((lStraight/2+1):((lStraight/2) + lCorner)) = (1/RCorner)*(1-cos((1:lCorner)/lCorner*pi*2))/2;
    td.Curvature(((lStraight/2) + lCorner + 1):((lStraight/2) + lCorner + lStraight)) = 0;
    td.Curvature(((lStraight/2) + lCorner + lStraight + 1):((lStraight/2) + lCorner + lCorner2 + lStraight)) = (1/RCorner2)*(1-cos((1:lCorner2)/lCorner2*pi*2))/2;
    td.Curvature(((lStraight/2) + lCorner + lCorner2 + lStraight + 1):(lCorner + lCorner2 + 2*lStraight)) = 0;
    td.Curvature(end+1) = td.Curvature(1);
    %%
    % td.Curvature = simout;
    td.sLap = linspace(0, (length(td.Curvature)-1), nGrid);
    %  td.Distance = linspace(0, 400, 400)';
    td.Curvature = interp1(0 : (length(td.Curvature)-1), td.Curvature, td.sLap);    
    td.curv = td.Curvature;
    % td.Curvature = interp1(Dist, td.Curvature, td.Distance);
    td.aYaw = cumsum(td.Curvature.*gradient(td.sLap));
    td.dAngle = gradient(td.aYaw);
    td.xCar = cumsum(cos(td.aYaw).*gradient(td.sLap));
    td.yCar = cumsum(sin(td.aYaw).*gradient(td.sLap));
    td.xCarLeft = td.xCar + TrackWidth/2*gradient(td.yCar)./(gradient(td.xCar).^2 + gradient(td.yCar).^2).^0.5;
    td.yCarLeft = td.yCar - TrackWidth/2*gradient(td.xCar)./(gradient(td.xCar).^2 + gradient(td.yCar).^2).^0.5;
    td.xCarRight = td.xCar - TrackWidth/2*gradient(td.yCar)./(gradient(td.xCar).^2 + gradient(td.yCar).^2).^0.5;
    td.yCarRight = td.yCar + TrackWidth/2*gradient(td.xCar)./(gradient(td.xCar).^2 + gradient(td.yCar).^2).^0.5;

    td.TrackWidth = TrackWidth;

    problem.dsSystem.td = td;
    
else
    
    problem.dsSystem.generateTrack = generateTrack;
    
    switch problem.options.method
        case 'trapezoid'
            nGrid = problem.options.trapezoid.nGrid;
        case 'hermiteSimpson'
            nGrid = 2*problem.options.hermiteSimpson.nSegment + 1;
        otherwise
            error('Invalid method.');
    end
    
%     % Import track from data
%     trackData = readtable('Nurburgring');
%     
%     Track = struct();
%     td = struct();
% 
%     nPoints = 1500;
%     Track.sLap = trackData.CorrLapDist';
% 
%     Track.TrackWidth = 8;
%     
%     Track.sLap = linspace(0,Track.sLap(end),nPoints);
%     Track.curv = interp1(trackData.CorrLapDist',...
%                       trackData.MathInverseCornerRadius',...
%                       Track.sLap,'linear','extrap');
%                   
%     Track.curv_nosmooth = Track.curv;
%     
%     Track.d_sLap   = [diff(Track.sLap), Track.sLap(end) - Track.sLap(end-1)];
%     Track.aYaw  	= cumsum(Track.curv.*Track.d_sLap);
%     td.xCar 	= cumsum(Track.d_sLap.*cos(Track.aYaw));
%     td.yCar 	= cumsum(Track.d_sLap.*sin(Track.aYaw));
% 
%     % here we have the track centerline with initial and end not together
% 
%     % calculate delta y and delta x
%     deltaX = (td.xCar(1) - td.xCar(end))/nPoints;
%     deltaY = (td.yCar(1) - td.yCar(end))/nPoints;
%     % put these delta x and y into the "new" x and y position
% %     td.xCar = td.xCar - deltaX;
% %     td.yCar = td.yCar - deltaY;
%     td.xCar = cumsum(td.xCar-td.xCar*deltaX);
%     td.yCar = cumsum(td.yCar-td.yCar*deltaY);
% 
%     % calculates the delta heading angle
%     newaYaw = atan(deltaY/deltaX)/nPoints;
% 
%     % updates it in the "old" heading angle
%     Track.Angle = Track.aYaw - newaYaw;
%     % then we have the final x and y coordinate
%     Track.XCoord 	= cumsum(Track.d_sLap.*cos(Track.Angle));
%     Track.YCoord 	= cumsum(Track.d_sLap.*sin(Track.Angle));
% 
%     Track.XCoordLeft = Track.XCoord + Track.TrackWidth/2*gradient(Track.YCoord)./(gradient(Track.XCoord).^2 + gradient(Track.YCoord).^2).^0.5;
%     Track.YCoordLeft = Track.YCoord - Track.TrackWidth/2*gradient(Track.XCoord)./(gradient(Track.XCoord).^2 + gradient(Track.YCoord).^2).^0.5;
%     Track.XCoordRight = Track.XCoord - Track.TrackWidth/2*gradient(Track.YCoord)./(gradient(Track.XCoord).^2 + gradient(Track.YCoord).^2).^0.5;
%     Track.YCoordRight = Track.YCoord + Track.TrackWidth/2*gradient(Track.XCoord)./(gradient(Track.XCoord).^2 + gradient(Track.YCoord).^2).^0.5;
% 
%     td.TrackWidth = Track.TrackWidth;
%     td.sLap = linspace(Track.sLap(1),Track.sLap(end),nGrid);
%     td.curv = interp1(Track.sLap,Track.curv,td.sLap);
%     td.curv_nonsmooth = td.curv;
%     td.d_sLap   = [diff(td.sLap), td.sLap(end) - td.sLap(end-1)];
%     td.aYaw  	= interp1(Track.sLap,Track.Angle,td.sLap);
%     td.xCar 	= interp1(Track.sLap,Track.XCoord,td.sLap);
%     td.yCar 	= interp1(Track.sLap,Track.YCoord,td.sLap);
%     td.xCarLeft = interp1(Track.sLap,Track.XCoordLeft,td.sLap);
%     td.yCarLeft = interp1(Track.sLap,Track.YCoordLeft,td.sLap);
%     td.xCarRight = interp1(Track.sLap,Track.XCoordRight,td.sLap);
%     td.yCarRight = interp1(Track.sLap,Track.YCoordRight,td.sLap); 

    % BARCELONA

    load('BAR')
    
    
    td = struct();
    
    td.TrackWidth = Track.TrackWidth;
    
    td.sLap = linspace(Track.sLap(1),Track.sLap(end),nGrid);
    td.curv = interp1(Track.sLap,Track.rCurv,td.sLap);
    td.curv_nonsmooth = td.curv;
    td.d_sLap   = [diff(td.sLap), td.sLap(end) - td.sLap(end-1)];
    td.aYaw  	= interp1(Track.sLap,Track.Angle,td.sLap);
    td.xCar 	= interp1(Track.sLap,Track.XCoord,td.sLap);
    td.yCar 	= interp1(Track.sLap,Track.YCoord,td.sLap);
    td.xCarLeft = interp1(Track.sLap,Track.XCoordLeft,td.sLap);
    td.yCarLeft = interp1(Track.sLap,Track.YCoordLeft,td.sLap);
    td.xCarRight = interp1(Track.sLap,Track.XCoordRight,td.sLap);
    td.yCarRight = interp1(Track.sLap,Track.YCoordRight,td.sLap); 


    
    problem.dsSystem.td = td;
end

%% Plot the track definition

if bPlotTrackDef
    fnPlotTrackDefinition(td);
end

end

function fnPlotTrackDefinition(td)

f1 = figure(1);
clf(f1);

f1.Name = 'Track Definition';

figure;plotbrowser;
plot(td.sLap, td.curv);
hold on;
plot(td.sLap, td.curv, 'LineWidth', 2);
grid minor
% plot(td.sLap, td.curv, 'ob');
title('Track Curvature');
xlabel('Center-Line Distance [m]'); ylabel('Curvature [1/m]');
set(gca,'FontSize',15)


figure;plotbrowser;
plot(td.sLap, td.aYaw);
title('Vehicle yaw angle');
xlabel('Center-Line Distance [m]'); ylabel('\theta [deg]');
set(gca,'FontSize',15)

figure;plotbrowser;
plot(td.sLap, td.xCar, 'b');
hold on;
plot(td.sLap, td.yCar, 'r');
title('X and Y-Coordinate against Center-Line Distance');
legend('X', 'Y');
xlabel('Center-Line Distance [m]'); ylabel('X and Y-Coordinate');
set(gca,'FontSize',15)

figure;plotbrowser;
plot(td.xCar, td.yCar, 'black--',td.xCarLeft, td.yCarLeft,'black',td.xCarRight, td.yCarRight,'black' );
hold on
plot(0,0,'>','LineWidth',3)
% xlim([min(Track.xCar)-50 max(Track.xCar)+50])
% title('Circuit de Barcelona-Catalunya','Interpreter','latex');
xlabel('X-Coordinate [m]','Interpreter','latex'); ylabel('Y-Coordinate [m]','Interpreter','latex');
daspect([1.5 1.15 10])
set(gca,'FontSize',15)
legend('Center-Line','Left Bound','Right Bound','Starting','Interpreter','latex','location','best')
% subplot(2,2,4)
% plot(td.xCar, td.yCar);
% title('Car position');
% xlabel('xCar'); ylabel('yCar');

% axis equal

end
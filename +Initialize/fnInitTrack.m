function [problem] = fnInitTrack(problem, bPlotTrackDef, generateTrack)

if generateTrack

    % just stores the boolean for creating the oval track from math in the problem struct.
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

    % Here a simple oval track is created from geometry 
    
    td = struct();

    td.Distance = 0;
    % set the length of the two straights 
    lStraight = 200;

    % defines the corner radius of the oval circuit 
    RCorner = 50;    
    RCorner2 = 50;

    lCorner  = pi*RCorner*2;
    lCorner2 = pi*RCorner2*2;
    TrackWidth = 10;
    
    % here is where the processing to create the track is done, we start
    % from the main straight and finishes with the same straight but after
    % passing through all the circuit 

    td.Curvature(1:(lStraight/2)) = 0;
    td.Curvature((lStraight/2+1):((lStraight/2) + lCorner)) = (1/RCorner)*(1-cos((1:lCorner)/lCorner*pi*2))/2;
    td.Curvature(((lStraight/2) + lCorner + 1):((lStraight/2) + lCorner + lStraight)) = 0;
    td.Curvature(((lStraight/2) + lCorner + lStraight + 1):((lStraight/2) + lCorner + lCorner2 + lStraight)) = (1/RCorner2)*(1-cos((1:lCorner2)/lCorner2*pi*2))/2;
    td.Curvature(((lStraight/2) + lCorner + lCorner2 + lStraight + 1):(lCorner + lCorner2 + 2*lStraight)) = 0;
    td.Curvature(end+1) = td.Curvature(1);

    td.sLap = linspace(0, (length(td.Curvature)-1), nGrid);
    td.Curvature = interp1(0 : (length(td.Curvature)-1), td.Curvature, td.sLap);    
    td.curv = td.Curvature;
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
    
    % just stores the boolean for re-constructing the track from data in
    % the problem struct.
    problem.dsSystem.generateTrack = generateTrack;
    

    % gets the number of discretization points as choosen previosly outside
    % this function
    switch problem.options.method
        case 'trapezoid'
            nGrid = problem.options.trapezoid.nGrid;
        case 'hermiteSimpson'
            nGrid = 2*problem.options.hermiteSimpson.nSegment + 1;
        otherwise
            error('Invalid method.');
    end

    % reads the telemetry data used to re-construct the centerline
    data = readmatrix("+TelemetryData\PaulRicard.csv");
    
   
    % define the track width
    track_width = 10; %m
    
    % calls the function that generates the track
%     Track = PreProcessing.generate_track_from_telemetry(data,nGrid,track_width);
    Track = PreProcessing.load_track_from_GPS("+TrackGPSData/Sakhir.csv",nGrid);

    
    problem.dsSystem.td = Track;
end

%% Plot the track definition

if bPlotTrackDef
    fnPlotTrackDefinition(problem.dsSystem.td);
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
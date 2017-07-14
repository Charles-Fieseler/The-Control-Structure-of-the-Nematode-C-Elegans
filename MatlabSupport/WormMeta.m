function [ dat ] = WormMeta( filename, settings )
%% Worm Metadata
% Given a .csv file of the worm's movement, returns overall
% statistics for the movement.
%
% INPUTS
%   filename - csv file of the worms movement
%   settings - struct that contains the visualization settings.
%       verbose = false; %Displays extra information
%       toPlot = true; %Whether or not to plot the information (regardless, returns it in a struct)
%       modulo2pi = true; %If calculating angles, whether or not to apply mod(2pi) to them
%       moduloLower = -pi/2; %The lower range of the mod(2pi) turning angle
%
% OUTPUTS - 
%   dat - struct which contains all the calculated metadata and the
%       user-passed settings this was run with. The fields are:
%
%       .tspan - an array of the times
%
%       .barVelocitiesRaw - velocities of the body segments; each side
%       saved separately (Ventral and Dorsal)
%       .barVelocities - above, averaged over side
%       .barAngVelocities - angular velocities of above
%       .barVelocitiesMax - maximum body segment velocity
%       .barVelocitiesMean - mean body segment velocity
%
%       .COMVelocityRaw - SAME AS ABOVE but for center of mass
%       .COMVelocity
%       .COMAngVelocity
%       .COMVelocityMax
%       .COMVelocityMean
%
%       .barKinetic - (linear) kinetic energy of the body segments
%       .barAngKinetic - (angular) kinetic energy of the body segments
%       .totKinetic - sum of above
%
%       .steadyStateFit - straight line fit to the last 3 seconds of motion
%       .steadyStateX - mean COM X positions that were used for the fit
%       .steadyStateY - mean COM Y positions that were used for the fit
%
%       .originalAngleRaw - angle that the worm turned, as compared with
%       the heading for default settings
%       .originalAngle - above, but modulo 2pi (if setting is checked;
%       otherwise, is the same)
%
%
% EXAMPLES
%
%   %%Uses data from a simulation with the original settings
%   dat = WormMeta('../Model/simdata_original.csv')
%
%
% Dependencies
%   Other m-files required: (updated on 10-Jul-2017)    
%         MATLAB (version 9.1)
%         v2struct.m
%         simdata_original_metadata.mat
%
%   See also: WormMetaDataObj.m
%
%
% Author: Charles Fieseler
% University of Washington, Dept. of Physics
% Email address: charles.fieseler@gmail.com
% Website: coming soon
% Created: 10-Jul-2017
%========================================


%% Initialize with defaults

defaults = struct; %The default values
defaults.verbose = false; %Displays extra information
defaults.toPlot = true; %Whether or not to plot the information (regardless, returns it in a struct)
defaults.modulo2pi = true; %If calculating angles, whether or not to apply mod(2pi) to them
defaults.moduloLower = -pi/2; %The lower range of the mod(2pi) turning angle


if ~exist('settings','var')
    settings = struct;
end

namesS = fieldnames(settings);
namesD = fieldnames(defaults);

for j = 1:length(namesS)
    n = namesS{j};
    
    if max(strcmp(n,namesD)) > 0 %Check to see if the given setting is used
        defaults.(n) = settings.(n);
    else
        fprintf('Warning: "%s" setting not used\n',n)
    end
end

[verbose, toPlot, modulo2pi, moduloLower] ...
    = v2struct(defaults); %Unpacks the struct into variables


%% Import worm data from *.csv file
if isempty(strfind(filename,'.'))
    filename = [filename '.csv'];
end

if isempty(strfind(filename,'.mat'))
    
    if exist(['../Simdata/' filename],'file')
        if verbose
            disp('Importing data from ../Simdata')
        end
        
        data = importdata(['../Simdata/' filename]);
        
    elseif exist(['../Model/' filename],'file')
        if verbose
            disp('Importing data from ../Model')
        end
        
        data = importdata(['../Model/' filename]);
        
    elseif exist(['./' filename],'file')
        if verbose
            disp('Importing data from the current directory')
        end
        
        data = importdata(['./' filename]);
        
    end
else
    data = importdata(filename); %This is then a struct; the user should input the full filename
    data = data.yval.';
    data = [ zeros(size(data,1),1), data ];
end
if size(data,1) <= 1
    error('%i frame(s) stored in the data file (need more than 1).',size(data,1));
end

%==========================================================================


%% Initialize variables and get COM

Sz = size(data);
Nt = round(Sz(1));
Dt = data(2,1) - data(1,1);
Duration = data(end,1);
%Save the time span
dat.tspan = data(:,1);

Nbar = (Sz(2)-1)/3;
NSEG = Nbar-1;
D = 80e-6;

CoM = zeros(Nt,Nbar,3);

for i = 1:Nt
    frame = i;
    for j = 1:Nbar %Extract COM data from a single long column vector per frame
        CoM(i,j,1) = data(frame,1 + (j-1)*3 + 1);
        CoM(i,j,2) = data(frame,1 + (j-1)*3 + 2);
        CoM(i,j,3) = data(frame,1 + (j-1)*3 + 3);
    end
end


%==========================================================================


%% Extract velocity from COM

%Get the differential distances, squared
dx = diff(CoM(:,:,1),1);
dx2 = dx.^2;
dy = diff(CoM(:,:,2),1);
dy2 = dy.^2;
dtheta = diff(CoM(:,:,3),1);

%Add, average over bars, and divide by dt to get the velocity
rawBarVels = sqrt(dx2 + dy2)./Dt;
barVels = sign(dx).*rawBarVels;
barAngVels = dtheta./Dt;

%Save
dat.barVelocitiesRaw = rawBarVels;
dat.barVelocities = mean(barVels,2);
dat.barAngVelocities = mean(barAngVels,2);

%Also save the max and mean
dat.barVelocitiesMax = max(abs(dat.barVelocities));
dat.barVelocitiesMean = mean(abs(dat.barVelocities));


%Next, get the overall CoM velocity
dx = diff( mean(CoM(:,:,1),2) );
dx2 = dx.^2;
dy = diff( mean(CoM(:,:,2),2) );
dy2 = dy.^2;
dtheta2 = diff( mean(CoM(:,:,3),2) );

%Add, divide by Dt, save
rawVelCoM = sqrt(dx2 + dy2)./Dt;
velCoM = sign(dx).*rawVelCoM;
angvelCoM = dtheta2./Dt;
dat.COMVelocityRaw = rawVelCoM;
dat.COMVelocity = velCoM;
dat.COMAngVelocity = angvelCoM;

%Also save the max and mean
dat.COMVelocityMax = max(abs(velCoM));
dat.COMVelocityMean = mean(abs(velCoM));

%==========================================================================


%% Calculate Kinetic energy of body
% This is a rough calculation just using the speeds of the segments. A more
% realistic calculation would involve integrating the forces produced by
% the muscles over time, as the passive body forces are a large part of the
% energy expenditure of the worm

%Calculate moments of inertia of the segments, approximating them as
%cylinders of different radii
numSegs = 49;
%Mass of a worm segment
mSeg = 1e-6/numSegs; %units: micro grams; from http://bmcecol.biomedcentral.com/articles/10.1186/1472-6785-9-14

%segs = 1:numSegs;
%R = 40e-6; %Max radius; originally micrometers
%bodyRadii = R*sin(acos( (segs-numSegs/2-1)/(numSegs/2+0.2) ));
lSeg = (1e-3)/numSegs;
inertia = (mSeg*lSeg.^2)/12;

%Linear and angular kinetic energy
barKinetic = (1/2)*sum(mSeg*barVels.^2,2);
barAngKinetic = (1/2)*sum(inertia.*barAngVels.^2,2);
totKinetic = barKinetic + barAngKinetic;

%Save
dat.barKinetic = barKinetic;
dat.barAngKinetic = barAngKinetic;
dat.totKinetic = totKinetic;

%==========================================================================


%% Plot, if the option is 'true'

if toPlot
    %Plot mean speeds
    figure;
    subplot(2,1,1)
    plot(dat.tspan(1:end-1),dat.barVelocities);
    xlabel('time (s)');
    ylabel('velocity (m/s)');
    legend({filename},'Interpreter','none');
    title('Average velocity of each body segment')
    
    subplot(2,1,2)
    plot(dat.tspan(1:end-1),dat.COMVelocity);
    xlabel('time (s)');
    ylabel('velocity (m/s)');
    legend({filename},'Interpreter','none');
    title('Average velocity of the overall COM')
    
    %Plot kinetic energy
    %     figure;
    %     plot(dat.tspan(1:end-1),dat.totKinetic)
    %     xlabel('time (s)')
    %     ylabel('Kinetic energy (J)')
    %     title('Kinetic energy of the worm')
end
%==========================================================================


%% Fit the last 3 seconds to a line

%We want to calculate some things in relation to the original
%simulation statistics, so let's load that data if it exists
if ~exist('simdata_original_metadata.mat','file') && ...
        ~strcmp(filename,'../Model/simdata_original.csv')
    
    orig_dat = WormMeta('../Model/simdata_original.csv');
    save('simdata_original_metadata.mat','orig_dat');
    
elseif ~strcmp(filename,'../Model/simdata_original.csv')
    load simdata_original_metadata.mat
end

dat = fitSteadyState(dat, CoM);

%==========================================================================


%% Private functions

    function dat = fitSteadyState(dat, CoM)
        %This function calculates the steady state velocity and the times
        %that the worm is in that state. Note that the worm is always
        %undulating, so there is that "noise" in the velocity data
        
        %Assume the steady state lasts for at least the last 3 seconds of
        %the data
        SSduration = 3.0;
        SSwhich = dat.tspan>dat.tspan(end)-SSduration;
        %         meanVel = mean(dat.COMVelocityRaw(SSwhich));
        %Get the max deviation in the assumed steady state
        %         maxDiff = max(dat.COMVelocityRaw(SSwhich)-meanVel);
        
        meanCOM = [mean(CoM(:,:,1),2), mean(CoM(:,:,2),2)];
        xVals = meanCOM(SSwhich,1);
        yVals = meanCOM(SSwhich,2);
        [linFit, linErr] = polyfit(xVals,yVals,1);
        %However, the data might be nearly vertical, in which case we
        %should switch X and Y and then 
        [switchFit, switchErr] = polyfit(yVals,xVals,1);
        if (linErr.normr/switchErr.normr)>1
            %If the switched fit error is better, then invert the fit back
            %to the normal style of x as the independent variable
            linFit = [1 -switchFit(2)]./switchFit(1);
        end
        
        dat.steadyStateFit = linFit;
        dat.steadyStateX = meanCOM(SSwhich,1);
        dat.steadyStateY = meanCOM(SSwhich,2);
        
        if toPlot
            figure
            plot(meanCOM(SSwhich,1),polyval(myFit,meanCOM(SSwhich,1)))
        end

        
        %Compare the final angles of the steady state trajectories
        %   Different if we are processing the original data file, of
        %   course
        if ~strcmp(filename,'../Model/simdata_original.csv')
            myDirection = dat.steadyStateX(end)-dat.steadyStateX(1);
            changeInAngle = -orig_dat.originalAngle + ...
                atan(linFit(1)) + ...
                double(myDirection<0).*pi;
            %Want the range to be -pi/2 to 3pi/2 by default
            dat.changeInAngleRaw = changeInAngle;
            if modulo2pi
                dat.changeInAngle = ...
                    mod(changeInAngle-moduloLower,2*pi)+moduloLower;
            else
                dat.changeInAngle = changeInAngle;
            end
            if verbose
                fprintf('The change in angle is: %f\n',dat.changeInAngle)
            end
        else
            %If the worm is going to the left, then we add pi
            myDirection = dat.steadyStateX(end)-dat.steadyStateX(1);
            originalAngle = atan(linFit(1)) + ...
                double(myDirection<0).*pi;
            dat.originalAngleRaw = originalAngle;
            dat.originalAngle = mod(originalAngle+pi,2*pi)-pi;
        end
        
    end
%==========================================================================

end


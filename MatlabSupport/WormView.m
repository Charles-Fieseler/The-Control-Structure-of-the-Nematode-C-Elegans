function [dat] = WormView(filename, settings)
%% Program to visualise the output of Worm.cpp
%Single-file visualization program for c elegans simulations
%
%
% INPUTS
%   filename - data file to use for visualization
%   settings - struct that contains the visualization settings. Field names
%              should correspond to the property being set:
%         makeMovie = false; %Exports a movie (.mov)
%         verbose = false; %Displays extra information about the program execution
%         usingLegend = false; %Displays a legend on the figures
%         usingObjects = false; %Sometimes there are obstacles that the worm might hit; will be in a separate file
%         usingCOM = false; %To display Center Of Mass
%         usingFullWorm = true;
%         usingMetaDataFunction = false; %Metadata regarding stimuli might be saved
%         usingMetaDataExport = false; %Save the metadata figures
%         usingMetaDataNormalize = false; %If a negative wave-like stimulus, flip it and make the baseline 0
%         usingSpline = false; %This is an attempt at fitting the long-term trajectory of the animal
%         pauseLength = 0.0; %Duration to pause after every timestep
%         pauseAt = 0.0; %Pause here until user presses "enter"
%         startAtPause = false; %Start at the above pause value
%         quitAtPause = false; %Quit at the above pause value
%         isFolder = false; %Can pass a folder to this function; in that case it recursively calls itself for each file
%         downsample = 1; %Speed up visualization by throwing out data
%
% OUTPUTS - 
%   dat - data struct that contains the figures produced, the metadata (if
%   used), and the settings struct passed by the user
%
% EXAMPLES:
%
%   %%Normal forward motion
%   WormView('../Simdata_N_SR_2nd/simdata_N_SR_6_Med_0.9.csv')
%
%   %%Simple omega turn
%   WormView('../Simdata_OmWave_Jac/simdata_Jac_K_PE1.00.csv')
%
%
% Dependencies
%   Other m-files required: (updated on 10-Jul-2017)    
%         MATLAB (version 9.1)
%         v2struct.m
%         WormMeta.m
%         simdata_original_metadata.mat
%
%   See also: Based on the original visualization function from the paper:
%   Boyle, Jordan H., Stefano Berri, and Netta Cohen. 
%   “Gait Modulation in C. Elegans: An Integrated Neuromechanical Model.” 
%   Frontiers in Computational Neuroscience 6 (2012): 10. PMC. Web. 10 July 2017.
%
%
% Author: Charles Fieseler; (see above)
% University of Washington, Dept. of Physics
% Email address: charles.fieseler@gmail.com
% Website: coming soon
% Created: 10-Jul-2017
%========================================


%% Initialize with defaults

defaults = struct; %The default values
defaults.makeMovie = false; %Exports a movie (.mov)
defaults.verbose = false; %Displays extra information about the program execution
defaults.usingLegend = false; %Displays a legend on the figures
defaults.usingObjects = false; %Sometimes there are obstacles that the worm might hit; will be in a separate file
defaults.usingCOM = false; %To display Center Of Mass
defaults.usingFullWorm = true;
defaults.usingMetaDataFunction = false; %Metadata regarding stimuli might be saved
defaults.usingMetaDataExport = false; %Save the metadata figures
defaults.usingMetaDataNormalize = false; %If a negative wave-like stimulus, flip it and make the baseline 0
defaults.usingSpline = false; %This is an attempt at fitting the long-term trajectory of the animal
defaults.pauseLength = 0.0; %Duration to pause after every timestep
defaults.pauseAt = 0.0; %Pause here until user presses "enter"
defaults.startAtPause = false; %Start at the above pause value
defaults.quitAtPause = false; %Quit at the above pause value
defaults.isFolder = false; %Can pass a folder to this function; in that case it recursively calls itself for each file
defaults.downsample = 1; %Speed up visualization by throwing out data

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

[makeMovie, verbose, usingLegend, usingObjects, usingCOM,...
    usingFullWorm, usingMetaDataFunction, usingMetaDataExport,...
    usingMetaDataNormalize, usingSpline, pauseLength, pauseAt,...
    startAtPause, quitAtPause, isFolder, downsample] ...
    = v2struct(defaults); %Unpacks the struct into variables

dat = struct;
dat.settings = settings;

%% If 'isFolder' option is checked, run recursively
if isFolder
    settings.isFolder = false;
    allFileNames = struct2cell(dir(filename)); %here, 'filename' is a folder
    allFileNames = allFileNames(1,:); %Only care about the names
    allSimdataNames = allFileNames(...
        cellfun(@(x) ~isempty(strfind(x,'simdata')),allFileNames) );
    
    if filename(end)~='/'
        filename = [filename '/'];
    end
    
    fprintf('Running WormView on %d files.\n',length(allSimdataNames));
    for jName=1:length(allSimdataNames)
        close all;
        WormView([filename allSimdataNames{jName}],settings);
        pause
    end
    return;
end

%==========================================================================


%% Import worm data from *.csv file
if isempty(strfind(filename,'.'))
    filename = [filename '.csv'];
end

if isempty(strfind(filename,'.mat'))

    if exist(filename,'file')
        data = importdata(filename);
        
        if usingMetaDataExport
            try
                metaData = importdata(replace(filename,'simdata','mdata'));
                dat.metaData = metaData;
            catch
                warning('Failed to import metaData\n; turning off usingMetaDataExport option')
                usingMetaDataExport = false;
            end
        end
    
    elseif exist(['../Simdata/' filename],'file')
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

%Downsample the data
if downsample>1
    data = data(1:downsample:end,:);
end

%% Check if objects.csv exists and if so, import object data
%NOTE: if not using objects, you must delete any objects.csv file generated
%during previous runs.
if usingObjects
    Objects = importdata('../Model/objects.csv');
    Sz = size(Objects);
    Nobj = Sz(1);
end

Sz = size(data);
Nt = round(Sz(1));
Dt = data(2,1) - data(1,1);
Duration = data(end,1);
tspan = 0:Dt:Duration;

if pauseAt>0.0
    %Translate the time to pause at into an index
    pauseAtIndex = find(abs(tspan-pauseAt)<Dt,1);
else
    pauseAtIndex = -1;
end

Nbar = (Sz(2)-1)/3;
NSEG = Nbar-1;
D = 80e-6;

%Commented out values were in the original visualizer, but not used
CoM = zeros(Nt,Nbar,3);
% CoMplot = zeros(Nt,2);
Dorsal = zeros(Nbar,2);
Ventral = zeros(Nbar,2);
% Midline = zeros(Nbar,2);
% act_D = zeros(Nt,Nbar-1);
% act_V = zeros(Nt,Nbar-1);
% L_D = zeros(Nt,Nbar-1);
% L_V = zeros(Nt,Nbar-1);
% X = zeros(Nt,25);
% Y = zeros(Nt,25);

if usingMetaDataFunction
    metaData = WormMeta(filename);
    dat.metaData = metaData;
end
%---------------------------------------------
% Plot the metadata (e.g. stretch receptor weights)
%---------------------------------------------
if usingMetaDataExport
    dat.metaFigure = figure('DefaultAxesFontSize',14);
    tspan = metaData.data(:,1);
    %---------------------------------------------
    % Just show the wave of suppression if TRUE
    %--------------------------------------------- 
    yspan = 1:(size(metaData.data,2)-1);
    if usingMetaDataNormalize
        
    end
    
    metaDataVals = metaData.data(:,2:end);
    colorVals = metaDataVals;
    if pauseAt>0.0
        colorIndices = pauseAtIndex-3:pauseAtIndex+3;
        %Change the color of the time slice we're pausing at, and some
        %surrounding slices (for clarity)
        colorVals(colorIndices,:) = ...
            max(max(colorVals))*2*ones(length(colorIndices),length(yspan));
    end
    surf(yspan,tspan,metaDataVals,colorVals);
    shading interp
    colormap hot
    ylabel('Time')
    xlabel('Segment number')
    zlabel(metaData.textdata{1}(5:end)) %The first 5 characters should be 'Time '
    view(-67,34);
    drawnow;
end

%---------------------------------------------
% Make the figure for the movie visualization
%---------------------------------------------
%NOTE: Create figure then check box to fill axes!
dat.wormFigure = figure('DefaultAxesFontSize',14,'Position',[1 1 824 588]);
XYratio = 1.3333;   %This is the ratio of Xrange/Yrange required to give equal axes

R = D/2.0*abs(sin(acos(((0:Nbar-1)-NSEG./2.0)./(NSEG/2.0 + 0.2))));

for i = 1:Nt
    frame = i;
    for j = 1:Nbar
        CoM(i,j,1) = data(frame,1 + (j-1)*3 + 1);
        CoM(i,j,2) = data(frame,1 + (j-1)*3 + 2);
        CoM(i,j,3) = data(frame,1 + (j-1)*3 + 3);
    end
end

if usingCOM || usingSpline
    meanCOM = [mean(CoM(:,:,1),2), mean(CoM(:,:,2),2)];
    if usingSpline
        xSpline = linspace(min(meanCOM(:,1)),max(meanCOM(:,1)));
        ySpline = spline(meanCOM(1:2:end,1),meanCOM(1:2:end,2));
    end
end

%% Do visualization

%---------------------------------------------
% Movie settings and figure size
%---------------------------------------------
if makeMovie
    movName = [filename(1:end-4) '.avi'];
    filesepInd = strfind(movName, '/');
    videoObj=VideoWriter(['../Movies/' movName(filesepInd(end)+1:end)]);
    videoObj.FrameRate=60;
    open(videoObj);
end

simMaxX = max(max(CoM(:,:,1))) + 0.1e-3;
simMinX = min(min(CoM(:,:,1))) - 0.1e-3;
simMaxY = max(max(CoM(:,:,2)));
simMinY = min(min(CoM(:,:,2)));
simXrange = simMaxX-simMinX;
simYrange = simMaxY-simMinY;

if simXrange >= simYrange*XYratio
    plotMaxX = simMaxX;
    plotMinX = simMinX;
    Yerr = simXrange/XYratio - simYrange;
    plotMaxY = simMaxY + Yerr/2;
    plotMinY = simMinY - Yerr/2;
else
    plotMaxY = simMaxY;
    plotMinY = simMinY;
    Xerr = simYrange*XYratio - simXrange;
    plotMaxX = simMaxX + Xerr/2;
    plotMinX = simMinX - Xerr/2;
end


%---------------------------------------------
% Plot as a function of time
%---------------------------------------------
if startAtPause
    iStart = pauseAtIndex;
else
    iStart = 1;
end

for i = iStart:Nt
    
    %Plot objects if, if they are present
    if usingObjects
        angles = (0:0.05:1).*(2*pi);
        outline = zeros(length(angles),2);
        for j = 1:Nobj
            outline(:,1) = Objects(j,1) + cos(angles).*Objects(j,3);
            outline(:,2) = Objects(j,2) + sin(angles).*Objects(j,3);
            plot(outline(:,1),outline(:,2),'b','linewidth',2)
            hold on
        end
    end
    
    if usingCOM
        %note: 'i' is the current timestep
        plot(meanCOM(i,1),meanCOM(i,2),'ro','LineWidth',3);
        hold on
        if usingSpline
            plot(xSpline,ppval(ySpline,xSpline))
        end
    end
    
    if usingMetaDataFunction
        plot(metaData.steadyStateX,...
            polyval(metaData.steadyStateFit,metaData.steadyStateX),...
            'LineWidth',3);
        hold on
    end
    
    if usingFullWorm
        
        for j = 1:Nbar
            dX = R(j)*cos(CoM(i,j,3));
            dY = R(j)*sin(CoM(i,j,3));
            Dorsal(j,1) = CoM(i,j,1) + dX;
            Dorsal(j,2) = CoM(i,j,2) + dY;
            Ventral(j,1) = CoM(i,j,1) - dX;
            Ventral(j,2) = CoM(i,j,2) - dY;
        end
        
        %Plot worm
        plot(Dorsal(:,1),Dorsal(:,2),'k','linewidth',4)
        hold on
        plot(Ventral(:,1),Ventral(:,2),'k','linewidth',4)
        plot([Ventral(1,1) Dorsal(1,1)],[Ventral(1,2) Dorsal(1,2)],'k','linewidth',4)
        plot([Ventral(end,1) Dorsal(end,1)],[Ventral(end,2) Dorsal(end,2)],'k','linewidth',4)
        if usingLegend
            legend({filename},'interpreter','none');
        end
        hold off
        
    end
    
    set(gca,'xticklabel','','yticklabel','','ytick',[],'xtick',[])
    xlim([plotMinX plotMaxX])
    ylim([plotMinY plotMaxY])
    set(gcf,'paperpositionmode','auto')
    title(sprintf('Time: %f',data(i,1)));
    
    drawnow;
    if makeMovie
        frame = getframe;
        writeVideo(videoObj,frame);
    end
        
    if pauseLength>0.0
        pause(pauseLength)
    elseif pauseLength<0.0
        pause;
    end
    if i==pauseAtIndex
        pause;
        if quitAtPause
            return
        end
    end
    
end

legend({filename},'interpreter','none');

if makeMovie
    close(videoObj);
end

end





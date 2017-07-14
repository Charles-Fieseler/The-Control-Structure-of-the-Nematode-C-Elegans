classdef WormMetaDataObj < handle
    %Worm meta data object
    %   This class allows access to many different ways of visualizing the
    %   worm simulation data, including metadata for an entire folder of
    %   simulations
    %
    %
    % INPUTS (to constructor)
    %   folderName - Name of folder to take data from
    %   indVarName - Name of independent variables (1 or 2)
    %   settings   - A struct of various plotting and importing settings,
    %                which should be passed as a struct with the proper
    %                field name:
    %         verbose = true;
    %         customClassVars = {}; %User defined variables that group data runs together if present in the parameter file
    %         useTable = false; %Display a table of parameters under the plot
    %         useScatter = true; %Plot using a scatterplot
    %         useSurface = false; %Plot using surf()
    %         useContour = false; %Plot a scatter plot with a contour plot under it
    %         useErrorBars = false; %Plot using errorbar()
    %         toPlotOnConstruct = true; %Plots all classes as they are finished
    %         whichVar = 'changeInAngle'; %Which variable to plot
    %         toExit = true; %To abort on failure to find independent variable
    %         ignoreVar = {}; %Which variables to ignore when sorting data runs into different classes
    %         modulo2pi = true; %To plot the change in angle modulo 2pi or not
    %         moduloLower = -pi/2; %The lower range of the mod(2pi) turning angle
    %
    %======================================================================
    % OUTPUTS -
    %   obj - WormMetaData object with properties described in the
    %   declaration and methods
    %
    %   PUBLIC METHODS:
    %---------------------------------------------
    % function obj = plotClass(obj, classNum)
    %---------------------------------------------
            %Plot data from many simulations on one figure
            %   Uses the classes generated from the createClasses method,
            %   and only plots a single one, corresponding to "classNum"
            %
            %   The plot is by default interactive, and clicking on a data
            %   point shows the corresponding simulation (i.e. the worm's
            %   behavior)
            %
            %   Uses the plot settings in obj.plotSettings, which can be
            %   changed using obj.changePlotSetting.
            %   The dependent variable(s) are in obj.plotVar, which can be
            %   changed using obj.changePlotVar. The possible dependent
            %   variables are listed in obj.fieldsToGet.
            %   The independent variable(s) are in obj.indVarName, which
            %   can be changed using obj.changeIndVar.
    %       
    %---------------------------------------------
    % function obj = plotAllClasses(obj)
    %---------------------------------------------
            %Plots all classes, using obj.plotClass
    %
    %---------------------------------------------
    % function obj = changeIndVar(obj,indVarName)
    %---------------------------------------------
            %Changes the independent variable to be plotted
            %   This variable should be one of the parameters found in the
            %   parameter files created along with the simulation data.
            %
            %   Note: the classes are sorted according to different
            %   non-independent parameters, so within a given class, only
            %   the ignoreVar (ignored variable(s)) and the independent
            %   variables should vary.
            %
            %   The user can also create their own independent variable
            %   that takes as input some of the raw parameters in the file.
            %   See: obj.createIndVar
    %
    %---------------------------------------------
    % function obj = createIndVar(obj,userIndVarName,userIndVarFunc,varargin)
    %---------------------------------------------
            %Allows the user to create their own independent variable with
            %their own name, indVarName, and defined mathematically by the
            %user's function, indVarFunc, with a single argument that can
            %perform operations on columns
            %   e.g. userIndVarFunc = @(x)x(:,1)+x(:,2)
            %
            %   The last argument(s) are the names of the independent
            %   variable(s) that are inputs into the user-defined function.
            %   Note that these should be what was used to create the
            %   classes, i.e. they should be in or equal to obj.indVarName
    %
    %---------------------------------------------
    % function obj = changePlotVar(obj,whichVar)
    %---------------------------------------------
            % Change the dependent variable to plot
    %
    %---------------------------------------------
    % function obj = changePlotSetting(obj,thisField,thisVal)
    %---------------------------------------------
            %Changes the plot setting to be used, which are explained in
            %the "defaults" section
            %   Some settings are mutually exclusive (useScatter,
            %   useSurface, useErrorBars), and they will be turned off in
            %   favor of the new preference
    %
    %======================================================================
    %   PRIVATE METHODS
    %---------------------------------------------
    % function obj = CreateClasses(obj)
    %---------------------------------------------
            %Separates the data files found in the given folder into
            %classes based on their corresponding parameter files.
            %
            %   The default behavior is to create a new class if a unique
            %   value of any parameter was found, excluding the ignoreVar
            %   and the indVarName. 
            %
            %   The user can optionally set useCustomVar=True and pass a
            %   set of customClassVars. In that case, data files will be
            %   sorted into classes according to presence of the 
            %   customClassVars in the parameter file.
    %
    %---------------------------------------------
    % function [redP, indV] = reducePrams(obj, thisPrams, whichVars, whichIgnoreVars)
    %---------------------------------------------
            %When saving the parameter classes, get rid of the filename and the
            %independent and ignored variable(s)
            %   Note: in the parameter files, the saved filename is one of
            %   the variables that is always cut in this function
    %
    %---------------------------------------------
    % function thisClass = getMetaData(obj, simdataName, thisIndVar, thisClass)
    %---------------------------------------------
            %Get the metadata and save it in the output struct "runClasses"
            %   Also saves a reference to the data file from which it came
            %
            %   Calls the function WormMeta.m, which by default calculates
            %   and returns 5 different values:
            %       'changeInAngle'    
            %       'barVelocitiesMax'
            %       'barVelocitiesMean'
            %       'COMVelocityMax'
            %       'COMVelocityMean'
    %
    %---------------------------------------------
    % function thisClass = cleanData(obj, thisClass)
    %---------------------------------------------
            %We want to put the data in a form that meshgrid can plot it in,
            %now that we know the ranges of the values
            %   Only need this if there are 2 independent variables
            %
            %   If this process fails, some plotting options will be
            %   unavailable and it will abort with a warning.
    %
    %---------------------------------------------
    % function plotThisPoint(obj,thisField,source,eventData)
    %---------------------------------------------
            %Callback function for the data plots
            %   When a data point is clicked on, the environmental
            %   simulation corresponding to that point will run
            %
            %   If multiple fields have been plotted, only the first layer
            %   is interactive
    %
    %---------------------------------------------
    % function ptKey = pt2key(obj,ptVec)
    %---------------------------------------------
            %Takes a vector (2 or 3 dimensions) and converts it to a
            %character vector, which is what's required for the
            %map.containers object
    %
    %---------------------------------------------
    % function createErrorBars(obj)
    %---------------------------------------------
            %This file takes the raw data and sorts it by the independent
            %variables, creating a new array of means and standard
            %deviations. 
            %   This only works if there is vertically
            %   stacked data for each independent variable.
    %
    %---------------------------------------------
    % function [h]=plot3d_errorbars(obj,x, y, z, ex, ey, ez, varargin)
    %---------------------------------------------
            %FROM THIS MATHWORKS QUESTION:
            %   http://stackoverflow.com/questions/23654296/multi-dimensional-2d-better-3d-scatter-plot-with-different-errorbars-in-matlab
            %Creates a 3d errorbar plot
    %
    %
    %
    %
    %======================================================================
    % EXAMPLES
    %
    %   %%Robustness of speed to changes in internal parameters
    %   % Varying 5 environmental variables by +-1%
    % 
    %    datEnvNSR=WormMetaDataObj('../Simdata_NSR_robust/',...
    %        {'N_SR'},...
    %        struct('ignoreVar',{{'DURATION','env_vars'}},...
    %        'toPlotOnConstruct',false,'useErrorBars',true));
    % 
    %    %Should be 2 classes
    %    datEnvNSR.changePlotVar('COMVelocityMean');
    %    datEnvNSR.plotAllClasses;
    %
    %
    %
    %   %%Shows dependence of turning angle on suppression wave timescale
    %   datFine=WormMetaDataObj('../Simdata_Omwave_FineScale/',...
    %         {'SR_B_TSTART','SR_B_TEND'},...
    %         struct('ignoreVar',{{'DURATION','NMJ_B_TEND','NMJ_B_TSTART'}},...
    %         'modulo2pi',false,...
    %         'toPlotOnConstruct',false));
    %   newVarFunc = @(x) x(:,2)-x(:,1);
    %   datFine.createIndVar('Wave duration',newVarFunc,'SR_B_TSTART','SR_B_TEND');
    %
    %   datFine.plotClass(4);
    %
    %
    %======================================================================
    % Dependencies
    %   .m files, .mat files, and MATLAB products required: (updated on 10-Jul-2017)                                        
    %             MATLAB (version 9.1)
    %             Statistics and Machine Learning Toolbox (version 11.0)
    %             v2struct.m
    %             WormMeta.m
    %             WormView.m
    %             simdata_original_metadata.mat
    %             MATLAB (version 9.1)
    %             Statistics and Machine Learning Toolbox (version 11.0)
    %
    %   See also: WormMeta.m, WormView.m
    %
    %
    % Author: Charles Fieseler
    % University of Washington, Dept. of Physics
    % Email address: charles.fieseler@gmail.com
    % Website: coming soon
    % Created: 25-Jan-2017
    %======================================================================
    
    properties (SetAccess = private)
        %From caller
        foldername  %Name of folder with data and parameter files
        settings    %Struct that contains the visualization settings
        indVarName  %Name of independent variables
        ignoreVar   %Variables that are ignored in the parameter files
        customClassVars %The variables that separate out the data into classes
        
        %Initialized in constructor
        runClasses  %Output file with all the data
        fieldsToGet %Fields to be gotten from metadata function
        useCustomClasses %True if making classes in a Jacobian style, i.e. different independent variables for different files
        
        %Initialized when classes are created
        allVars     %Names of variables found in parameter files
        runClassNum %Number of classes found
        
        %Plot settings
        plotSettings
    end
    
    properties (Access = private)
        %Settings from WormMetaPlot (plotting a folder)
        toExit      %Boolean for aborting if the independent variables aren't found
        verbose     %Displays extra text relating to execution
        
        %Initialized in constructor
        pramsFiles  %Names of the parameter files
        simdataFiles%Names of the simdata files
        numWarnings %The number of files that don't have the independent variable
        numIndVar   %Number of independent variables
        indexedFiles%Map of files indexed by independent variable position
        indexedData %Struct of ordered data points, separated by dependent variable
        importFinish%Boolean of whether or not the importing is finished
        currentImport%Index of current file being imported
    end
    
    methods
        %Constructor
        function obj = WormMetaDataObj(foldername,indVarName,settings)
            %% Initialize with defaults
            defaults = struct; %The default values
            defaults.verbose = true;
            defaults.customClassVars = {}; %User defined variables that group data runs together if present in the parameter file
            defaults.useTable = false; %Display a table of parameters under the plot
            defaults.useScatter = true; %Plot using a scatterplot
            defaults.useSurface = false; %Plot using surf()
            defaults.useContour = false; %Plot a scatter plot with a contour plot under it
            defaults.useErrorBars = false; %Plot using errorbar()
            defaults.toPlotOnConstruct = true; %Plots all classes as they are finished
            defaults.whichVar = 'changeInAngle'; %Which variable to plot
            defaults.toExit = true; %To abort on failure to find independent variable
            defaults.ignoreVar = {}; %Which variables to ignore when sorting data runs into different classes
            defaults.modulo2pi = true; %To plot the change in angle modulo 2pi or not
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
            
            [obj.verbose, obj.customClassVars, ...
                obj.plotSettings.useTable,...
                obj.plotSettings.useScatter,...
                obj.plotSettings.useSurface,...
                obj.plotSettings.useContour,...
                obj.plotSettings.useErrorBars,...
                obj.plotSettings.toPlotOnConstruct,...
                obj.plotSettings.whichVar, obj.toExit, obj.ignoreVar,...
                obj.plotSettings.modulo2pi,...
                obj.plotSettings.moduloLower] ...
                = v2struct(defaults); %Unpacks the struct into variables
            
            obj.settings = settings;
            obj.foldername = foldername;
            obj.indVarName = indVarName;
            obj.importFinish = false;
            obj.useCustomClasses = ~isempty(obj.customClassVars);
            %==========================================================================
            
            %% Sort the directory into prams and simdata
            
            if isempty(foldername)
                foldername = './';
            end
            
            myDir = struct2cell(dir(foldername));
            myDir = myDir(1,:); %We only care about the names
            
            assert(~isempty(myDir),'WormMetaDataObj:notFound',...
                'Error: folder not found on path.')
            
            obj.pramsFiles = myDir(cellfun(@(x)~isempty(x),strfind(myDir,'prams')));
            obj.simdataFiles = myDir(cellfun(@(x)~isempty(x),strfind(myDir,'simdata')));
            %==========================================================================
            
            %% Preprocessing
            
            if length(obj.pramsFiles) < length(obj.simdataFiles)
                error('Some parameter files are missing; expected %d and found %d',...
                    length(obj.simdataFiles),length(obj.pramsFiles))
            elseif length(obj.pramsFiles) > length(obj.simdataFiles)
                warning('Unused parameter files found in the folder: %s',foldername)
            end
            
            % We want to put the different parameter files into groups, and only plot
            % runs together if they differ only in the independent variable
            obj.runClasses = struct;
            obj.runClassNum = 1; %We'll number the categories
            
            %---------------------------------------------
            % Set fields we'll be saving to plot
            %---------------------------------------------
            obj.fieldsToGet = {'changeInAngle','barVelocitiesMax',...
                'barVelocitiesMean','COMVelocityMax','COMVelocityMean'};
            if isempty(find(cellfun(...
                    @(q) strcmp(obj.plotSettings.whichVar,q),...
                    obj.fieldsToGet), 1))
                error('Variable to plot in obj.plotSettings.whichVar (%s) not found',...
                    obj.plotSettings.whichVar)
            end
            %Initialize the "dictionary"
            for jField=1:length(obj.fieldsToGet)
                obj.indexedFiles.(obj.fieldsToGet{jField}) = containers.Map;
            end
            %==========================================================================
            
            
            %% Import and visualize data
            
            %Create the parameter classes in the folder
            obj.CreateClasses();
            
            %Plot it all using some basic settings
            if obj.plotSettings.toPlotOnConstruct
                obj.plotAllClasses();
            end
            
            %==========================================================================
        end
        
        function obj = plotClass(obj, classNum)
            %Plot data from many simulations on one figure
            %   Uses the classes generated from the createClasses method,
            %   and only plots a single one, corresponding to "classNum"
            %
            %   The plot is by default interactive, and clicking on a data
            %   point shows the corresponding simulation (i.e. the worm's
            %   behavior)
            %
            %   Uses the plot settings in obj.plotSettings, which can be
            %   changed using obj.changePlotSetting.
            %   The dependent variable(s) are in obj.plotVar, which can be
            %   changed using obj.changePlotVar. The possible dependent
            %   variables are listed in obj.fieldsToGet.
            %   The independent variable(s) are in obj.indVarName, which
            %   can be changed using obj.changeIndVar.
            
            f=figure('DefaultAxesFontSize',14);
            
            %Get the class we're plotting here
            dataFields = fieldnames(obj.runClasses);
            thisClassName = dataFields{classNum};
            thisClass = obj.runClasses.(thisClassName);
            
            %Set which field we're going to make interactive
            thisField = obj.plotSettings.whichVar;
            if iscell(thisField)
                for jField=length(thisField):-1:1
                    thisData(jField) = thisClass.(thisField{jField});
                end
            else
                thisData = thisClass.(thisField);
            end
            
            %---------------------------------------------
            % Different plots if 1 or 2 indVar (independent variables)
            %---------------------------------------------
            if obj.numIndVar==1
                if obj.plotSettings.useTable
                    subplot(2,1,1)
                end
                
                if iscell(thisField)
                    %This means we're plotting more than one metadata field
                    hold on
                    for jField=1:length(thisField)
                        disp('Only the first legend entry of points is interactive')
                        if jField==1
                            plot(thisData(jField).raw(:,1), thisData(jField).raw(:,2),...
                                'o','LineWidth',2,'ButtonDownFcn',...
                                @(src,evt) plotThisPoint(obj,thisField,src,evt),...
                                'UserData',thisClassName)
                        else
                            plot(thisData(jField).raw(:,1), thisData(jField).raw(:,2),...
                                'o','LineWidth',2,...
                                'UserData',thisClassName)
                        end
                    end
                else
                    plot(thisData.raw(:,1), thisData.raw(:,2),...
                        'o','LineWidth',2,'ButtonDownFcn',...
                        @(src,evt) plotThisPoint(obj,thisField,src,evt),...
                        'UserData',thisClassName)
                    
                    if obj.plotSettings.useErrorBars
                        errorbar(thisData.rawWErrs(:,1),...
                            thisData.rawWErrs(:,2),...
                            thisData.rawWErrs(:,3),'LineWidth',2);
                    end
                         
                end
                
                %Various labels, depending
                %   Note that we might be sending a cell array here
                if logical(sum(strcmp(thisField,'barVelocityMax')+...
                        strcmp(thisField,'barVelocityMean')))
                    ylabel('Mean(v_{segment}) (m/s)')
                elseif logical(sum(strcmp(thisField,'changeInAngle')))
                    ylabel('\Delta\Theta_{trajectory}')
                elseif logical(sum(strcmp(thisField,'COMVelocityMax')+...
                        strcmp(thisField,'COMVelocityMean')))
                    ylabel('v_{COM} (m/s)')
                end
                legend(thisField);
                xlabel(obj.indVarName,'interpreter','none')
                
                if obj.plotSettings.useTable
                    subplot(2,1,2)
                end
                
            elseif obj.numIndVar==2
                if obj.plotSettings.useTable
                    subplot(2,1,1)
                end
                
                if iscell(thisField)
                    %This means we're plotting more than one metadata field
                    hold on
                    for jField=1:length(thisField)
                        plot3(thisData(jField).raw(:,1), thisData(jField).raw(:,2),...
                            thisData(jField).raw(:,3), 'o',...
                            'LineWidth',2,'ButtonDownFcn',...
                            @(src,evt) plotThisPoint(obj,thisField,src,evt),...
                                'UserData',thisClassName)
                    end
                else
                    if obj.plotSettings.useScatter && obj.plotSettings.useContour
                        %Draw a scatter plot with a contour plot under it
                        if ~isfield(thisData,'xx')
                            warning('In cleaning step, data was unable to be reshaped; skipping surface plotting for this class')
                            close(f);
                            return
                        end
                        dotSizes = 40*(tanh(1./abs(pi-thisData.raw(:,3)))+...
                            tanh(1./abs(pi+thisData.raw(:,3))));
                        contour(thisData.xx,thisData.yy,thisData.zz);
                        hold on;
                        scatter3(thisData.raw(:,1),thisData.raw(:,2),...
                            thisData.raw(:,3)+1,...
                            dotSizes,dotSizes,'LineWidth',2,'ButtonDownFcn',...
                            @(src,evt) plotThisPoint(obj,thisField,src,evt),...
                                'UserData',thisClassName)
                        colormap cool
                        view(0,90)
                    elseif obj.plotSettings.useScatter
                        %We want the dots to be large if it is similar to
                        %an omega turn, i.e. a change in angle of pi or -pi
                        dotSizes = 40*(tanh(1./abs(pi-thisData.raw(:,3)))+...
                            tanh(1./abs(pi+thisData.raw(:,3))));
                        scatter3(thisData.raw(:,1),thisData.raw(:,2),thisData.raw(:,3),...
                            dotSizes,dotSizes,'LineWidth',2,'ButtonDownFcn',...
                            @(src,evt) plotThisPoint(obj,thisField,src,evt),...
                                'UserData',thisClassName)
                        colormap cool
                        view(0,90)
                        
                        if obj.plotSettings.useErrorBars
                            %Uses a non-matlab plotter
                            obj.plot3d_errorbars(thisData.rawWErrs(:,1),...
                                thisData.rawWErrs(:,2),...
                                thisData.rawWErrs(:,3),...
                                zeros(size(thisData.rawWErrs(:,4))),...
                                zeros(size(thisData.rawWErrs(:,4))),...
                                thisData.rawWErrs(:,4),...
                                'ro', 'LineWidth',5)
                        end
                        
                    elseif obj.plotSettings.useSurface
                        if ~isfield(thisData,'xx')
                            warning('In cleaning step, data was unable to be reshaped; skipping surface plotting for this class')
                            close(f);
                            return
                        end
                        surf(thisData.xx,thisData.yy,thisData.zz,...
                                'UserData',thisClassName)
                        shading interp
                        colorbar
                    else
                        plot3(thisData.raw(:,1),thisData.raw(:,2),thisData.raw(:,3),...
                            'o','LineWidth',3,'ButtonDownFcn',...
                            @(src,evt) plotThisPoint(obj,thisField,src,evt),...
                                'UserData',thisClassName);
                    end
                end
                
                %Various labels, depending
                if logical(sum(strcmp(thisField,'barVelocityMax')+...
                        strcmp(thisField,'barVelocityMean')))
                    zlabel('Mean(v_{segment}) (m/s)')
                elseif logical(sum(strcmp(thisField,'changeInAngle')))
                    zlabel('\Delta\Theta_{trajectory}')
                elseif logical(sum(strcmp(thisField,'COMVelocityMax')+...
                        strcmp(thisField,'COMVelocityMean')))
                    zlabel('v_{COM} (m/s)')
                end
                legend(thisField);
                xlabel(obj.indVarName{1},'interpreter','none')
                ylabel(obj.indVarName{2},'interpreter','none')
            end
            
            %Display the other variables in a table
            if obj.plotSettings.useTable
                t=uitable(f,'Units','normalized','Position',[0.1 0.1 0.8 0.2]);
                classNames = cell(size(thisClass.Parameters.textdata,1),1);
                [classNames{:}] = thisClass.Parameters.textdata{:,1}; %Catch all the variable names
                
                t.ColumnName = classNames;
                t.Data = thisClass.Parameters.data.';
            end
            
            title(sprintf('Category %d in folder %s',...
                classNum,obj.foldername),'interpreter','none');
            
            %Pass the figure handle to the user
            thisClass.(sprintf('ClassFig_%d',classNum)) = f;
            obj.runClasses.(sprintf('Class_%d',classNum)) = thisClass;
            drawnow;
        end
        
        function obj = plotAllClasses(obj)
            %Plots all classes, using obj.plotClass
            
            dataFields = fieldnames(obj.runClasses);
            for jField1=1:length(dataFields)
                obj.plotClass(jField1);
            end
        end
        
        function obj = changeIndVar(obj,indVarName)
            %Changes the independent variable to be plotted
            %   This variable should be one of the parameters found in the
            %   parameter files created along with the simulation data.
            %
            %   Note: the classes are sorted according to different
            %   non-independent parameters, so within a given class, only
            %   the ignoreVar (ignored variable(s)) and the independent
            %   variables should vary.
            %
            %   The user can also create their own independent variable
            %   that takes as input some of the raw parameters in the file.
            %   See: obj.createIndVar
            
            if ~prod(strcmp(indVarName,obj.indVarName))
                %Change the independent variables and remake the classes
                obj.indVarName = indVarName;

                %Get rid of the old fields
                obj.runClasses = rmfield(obj.runClasses,fieldnames(obj.runClasses));
                obj.runClassNum = 1;

                obj.CreateClasses;
            else
                disp('Independent variable left unchanged')
            end
        end
        
        function obj = createIndVar(obj,userIndVarName,userIndVarFunc,varargin)
            %Allows the user to create their own independent variable with
            %their own name, indVarName, and defined mathematically by the
            %user's function, indVarFunc, with a single argument that can
            %perform operations on columns
            %   e.g. userIndVarFunc = @(x)x(:,1)+x(:,2)
            %
            %   The last argument(s) are the names of the independent
            %   variable(s) that are inputs into the user-defined function.
            %   Note that these should be what was used to create the
            %   classes, i.e. they should be in or equal to obj.indVarName
            
            numUserArgs = length(varargin);
            if iscell(userIndVarName)
                numUserVars = length(userIndVarName);
            elseif ischar(userIndVarName)
                numUserVars = 1;
            else
                error('Input cell array or string for the desired variable name(s)')
            end
            obj.numIndVar = numUserVars;
            
            %Check that everything is a good name
            for jV = 1:numUserArgs
                assert(ismember(varargin{jV},obj.allVars),...
                    'Not a valid argument; name not in variable list.\n');
            end
            
            obj.allVars = [obj.allVars; {userIndVarName}];
            
            %If we already have the correct independent variables, we don't
            %have to remake the classes
            if isequal(varargin,obj.indVarName) ||...
                    isequal(varargin{1},obj.indVarName)
                
                %Overwrite all the classes old independent variables with
                %the new one (a column)
                for jF = 1:length(obj.fieldsToGet)
                    thisField = obj.fieldsToGet{jF};
                    
                    for jC = 1:obj.runClassNum-1
                        thisClass = sprintf('Class_%d',jC);
                        oldDat = obj.runClasses.(thisClass).(thisField).raw;
                        newIndVar = userIndVarFunc(oldDat(:,1:numUserArgs));
                        newDat = [newIndVar oldDat(:,end)];
                        
                        %Need to change the dictionary separately
                        for jV = 1:size(newIndVar,1)
                            newkey = obj.pt2key(newDat(jV,:));
                            oldkey = obj.pt2key(oldDat(jV,:));
                            oldFile = ...
                                obj.runClasses.(thisClass).(thisField).indexedFiles(oldkey);
                            remove(obj.runClasses.(thisClass).(thisField).indexedFiles,oldkey);
                            obj.runClasses.(thisClass).(thisField).indexedFiles(newkey) =...
                                oldFile;
                        end
                        %Note that we keep the dependent variable the same,
                        %i.e. the last column of the raw data
                        obj.runClasses.(thisClass).(thisField).raw = newDat;
                        
%                         obj.runClasses.(thisClass).(thisField).raw = ...
%                             [newIndVar ...
%                             obj.runClasses.(thisClass).(thisField).raw(:,end)];
                        
%                         if numUserVars == 2
%                             obj.indexedFiles.(thisField) = ...
%                                 [obj.indexedFiles.(thisField); ...
%                                 mat2cell(newDat,ones(sz))];
%                         elseif numUserVars == 1
%                             obj.indexedFiles.(thisField) = ...
%                                 [obj.indexedFiles.(thisField); ...
%                                 mat2cell([newDat zeros(sz)],ones(sz))];
%                         else
%                             error('More than 2 independent variables not supported')
%                         end
                    end                    
                end
                
            else
                error('Declare the variables to be used in the user-defined function as ''indVar'' first')
%                 for jC = 1:obj.runClassNum-1
%                     thisClass = obj.runClasses(jC);
%                     dat = thisClass.(thisField).raw;
%                 end
%                 
%                 %We only want to remake the new variable, not all of them
%                 oldFields = obj.fieldsToGet;
%                 obj.fieldsToGet = userIndVarName;
% 
%                 obj.fieldsToGet = [oldFields userIndVarName];
            end
            

            obj.numIndVar = 1;
            obj.indVarName = userIndVarName;
            
        end
        
        function obj = changePlotVar(obj,whichVar)
            % Change the dependent variable to plot
            
            assert(logical(...
                prod(ismember(whichVar,obj.fieldsToGet))),...
                'Not a valid dependent variable. Options are in obj.fieldsToGet')
            
            obj.plotSettings.whichVar = whichVar;
            if obj.plotSettings.toPlotOnConstruct
                obj.plotAllClasses();
            end
        end
        
        function obj = changePlotSetting(obj,thisField,thisVal)
            %Changes the plot setting to be used, which are explained in
            %the "defaults" section
            %   Some settings are mutually exclusive (useScatter,
            %   useSurface, useErrorBars), and they will be turned off in
            %   favor of the new preference
            
            assert(logical(...
                sum(ismember(thisField,fieldnames(obj.plotSettings)))),...
                'Not a valid plot setting. Options are the fields of obj.plotSettings')
            obj.plotSettings.(thisField) = thisVal;
            %Don't want to have competing plot settings
            if strcmp(thisField,'useScatter') && thisVal==true
                obj.plotSettings.useSurface = false;
            elseif strcmp(thisField,'useSurface') && thisVal==true
                obj.plotSettings.useScatter = false;
            end
            if strcmp(thisField,'useErrorBars') && thisVal==true
                obj.createErrorBars();
            end
            
            if obj.plotSettings.toPlotOnConstruct
                obj.plotAllClasses();
            end
        end
        
    end
    
    methods (Access = private)

        %---------------------------------------------
        % Separates the files into parameter classes
        %---------------------------------------------
        function obj = CreateClasses(obj)
            %Separates the data files found in the given folder into
            %classes based on their corresponding parameter files.
            %
            %   The default behavior is to create a new class if a unique
            %   value of any parameter was found, excluding the ignoreVar
            %   and the indVarName. 
            %
            %   The user can optionally set useCustomVar=True and pass a
            %   set of customClassVars. In that case, data files will be
            %   sorted into classes according to presence of the 
            %   customClassVars in the parameter file.
            
            %---------------------------------------------
            % Make sure either one independent variable or two
            %---------------------------------------------
            if isa(obj.indVarName,'cell')
                switch length(obj.indVarName)
                    case 1
                        obj.indVarName = obj.indVarName{1}; %Just want the string
                        obj.numIndVar = 1;
                    case 2
                        %We'll check later to see if they're valid variables
                        obj.indVarName = obj.indVarName;
                        obj.numIndVar = 2;
                    otherwise
                        if ~obj.useCustomClasses
                            error('More than two independent variables not supported... Did you mean to set "useCustomClasses" to true?')
                        else
                            %Even though there are multiple independent
                            %variables, we'll be using them one at a time
                            obj.numIndVar = 1;
                        end
                end
            elseif ~isa(obj.indVarName,'char')
                error('Datatype %s not supported for indVar; input a cell array or string',...
                    class(obj.indVarName))
            else
                obj.numIndVar=1;
            end
            
            %---------------------------------------------
            % Separate the files and import them
            %---------------------------------------------
            numFiles = length(obj.simdataFiles);
            if obj.verbose
                approxTime = round(numFiles/100)/10;
                fprintf('Found %d files; may take ~%.1f minute(s)\n',...
                    numFiles, approxTime);
            end
            for iFile = 1:numFiles
                obj.currentImport = iFile;
                
                %Note that both lists are sorted, but we should check
                thisSimdataName = obj.simdataFiles{iFile};
                thisPramsName = obj.pramsFiles{iFile};
                
                assert(strcmp(thisSimdataName(8:end),thisPramsName(6:end)),...
                    'Simdata (%s) and Parameter file (%s) names do not match.',...
                    thisSimdataName,thisPramsName);
                
                %Import the parameter file
                thisFullName = [obj.foldername thisPramsName];
                thisPrams = importdata(thisFullName);
                thisPramsTable = readtable(thisFullName,...
                    'ReadVariableNames',false);
                if ~isequal( size(thisPramsTable),size(thisPrams) )
                    %This means we have non-numeric data in the file
                    thisPrams.textdata = cell(size(thisPramsTable));
                    thisPrams.data = zeros(size(thisPramsTable,1)-1,1);
                    %First read the filename
                    thisPrams.textdata{1,1} = ...
                        thisPramsTable{1,1}{1};
                    
                    for jField =2:size(thisPramsTable,1)
                        %Read out the variable names as strings, skipping
                        %the filename
                        thisPrams.textdata{jField,1} = ...
                            thisPramsTable{jField,1}{1};
                        %Read out the variable names as numbers if
                        %possible, but we're only here if a boolean messed
                        %us up
                        
                        thisPrams.data(jField-1) = ...
                            str2double(cell2mat(thisPramsTable{jField,2}));
                        if isequaln(thisPrams.data(jField-1),NaN)
                            if strcmp(thisPramsTable{jField,2}{1},'true')
                                thisPrams.data(jField-1) = 1;
                            elseif strcmp(thisPramsTable{jField,2}{1},'false')
                                thisPrams.data(jField-1) = 0;
                            else
                                error('Unparsable data in parameter file %s\n',...
                                    thisFullName)
                            end
                        end
                    end
                end
                
                %Check the first column to see if we've actually used the user's
                %'indVarName'
                thisVars = thisPrams.textdata(:,1);
                
                %Save the full parameter list of the first file (all the
                %files should have the same parameters, unless the user 
                %created custom classes)
                if isempty(obj.allVars)
                    obj.allVars = thisVars;
                end
                
                %Get the index of the independent variable(s) in this data
                %file
                if obj.numIndVar==1
                    whichVars = contains(thisVars,obj.indVarName);
                elseif obj.numIndVar==2
                    whichVars = cellfun(...
                        @(x) logical(strcmp(x,obj.indVarName{1})+...
                        strcmp(x,obj.indVarName{2})),...
                        thisVars);
                end  
                if obj.useCustomClasses
                    assert(~isempty(obj.customClassVars),...
                        'Must input class variables in setting "customClassVars"')
                    assert(isempty(find(contains(obj.customClassVars,obj.indVarName), 1)),...
                        '"customClassVars" shouldn''t overlap with the independent variable')
                    whichClassVars = contains(thisVars,obj.customClassVars);
                    foundVars = length(find(whichClassVars));
                    if foundVars>1
                        firstVar = find(whichClassVars,1);
                        whichClassVars = zeros(size(whichClassVars));
                        whichClassVars(firstVar) = 1;
                        warning('Found multiple independent variables, using just the first one')
                    elseif foundVars == 0
                        error('No custom class variables found in data file:\n %s',...
                            thisPramsName)
                    end
                    %To identify the class
                    thisClassVar = thisVars(whichClassVars);
                    thisClassVar = thisClassVar{1};
                else
                    thisClassVar = '';
                end
                
                whichIgnoreVars = zeros(size(whichVars));
                if ~isempty(obj.ignoreVar)
                    %Check for special strings which mean a hard-coded list
                    %of ignored variables
                    ignoreInd = contains(obj.ignoreVar,'env_vars');
                    if ~isempty(find(ignoreInd,1))
                        env_vars = ...
                            {'L_SEG','K_PE','D_PE','AE_PE_RATIO',...
                            'K_DE','D_DE','T_MUSCLE'};
                        obj.ignoreVar = ...
                            [env_vars obj.ignoreVar(~ignoreInd)];
                    end
                    
                    %Function to compare the ignored variables to the ones
                    %we're using
                    f = @(var)cellfun(@(x)strcmp(x,var),thisVars);
                    %First get rid of the independent variables, then the ones we want
                    %to ignore
                    for jVar = 1:length(obj.ignoreVar)
                        whichIgnoreVars = whichIgnoreVars+f(obj.ignoreVar{jVar});
                    end
                end
                
                %We might only have found one of the two requested variables
                if length(find(whichVars))~=obj.numIndVar
                    if isempty(obj.numWarnings)
                        obj.numWarnings = 0;
                    else
                        obj.numWarnings = obj.numWarnings+1;
                    end
                    if obj.numWarnings<1
                        warning(['Did''t find the independent variable in the parameter file for %s\n'...
                            'Skipping plotting step and suppressing further warnings.\n'],...
                            thisSimdataName)
                    elseif obj.toExit
                        disp('Exiting and returning variables found in out.allVars')
                        break;
                    end
                else
                    
                    %Get the reduced data file, which doesn't include the independent
                    %variable, the filename, or any ignored variables
                    %   As well as the independent variable(s) by them/itself
                    [thisRedPrams, thisIndVar] = ...
                        reducePrams(obj, thisPrams,whichVars,whichIgnoreVars);
                    
                    %Next, search to see which data run class we're in
                    runClassNames = fieldnames(obj.runClasses);
                    %Get rid of the saved figure fields
                    runClassNames(cellfun(@(x) isempty(strfind(x,'Fig')),runClassNames));
                    
                    %Then we have the independent variable, and we want to see if the
                    %other parameters are the same as a previously read data run
                    for jClass = 1:length(runClassNames)
                        thisName = runClassNames{jClass};
                        checkClassParams = obj.runClasses.(thisName).Parameters;
                        %Only use the next if useCustomClasses=true
                        checkClassVar = obj.runClasses.(thisName).ClassVar;
                        if isequal(thisRedPrams,checkClassParams) ||...
                                (obj.useCustomClasses && strcmp(thisClassVar,checkClassVar))
                            %Get the number on the end of the class name; that's the
                            %figure we want to plot on
                            thisClassNum = str2double(thisName(strfind(thisName,'_')+1:end));
                            thisClassName=sprintf('Class_%d',thisClassNum);
                            %Save the parameters
                            obj.runClasses.(thisClassName).Parameters = thisRedPrams;
                            %Get the data, which will be returned in various data fields
                            obj.runClasses.(thisClassName) = ...
                                obj.getMetaData(thisSimdataName, ...
                                thisIndVar, obj.runClasses.(thisClassName));
                            break
                        end
                        
                        if jClass==length(runClassNames)
                            obj.runClassNum = obj.runClassNum + 1;
                            %This means we need to create a new class
                            thisClassName=sprintf('Class_%d',obj.runClassNum);
                            %Save the parameters
                            obj.runClasses.(thisClassName).Parameters = thisRedPrams;
                            obj.runClasses.(thisClassName).ClassVar = thisClassVar;
                            %Get the data, which will be returned in various data fields
                            obj.runClasses.(thisClassName) = ...
                                obj.getMetaData(thisSimdataName, ...
                                thisIndVar, obj.runClasses.(thisClassName));
                        end
                    end
                    
                    if isempty(runClassNames)
                        %This means we need to create the first class of its type
                        thisClassName=sprintf('Class_%d',obj.runClassNum);
                        %Save the parameters
                        obj.runClasses.(thisClassName).Parameters = thisRedPrams;
                        obj.runClasses.(thisClassName).ClassVar = thisClassVar;
                        %Get the data, which will be returned in various data fields
                        obj.runClasses.(thisClassName) = ...
                            obj.getMetaData(thisSimdataName, ...
                            thisIndVar, obj.runClasses.(thisClassName));
                    end
                end
            end
            
            obj.importFinish = true;
            
            %---------------------------------------------
            % Clean the data
            %---------------------------------------------
            dataFields = fieldnames(obj.runClasses);
            for jField1=1:length(dataFields)
                obj.runClasses.(dataFields{jField1}) = ...
                    obj.cleanData(obj.runClasses.(dataFields{jField1}));
            end
            
            if obj.plotSettings.useErrorBars
                obj.createErrorBars;
            end
        end
        
        function [redP, indV] = reducePrams(obj, thisPrams, whichVars, whichIgnoreVars)
            %When saving the parameter classes, get rid of the filename and the
            %independent and ignored variable(s)
            %   Note: in the parameter files, the saved filename is one of
            %   the variables that is always cut in this function
            
            if size(thisPrams.data,1)>2
                redP.textdata = thisPrams.textdata(~(whichVars+whichIgnoreVars),:);
                redP.textdata = redP.textdata(2:end,:);
                redP.data = thisPrams.data(~(whichVars(2:end)+whichIgnoreVars(2:end)),:);
            else
                %This means there was only the filename and the independent
                %variable
                redP.textdata = '';
                redP.data = '';
            end
            
            %Get just the independent variable(s)
            if obj.numIndVar==1
                indV = thisPrams.data(whichVars(2:end),:);
            else %The variables might be out of order!
                indV = zeros(obj.numIndVar,1);
                for jVar = 1:obj.numIndVar
                    indV(jVar) = ...
                        thisPrams.data(ismember(obj.allVars(2:end),obj.indVarName{jVar}),:);
                end
            end
        end
        
        function thisClass = getMetaData(obj, simdataName, thisIndVar, thisClass)
            %Get the metadata and save it in the output struct "runClasses"
            %   Also saves a reference to the data file from which it came
            %
            %   Calls the function WormMeta.m, which by default calculates
            %   and returns 5 different values:
            %       'changeInAngle'    
            %       'barVelocitiesMax'
            %       'barVelocitiesMean'
            %       'COMVelocityMax'
            %       'COMVelocityMean'
            
            %Get the metadata for this data run
            fullFileName = [obj.foldername simdataName];
            dat = WormMeta(fullFileName,...
                struct('toPlot',false,...
                'modulo2pi',obj.plotSettings.modulo2pi,...
                'moduloLower',obj.plotSettings.moduloLower));
            
            %Get 5 types of data (might not plot all)
            %Also, create dictionary of correspondence between data files
            %and data points
            for jField=1:length(obj.fieldsToGet)
                thisField = obj.fieldsToGet{jField};
                if ~isfield(thisClass,thisField)
                    thisClass.(thisField).raw = [];
                    thisClass.(thisField).indexedData = [];
                    thisClass.(thisField).indexedFiles = containers.Map;
                end
                %Append the data from this file
                %   Each column is an independent variable
                %Also, build the struct of filenames that correspond to
                %each data point, either x and y or x, y, and z
                if obj.numIndVar==1
                    thisCoord = [thisIndVar dat.(thisField)];
                    thisClass.(thisField).raw = ...
                        [thisClass.(thisField).raw; thisCoord];
                elseif obj.numIndVar==2
                    thisCoord = [thisIndVar(1) thisIndVar(2) dat.(thisField)];
                    thisClass.(thisField).raw = ...
                        [thisClass.(thisField).raw; thisCoord];
                end
                
                %Dictionary file for plot interactivity
                thisClass.(thisField).indexedData =...
                    [thisClass.(thisField).indexedData; thisCoord];
                thisClass.(thisField).indexedFiles(obj.pt2key(thisCoord)) = ...
                    fullFileName;
            end
            
            
        end
        
        function thisClass = cleanData(obj, thisClass)
            %We want to put the data in a form that meshgrid can plot it in,
            %now that we know the ranges of the values
            %   Only need this if there are 2 independent variables
            %
            %   If this process fails, some plotting options will be
            %   unavailable and it will abort with a warning.
            
            alreadyWarned = false;
            
            for jField=1:length(obj.fieldsToGet)
                thisField = obj.fieldsToGet{jField};
                dat = thisClass.(thisField).raw;
                %Want to redo this field as three arrays, not one
                
                %If only one independent variable, just save it without
                %reshaping
                if obj.numIndVar==1
                    continue
                end
                
                %Get the number of different x and y parameters
                sz = ones(1,2);
                %Can't just use "unique" because there might be
                %repetitions... so we need to find when the parameters
                %start to repeat, and that is the size we want
                %   the diff() function returns 1 when we switch from 0 to
                %   1, i.e. return to the first value of dat
                sz1 = find(diff( dat(:,1)==dat(1,1) )==1,1);
                sz2 = find(diff( dat(:,2)==dat(1,2) )==1,1);
                if ~isempty(sz1)
                    sz(1) = sz1;
                    sz(2) = length(dat(:,1))/sz1;
                elseif ~isempty(sz2)
                    %This means we only had a consecutive run of the first
                    %value in the column
                    sz(2) = sz2;
                    sz(1) = length(dat(:,1))/sz2;
                else
                    %This means that neither of the parameters were in a
                    %repeating block, i.e. the data wasn't taken in a
                    %square... Next check if the data was taken as a
                    %trapezoid (or a square but not in these variables)
                    
                    %Find consecutive values that are the same
                    sz1 = find(diff( dat(:,1) ),1);
                    sz2 = find(diff( dat(:,2) ),1);
                    if ~isempty(sz1)
                        sz(1) = sz1;
                        sz(2) = length(dat(:,1))/sz1;
                    elseif ~isempty(sz2)
                        %This means we only had a consecutive run of the first
                        %value in the column
                        sz(2) = sz2;
                        sz(1) = length(dat(:,1))/sz2;
                    else
                        %This means we don't have rectangular data
                        if ~alreadyWarned
                            warning('Data not collected in a format that surf() can use for field %s',thisField)
                            alreadyWarned=true;
                        end
                        %obj.plotSettings.useSurface = 0;
                    end
                end
                
                try
                    thisClass.(thisField).xx = reshape(dat(:,1),sz);
                    thisClass.(thisField).yy = reshape(dat(:,2),sz);
                    thisClass.(thisField).zz = reshape(dat(:,3),sz);
                catch
                    if ~alreadyWarned
                        warning(['Failed to reshape data.'...
                            'Check the independent and ignored variables'])
                        alreadyWarned = true;
                    end
                end
            end
        end
        
        function plotThisPoint(obj,thisField,source,eventData)
            %Callback function for the data plots
            %   When a data point is clicked on, the environmental
            %   simulation corresponding to that point will run
            %
            %   If multiple fields have been plotted, only the first layer
            %   is interactive
            
            if iscell(thisField)
                %Only the first cell layer is interactive
                thisField = thisField{1};
            end
            thisClass = source.UserData;            
            thisFile = obj.runClasses.(thisClass).(thisField).indexedFiles(...
                obj.pt2key(eventData.IntersectionPoint));
            
            options = struct('usingCOM',false,'downsample',1);
            WormView(thisFile,options);
        end
        
        function ptKey = pt2key(obj,ptVec)
            %Takes a vector (2 or 3 dimensions) and converts it to a
            %character vector, which is what's required for the
            %map.containers object
            
            assert(length(ptVec)==2 || length(ptVec)==3,...
                'Can only convert 2d or 3d vectors to keys, which should be length=numIndVar')
            
            if length(ptVec)==2
                ptVec = [ptVec 0 ];
            end
            
            ptKey = '';
            
            for jL = 1:length(ptVec)
                ptKey = [ptKey sprintf('%.3e',ptVec(jL))];
            end
        end
        
        function createErrorBars(obj)
            %This file takes the raw data and sorts it by the independent
            %variables, creating a new array of means and standard
            %deviations. 
            %   This only works if there is vertically
            %   stacked data for each independent variable.
            
            cNames = fieldnames(obj.runClasses);
            
            for jC = 1:length(cNames)
                thisClass = obj.runClasses.(cNames{jC});
                fNames = obj.fieldsToGet;
                for jF = 1:length(fNames)
                    dat = thisClass.(fNames{jF}).raw;
                    %Get all possible singles or pairs of ind variables
                    if obj.numIndVar>1
                        allInds = unique(dat(:,1:end-1),'rows');
                        sz = [size(allInds,1) 1];
                        allMeans = zeros(sz); 
                        allErrs = zeros(sz);
                        
                        for jI = 1:sz(1)
                            %Get all the data points that are stacked on
                            %this independent variable coordinate pair
                            thisDat = ...
                                dat( logical( prod( ...
                                dat(:,1:end-1)==allInds(jI,:),2)),end);
                            allMeans(jI) = mean(thisDat);
                            if isnan(allMeans(jI))
                                warning('Removing nan''s in the data')
                                allMeans(jI) = nanmean(thisDat);
                                allErrs(jI) = nanstd(thisDat);
                            else
                                allErrs(jI) = std(thisDat);
                            end
                        end
                        
                    else
                        allInds = unique(dat(:,1));
                        allMeans = zeros(size(allInds)); 
                        allErrs = zeros(size(allInds));
                        
                        for jI = 1:length(allInds)
                            %Get all the data points that are stacked on
                            %this independent variable coordinate
                            thisDat = dat(dat(:,1)==allInds(jI),end);
                            allMeans(jI) = mean(thisDat);
                            allErrs(jI) = std(thisDat);
                        end
                    end
                    
                    obj.runClasses.(cNames{jC}).(fNames{jF}).rawWErrs = ...
                        [allInds allMeans allErrs];
                    
                end
            end
            
        end
        
        function [h]=plot3d_errorbars(obj,x, y, z, ex, ey, ez, varargin)
            %FROM THIS MATHWORKS QUESTION:
            %   http://stackoverflow.com/questions/23654296/multi-dimensional-2d-better-3d-scatter-plot-with-different-errorbars-in-matlab
            %Creates a 3d errorbar plot
            
            % create the standard 3d scatterplot
            hold off;
            h=plot3(x, y, z, varargin{:});
            hold on
            
            % now draw the vertical errorbar for each point
            for i=1:length(x)
                xV = [x(i); x(i)];
                yV = [y(i); y(i)];
                zV = [z(i); z(i)];
                
                xMin = x(i) + ex(i);
                xMax = x(i) - ex(i);
                yMin = y(i) + ey(i);
                yMax = y(i) - ey(i);
                zMin = z(i) + ez(i);
                zMax = z(i) - ez(i);
                
                xB = [xMin, xMax];
                yB = [yMin, yMax];
                zB = [zMin, zMax];
                
                % draw error bars
                h=plot3(xV, yV, zB, '-k');
                set(h, 'LineWidth', 2);
                h=plot3(xB, yV, zV, '-k');
                set(h, 'LineWidth', 2);
                h=plot3(xV, yB, zV, '-k');
                set(h, 'LineWidth', 2);
            end
        end
        
    end
end


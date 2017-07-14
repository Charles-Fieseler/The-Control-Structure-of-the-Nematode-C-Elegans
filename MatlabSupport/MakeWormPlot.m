%% MakeWormPlot
% Run this to produce figure plots
%   This function pauses after each plot to allow the user to determine
%   figure angles, etc. if the "toCheck" option is true
%   If "toSave" is true, then the program will save these figures as .tikz
%   (a latex format) and .png. Most of the figures are saved with and
%   without Matlab text.
%
%   Note about dependencies: most of these are from the "matlab2tikz"
%   function, and they are packaged together.
%
% OUTPUTS -
%   OUTPUT1 - Omega turn plots
%
% Dependencies
%   Other m-files required: (updated on 10-Jul-2017)                            
%         MATLAB (version 9.1)
%         Statistics and Machine Learning Toolbox (version 11.0)
%         v2struct.m
%         m2tInputParser.m
%         matlab2tikz.m
%         errorUnknownEnvironment.m
%         getEnvironment.m
%         guitypes.m
%         isAxis3D.m
%         isVersionBelow.m
%         m2tUpdater.m
%         m2tstrjoin.m
%         versionArray.m
%         versionString.m
%         WormMeta.m
%         WormMetaDataObj.m
%         WormView.m
%         simdata_original_metadata.mat
%   Subfunctions:
%
%   See also: WormView.m
%
%
% Author: Charles Fieseler
% University of Washington, Dept. of Physics
% Email address: charles.fieseler@gmail.com
% Website: coming soon
% Created: 13-Feb-2017
%========================================


%% Check a file by hand, if not already done
%When does the interesting behavior start and end?

toCheck = false;
toSave = false;

if toCheck
    filename = ...
        '../Model_Parallel/simdata_om_tst_2.8_tend7.3_NMJw0.2_tyr.csv';
    
    options = struct(...
        'pauseLength',-1,...%This means we pause at every frame
        'downsample',5);
    
    dat = WormView(filename,options);
    return
end

%==========================================================================



%% Make and save different slices of an omega turn
%Using ~5 frames from a single data run
close all;

options = struct(...
    'usingMetaDataExport',true,...
    'pauseAt',2.5,...
    'quitAtPause',true,...
    'startAtPause',true);

%Note that these times are decided by hand from the above section
pauseTimes = 2.5:1.0:8.5;
%This data run is also decided by hand (i.e. it's pretty)
filename = ...
    '../Simdata_Tyramine_Real/simdata_om_tst_2.8_tend7.3_NMJw0.2_tyr.csv';

for iPause = 1:length(pauseTimes)
    options.pauseAt = pauseTimes(iPause);
    dat=WormView(filename,options);
    if toSave
        saveFilename = ...
            sprintf('tst_2.8_tend7.3_pause%.1f.jpg',pauseTimes(iPause));
         saveas(dat.metaFigure,['./PaperPlots/metafig_' saveFilename '.jpg']);
         saveas(dat.wormFigure,['./PaperPlots/wormfig_' saveFilename '.jpg']);
        %Also save a text-free version
        title(dat.wormFigure.CurrentAxes,'');
        saveas(dat.wormFigure,['./PaperPlots/wormfig_' saveFilename '_noText.jpg']);
    end
end

%---------------------------------------------
% Make a wave-moving plot
%---------------------------------------------
tspan = dat.metaData.data(:,1);
Dt = tspan(2)-tspan(1);
yspan = 1:(size(dat.metaData.data,2)-1);

%Get the time indices where we'll be pausing
for jPause=1:length(pauseTimes)
    whichPause(jPause) = ...
        find(abs(tspan-pauseTimes(jPause))<Dt,1);
end

%---------------------------------------------


%% Get the wave, but get rid of the baseline slope and normalize
%---------------------------------------------
metaDataVals = abs(dat.metaData.data(:,2:end)-dat.metaData.data(1,2:end));
metaDataVals = metaDataVals./max(metaDataVals);
%metaDataVals = metaDataVals(whichPause);
%metaDataVals = metaData.data(:,2:end);
colorVals = metaDataVals;
%Plot the waveform at the different pause times, but each with
%different colors
f = figure('DefaultAxesFontSize',14);
hold on
for jPause=1:length(pauseTimes)
    thisPause = whichPause(jPause);
    h{jPause} = waterfall(yspan,tspan(thisPause),metaDataVals(thisPause,:));
    set(h{jPause},'LineWidth',3)
    set(h{jPause},'FaceAlpha',0.5)
    %set(h{jPause},'EdgeColor',[thisPause/max(whichPause);0;0])
end
ylabel('Time')
xlabel('Segment number')
zlabel(dat.metaData.textdata{1}(5:end)) %The first 5 characters should be 'Time '
view(17,57);
drawnow;
if toSave
    saveas(f,['./PaperPlots/metafig_slices' '.jpg']);
end

%==========================================================================


%% Make an interpolated plot
%   Pretty much the same as above, but highlight the paused slices instead
%   of isolating them
f = figure('DefaultAxesFontSize',14);
surfl(yspan,tspan,metaDataVals);
xlim([1 12]);
ylim([2 8]);
alpha(0.2);
hold on
for jPause=1:length(pauseTimes)
    thisPause = whichPause(jPause);
    h{jPause} = waterfall(yspan,tspan(thisPause),metaDataVals(thisPause,:));
    set(h{jPause},'LineWidth',3)
    set(h{jPause},'FaceAlpha',0.5)
    %set(h{jPause},'EdgeColor',[thisPause/max(whichPause);0;0])
end
shading interp
colormap hot
ylabel('Time')
xlabel('Segment number')
%zlabel(dat.metaData.textdata{1}(5:end)) %The first 5 characters should be 'Time '
zlabel('Fraction of suppression')
view(0,80);
title('Wave of Suppression on the Stretch Receptors')
drawnow;

pause

if toSave
    saveas(f,['./PaperPlots/metafig_slices_interp' '.jpg']);
end

%Clear all the text and save a raw version
title('')
ylabel('')
xlabel('')
zlabel('')
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
if toSave
    saveas(f,['./PaperPlots/metafig_slices_interp_noText' '.jpg']);
end


%==========================================================================


%% Make a floating 3d plot
%Use pcolor and a floating surfl

f2 = figure('DefaultAxesFontSize',14);
%xlim manual
surfl(yspan,tspan,metaDataVals);
xlim([1 12]);
ylim([2 8]);

surfl(yspan,tspan,metaDataVals+1); 
hold on
pcolor(yspan,tspan,metaDataVals);

shading interp
colormap(gray)
ylabel('Time')
xlabel('Segment number')
zlabel('Fraction of suppression')
view(-25,26);
title('Wave of Suppression on the Stretch Receptors')
drawnow;

pause

if toSave
    saveas(f2,['./PaperPlots/metafig_pcolor_interp' '.jpg']);
end

%Clear all the text and save a raw version
title('')
ylabel('')
xlabel('')
zlabel('')
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
if toSave
    saveas(f2,['./PaperPlots/metafig_pcolor_interp_noText' '.jpg']);
end


%==========================================================================


%% Save a large data run
%Using ~100 files
%   For now, just the change in angle (not speed or anything)

datFine=WormMetaDataObj('../Simdata_Omwave_FineScale/',...
    {'SR_B_TSTART','SR_B_TEND'},...
    struct('ignoreVar',{{'DURATION','NMJ_B_TEND','NMJ_B_TSTART'}},...
    'toPlotOnConstruct',false));

toPlot = 5;
datFine.plotClass(toPlot);
view(5,55)
saveas(datFine.runClasses.(sprintf('Class_%d',toPlot)).ClassFig_5,...
    ['./PaperPlots/Simdata_Omwave_FineScale_Class' num2str(toPlot) '.jpg']);
%==========================================================================


%% Number of stretch receptors
if ~exist('datNSR','var')
    datNSR=WormMetaDataObj('../Simdata_N_SR_2nd/',...
        {'N_SR','MEDIUM'},...
        struct('ignoreVar',{{'DURATION'}},...
        'toPlotOnConstruct',false));
end

%Average velocity
datNSR.changePlotVar('COMVelocityMean');
toPlot = 1;
datNSR.plotClass(toPlot);
f = datNSR.runClasses.(sprintf('Class_%d',toPlot)).ClassFig_1;

xlabel('Number of stretch receptors')
ylabel('Viscosity parameter')

view(2,22)

saveas(f,['./PaperPlots/Simdata_Omwave_NSR_Class' num2str(toPlot) '.jpg']);

%Also save a text-free version
xlabel('')
ylabel('')
zlabel('')
legend('')
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
title('')
saveas(datNSR.runClasses.(sprintf('Class_%d',toPlot)).ClassFig_1,...
    ['./PaperPlots/Simdata_Omwave_NSR_Class' num2str(toPlot) '_noText.jpg']);


%Now do a 2d version
datNSR.changeIndVar('N_SR');
toPlot = 8;
datNSR.plotClass(toPlot);
f = datNSR.runClasses.(sprintf('Class_%d',toPlot)).ClassFig_8;
xlabel('Stretch receptor length')
title('Optimized stretch receptor length')
pause

%Also export to .tex
matlab2tikz('figurehandle',f,'filename','./PaperPlots/N_SRplot.tex')
%==========================================================================



%% Difference in forward and backwards weightings
%---------------------------------------------
% NMJ weights
%---------------------------------------------
NSEG = 48;
xspan = 1:NSEG;

NMJ_weight_B = 0.7*(1-xspan*0.6/NSEG);
NMJ_weight_B(1) = NMJ_weight_B(1)/1.5;
NMJ_weight_A = 0.7*(1-(NSEG-xspan-1)*0.6/NSEG);
NMJ_weight_A(end) = NMJ_weight_A(end)/1.5;

f = figure('DefaultAxesFontSize',14);
subplot(2,1,1)
plot(xspan./NSEG,NMJ_weight_B,'LineWidth',2)
hold on
plot(xspan./NSEG,NMJ_weight_A,'LineWidth',2)

title('Different NMJ gradients required for forward vs. backward motion')
legend('B-class','A-class')
%xlabel('Body length')
ylabel('NMJ weighting')

%---------------------------------------------
% Stretch receptor weightings
%---------------------------------------------

N_units = 12;
N_seg_per_unit = NSEG/N_units;
xspan = 1:N_units;

SR_weight_A = 0.65*(0.4 + 0.08*(N_units-xspan-1))*(N_units/12.0)*(2.0/N_seg_per_unit);
SR_weight_B = (0.4 + 0.08*xspan)*(N_units/12.0)*(2.0/N_seg_per_unit);

%f2 = figure('DefaultAxesFontSize',14)
subplot(2,1,2)
plot(xspan./N_units,SR_weight_B,'LineWidth',2)
hold on
plot(xspan./N_units,SR_weight_A,'LineWidth',2)

title('Different SR gradients required for forward vs. backward motion')
%legend('B-class','A-class')
xlabel('Body length')
ylabel('SR weighting')
pause

%---------------------------------------------
% Save as tikz file
%---------------------------------------------

matlab2tikz('figurehandle',f,'filename','./PaperPlots/gradientplot.tex')


%---------------------------------------------
% Original C code
%---------------------------------------------
% 	for(int i = 0; i < NSEG; ++i){
% NMJ_weight_B[i] =  0.7*(1.0 - i * 0.6/NSEG) *(1.0+(env->NMJ_weight_B_factor)*doubleSigmoid(timenow,i*(env->NMJ_B_vel)+env->NMJ_B_tStart,i*(env->NMJ_B_vel)+env->NMJ_B_tEnd,env->NMJ_B_slope));
% 	}
%    	NMJ_weight_B[0] /= 1.5
% for(int i = 0; i < NSEG; ++i){
% 		NMJ_weight_A[i] =  0.7*(1.0 - (NSEG - i - 1) * 0.6/NSEG) *(1.0+(env->NMJ_weight_A_factor)*doubleSigmoid(timenow,i*(env->NMJ_A_vel)+env->NMJ_A_tStart,i*(env->NMJ_A_vel)+env->NMJ_A_tEnd,env->NMJ_A_slope));
% 	}
%    	NMJ_weight_A[NSEG-1] /= 1.5;
% for(int i = 0; i < N_units; ++i){	
% 		SR_weight_A[i] = (0.65* (1.0-doubleSigmoid(SRpos,i*(env->SR_A_vel)+env->SR_A_tStart,i*(env->SR_A_vel)+env->SR_A_tEnd,env->SR_A_slope)) )*(0.4 + 0.08*(N_units-i-1))*(N_units/12.0)*(2.0/N_seg_per_unit);
% 		SR_weight_B[i] = (0.65* (1.0-doubleSigmoid(SRpos,i*(env->SR_B_vel)+env->SR_B_tStart,i*(env->SR_B_vel)+env->SR_B_tEnd,env->SR_B_slope)) )*(0.4 + 0.08*i)*(N_units/12.0)*(2.0/N_seg_per_unit);		
%    	} 


%==========================================================================


%% Robustness studies with NSR
% Varying 5 environmental variables by +-1%

datEnvNSR=WormMetaDataObj('../Simdata_NSR_robust/',...
{'N_SR'},...
struct('ignoreVar',{{'DURATION','env_vars'}},...
'toPlotOnConstruct',false,'useErrorBars',true));

%Should be 2 classes
datEnvNSR.changePlotVar('COMVelocityMean');
datEnvNSR.plotAllClasses;

f1 = datEnvNSR.runClasses.Class_1.ClassFig_1;
f2 = datEnvNSR.runClasses.Class_2.ClassFig_2;


if toSave
    saveas(f1,'./PaperPlots/robust_NSR_1med_errbar.jpg');
    saveas(f2,'./PaperPlots/robust_NSR_0med_errbar.jpg');

    matlab2tikz('figurehandle',f1,'filename','./PaperPlots/robust_NSR_1med_errbar.tex')
    matlab2tikz('figurehandle',f2,'filename','./PaperPlots/robust_NSR_0med_errbar.tex')
end

%==========================================================================

%% Robustness
%Plotting different turn angle differences when varying internal parameters

datRobust=WormMetaDataObj('../Simdata_OmWave_Jac2/',...
{'SR_B_TSTART'},...
struct('ignoreVar',{{'DURATION','NMJ_B_TEND','NMJ_B_TSTART','SR_B_TEND'}},...
'toPlotOnConstruct',false,'useErrorBars',true,'modulo2pi',false,...
'customClassVars',{{'L_SEG','K_PE','D_PE','AE_PE_ratio','K_DE','D_DE','T_MUSCLE'}}));

%Plot them on the same error plot by hand
datMat = [];
varNames = {};
classNames = fieldnames(datRobust.runClasses);
thisF = 'changeInAngle'; %The field to plot

for j=1:length(classNames)
    thisC = classNames{j};
    thisDat = datRobust.runClasses.(thisC);
    datMat = [datMat; thisDat.(thisF).rawWErrs];
    varNames{j} = thisDat.ClassVar;
end

%Subtract the mean to plot a turning angle deviation
datMat(:,2) = datMat(:,2) - mean(datMat(:,2));

f = figure('DefaultAxesFontSize',14);
hold on;
for j=1:length(varNames)
    errorbar(j,...
        datMat(j,2),datMat(j,3),'MarkerSize',2,'LineWidth',2);
end
ylabel('Deviation of turning angle')
xlabel('Body parameters changed')
xticks(1:length(varNames))
xticklabels(varNames);
%set(gca,'TickLabelInterpreter','none')

if toSave
    saveas(f,'./PaperPlots/robust_env.jpg');

    matlab2tikz('figurehandle',f,'filename','./PaperPlots/robust_env.tex')
end

%==========================================================================

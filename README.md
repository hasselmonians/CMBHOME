## CMBHOME: A Custom MATLAB Class for Neural Data
<center>
Bill Chapman <sup>a</sup>, Andrew Bogaard<br>
Boston University, Department of Brain and Psychological Sciences
</center>

a: e-mail: wchapman@bu.edu

[TOC]

## Introduction:
This documentation for the MATLAB Toolbox, CMBHOME, is meant to be a reference sheet for the standard data structure and its functionatility. For more in depth information about the analysis methods, use the doc function within MATLAB, or view the comment blocks at the top of each .m file included in the Toolbox. Where possible we have also included citations of published work defining the methods for each function.

The overall layout of the Toolbox is:
+ ```+CMBHOME```
    * ```@Session```: Includes time varying behavioral (video) data.
    * ```@Spike```:   Includes timestamps for single-unit activity
    * ```@LFP```:     Incledes time-varying Local Field Potential data.

Throughout this tutorial we call the main CMBHOME Session object "root", however this variable name can be whatever you desire in your implementation. Fields and methods are referenced as root.fieldName or root.method.

For the most basic end-user applications, only ```Session``` objects will be referenced directly. More advanced users or those implementing new analyses may however reference ```Spike``` and ```LFP```.

## Installation
The entire package can be downloaded from GitHub at: github.com/wchapman/CMBHOME. 

Once the files have been downloaded from GitHub, you need to add the folder containing +CMBHOME to your MATLAB path. There are a couple of different methods to do this:
    1. GUI:
        Click   
    2. Commands:
        ```
        addpath('/path/to/local/git')
        ```


You can add this command to your startup.m file to execute each time that MATLAB starts up.

Once the package is in your path, type import CMBHOME.* to import the entire package to your workspace. You are now set up to use the toolbox. It is also safe to update the toolbox via GitHub without breaking your configurations.

## Importing Data

## Using the Toolbox
At it's most basic functionality, CMBHOME aims to organize behavioral/electrophysiological coming from different sources in to a common format. This is achieved during the import process. The imported data is now available as fields (usually vectors) in the ```Session``` object. However, there are additional functionalities to filter this data and apply some common analyses.

### Session State:
There are two properties of a session object which effect the state of the object. These are root.epoch and root.cel. Setting either one of these values changes the **state** of the session object. Once the state is changed, CMBHOME automatically updates the dynamic fields of Session objects. Using the state fields of session is one of the primary functionalities of CMBHOME.

#### Epochs
One of the main functionalities of CMBHOME is to create "epochs", or periods of time, to analyze seprately. 

For instance, you may want to analyze the first and last 100 seconds of the session. To do this, you would use:
```matlab
root.epoch = [0 100; ...
              root.b_ts(end)-100:root.b_ts(end)];
```
At this point, all of the time-variable fields in root would be a 2x1 cell array, with each cell containing data from a single epoch. This data can then be analyzed either seperately using ```cellfun(...)```, or in a continuized manner using ```CMBHOME.Utils.ContinuizeEpochs(...)```.

#### Active Cell
After importing spiking data into a session object, there are multiple listings in ```root.cells```. root.cel_XX fields are dynamically populated based on the epoch and the active cell. So, if you want to get the x spiking locations of the second cluster on tetrode 12, use:
```matlab
root.cel = [12 2];
root.cel_x
```
root.cel_XX fields will always be a cell array of MxN, where M is the number of epochs, and N is the number of cells currently set. 

## Visualize2
To perform many first-pass analyses, we have implemented ```root.Visualize2```. This method allows non-programmers to view some basic information about a cell without using the command line at all. 

Setting root.cel can be achieved from the drop-down menu (A), and changing the epoch can be achieved from the appropriate button (B). The desired analyses can be chosen from the tick box menu (C). Once all of this is set, press the update button (D) to generate a figure containing the results from all of the desired analyses.

|Field                      |Description                                   |
|---------------------------|----------------------------------------------|
|Title Information          |Summary information about the cell and session|
|Rate Map                   |Displays the normalized spatial rate map
|Polar Rate Map             |Displays the head direction tuning map
|Trajectory                 |Shows position (x,y) of the rat in black, and spike locations (red dots ) |
|ISI Histogram              | Histogram of the interspike interval         |
|Spike Time AutoCorrelation |                                              |
|Head Direction             |Plots head direction over time                |
|Spike Raster               |                                              |
|Rate Map Autocorrelogram   |                                              |
|Gridness Autocorrelogram   |                                              |
|Gridness 3 Score           |                                              |
|Gridness 2 Score           |                                              |
|Speed vs Firing Frequency  |Modulation of firing rate as a function of running speed, and summary stats|
|Intrinsic Frequency        |                                              |
|Precession (grid)          |Use only when grid cell suspected             |
|Precession (place)         |Use only when place cell suspected            |
|thetaModGamma              |Theta Modulation of Gamma Power               |
|waveform                   |Plots average waveform on each tetrode        |


## Tutorial
We now provide tutorials for some simple analyses that go beyond the functionality of ```Visualize2```.

### Example 1: Play a video of the session:
This is provided as a simple introduction to accessing the fields of a Session object. Once you load your object and rename it as root, understand and execute the lines below:

```matlab
% Plays a video of the trajectory
root.cel = root.cells(1,:); % watch the first cell
speed = 10; % play back at 10X
lookBack = 20;  % trail 20 seconds

% Fix the position and headdir to make sense
root = root.FixDir;
root = root.FixPos;

pauseTime = (1/root.fs_video);  % how long between each frame
stepSize = speed;

x = root.x; y= root.y; t = root.ts; hd = CMBHOME.Utils.ContinuizeEpochs(root.headdir);

cx = CMBHOME.Utils.ContinuizeEpochs(root.cel_x);    % Get all x positions at firing
cy = CMBHOME.Utils.ContinuizeEpochs(root.cel_y);    % Get all y positions at firing
ci = CMBHOME.Utils.ContinuizeEpochs(root.cel_i);    % Get all indices of firing
f1=figure('Position',[200 200 1024 800]);

lookBack = ceil(lookBack*root.fs_video);

%% Loop:
for i = stepSize+1:stepSize:length(t)
    % figure setup
    figure(f1);cla;hold on
    
    inds = max([i-lookBack, 1]):i;
    cinds = (ci<=i) & (ci>=inds(1));
    
    plot(x(inds),y(inds));
    
    if ~isempty(inds) 
        plot(cx(cinds),cy(cinds),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0]);    
    end
    
    
    % render
    xlim([min(root.x) max(root.x)]);ylim([min(root.y) max(root.y)])
    title(num2str(floor(t(i))))
    spks = numel(intersect((i-stepSize):i,ci))>0;
    
    pause(pauseTime)
    
end
```

clear all; close all;clc;dbstop if error;

CMBHPATH = 'D:\Dropbox\scripts\packages\CMBHOME_stable\';
PIPATH = 'D:\Dropbox\scripts\general_scripts\pass_index';

load('D:\Dropbox\scripts\jrc_scripts\Projects\CaitlinData\CMBobject_clamps-69 (1).mat')
% root = root.FixPos;
% root.spatial_scale=1;
% root.spatial_scale = 150/range(root.x);
%%
root.cel = [15 4];
root.active_lfp = 15;
root.pass_index
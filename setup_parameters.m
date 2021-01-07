% Setup_parameters for ambient noise processing
%
% NJA, 4/2/2016
% JBR, 6/16/2016

addpath('./functions/');
addpath('./functions/calc_Rayleigh_disp/');

%%% --- Paths to important files --- %%%
parameters.workingdir = pwd;
parameters.workingdir = [parameters.workingdir,'/'];

% parameters.datapath = ['/data/irma6/jrussel/YoungPacificORCA/SAC_1Hz/']; %'../nomelt_data_5sta/';
parameters.datapath = ['/Volumes/Russell_2TB/irma6/jrussel/YoungPacificORCA/SAC_50Hz_rmresp/']; %'../nomelt_data_5sta/';

parameters.PZpath = '../INSTRUMENT/';
parameters.ccfpath = '/Volumes/Russell_2TB/irma6/jrussel/YoungPacificORCA/AmbNoise/ccf_50Hz/ccf_raw/';
% parameters.ccfpath = './ccf/ccf_OBNW/';
parameters.figpath = [parameters.workingdir,'figs/'];
parameters.seis_path = [parameters.workingdir,'seismograms/'];
parameters.orientation_path = '/Users/russell/Lamont/PROJ_YoungPacificORCA/ORIENTATIONS/ORCA_orientations.txt';
[stalist, stalat, stalon, staz] = textread(['stations_all.txt'],'%s %f %f %f\n');
parameters.stalist = stalist;
parameters.stalat = stalat;
parameters.stalon = stalon;
parameters.staz = staz;
parameters.nsta = length(parameters.stalist);

%%% --- Parameters to build up gaussian filters --- %%%
parameters.min_width = 0.18;
parameters.max_width = 0.30;

%%% --- Parameters for initial processing --- %%%
parameters.dt = 1/50; %1; % sample rate
parameters.comp = 'BH'; % component
parameters.mindist = 0; % min. distance in kilometers

%%% --- Parameters for ccf_ambnoise --- %%%
parameters.winlength = 3; %hours
parameters.Nstart_sec = 50; % number of sections to offset start of seismogram

%%% --- Parameters for fitbessel --- %%%
parameters.npts = parameters.winlength*3600;

%%% --- Parameters for using Radon Transform picks --- %%%
parameters.path_LRT_picks = './mat-LRTdisp/LRT_picks/';

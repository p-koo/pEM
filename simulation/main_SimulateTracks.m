%-------------------------------------------------------
% This script generates a collection of synthetic particle
% trajectories that have different diffusive states,
% specified by their diffusion coefficient, static localization
% noise.  This simulation does not allow particle tracks to 
% transition between diffusive states within the trajectory. 
%
% Code written by:
% 	Peter Koo
%	Yale University, Department of Physics, New Haven, CT, 06511
%-------------------------------------------------------

clear all;
clc;
close all;

%% simulate tracks

filename = 'case1_500';     % name of output file
numTracks = 500;            % number of tracks to simulate
Dindex = [.01, .3, 1.2, 2.8];  % diffusion coefficients of each state (um^2s^-1)
Pindex = [.2, .3, .4, .1];     % population fraction of each state
Sindex = ones(1,4)*.04;     % static localization noise of each state (um)
dt = .032;                  % time duration between steps


variableLength = 0;     % sets whether particle tracks have variable lengths or not
params.meanN = 25;      % mean length of particle trajectory
if variableLength == 1    
    params.minN = 15;       % minimum length of particle trajectory
    params.maxN = 60;       % maximum length of particle trajectory
else
    params.minN = params.meanN;         % minimum length of particle trajectory
    params.maxN = params.meanN+1;       % maximum length of particle trajectory
end    

% generate particle tracks
simResults = GenerateParticleTracks(Dindex,Pindex,Sindex,numTracks,dt,params);
X = simResults.X;

% since all the tracks start at origin, create a rectangular grid
% to place each track to spatially separate
height = numTracks/20;        % um
width = numTracks/height;     % um
separation = .5;              % um
X_spatial = cell(numTracks,1);
k = 1;
for i = 1:width
    for j = 1:height
        x = X{k}(:,1)+i*separation;
        y = X{k}(:,2)+j*separation;
        X_spatial{k} = [x y];
        k = k + 1;
    end
end
X = X_spatial;

% save tracks in format for pEM
save(filename,'X');







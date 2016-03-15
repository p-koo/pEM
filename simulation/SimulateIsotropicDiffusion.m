function results = SimulateIsotropicDiffusion(N,numSim,D,sigma,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description: This function simulates isotropic brownian motion with
% motion blur and static localization noise.  motion blur is added by
% simulating substeps at 1 ms.
% 
% input:
%     N - number of steps for a given trajectory
%     numSim - number of simulations
%     D - diffusion coefficient to simulate in um^2/s
%     sigma - static localization uncertainty in um, usually around 0.05 
%     dt - time between steps in s (must be larger than .001)
%         
% output:
%     results - structure of the positions (N x numSim), i.e. each simulated
%               track is a column in the matrix
%       results.allTruePositionsX - true x positions including substeps
%       results.allTruePositionsY - true y positions including substeps
%       results.truePositionsX - true x positions at each dt time step
%       results.truePositionsY - true y positions at each dt time step
%       results.motionPositionsX - motion blurred x positions at each dt 
%       results.motionPositionsY - motion blurred y positions at each dt      
%       results.observedPositionsX - noisy x positions at each dt (motion 
%                                    blur + static localization noise)
%       results.observedPositionsY - noisy y positions at each dt      
%
% Written by: Peter Koo
% Mochrie Lab
% Dept. of Physics
% Yale University
% New Haven, CT 06511
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% simulate x and y displacements for 1 ms time steps
numSubSteps = dt*1000;
displacementsX = randn(N*numSubSteps,numSim)*sqrt(2*D*.001);
displacementsY = randn(N*numSubSteps,numSim)*sqrt(2*D*.001);

% get the true x displacements
allTruePositionsX = cumsum(displacementsX);
allTruePositionsY = cumsum(displacementsY);

% find true displacements for dt
truePositionsX = allTruePositionsX(1:numSubSteps:N*numSubSteps,:);
truePositionsY = allTruePositionsY(1:numSubSteps:N*numSubSteps,:);

% average over true positions to get the motion blur 
motionPositionsX = zeros(N,numSim);
motionPositionsY = zeros(N,numSim);
for i = 1:N
    range = i*numSubSteps-numSubSteps+1:i*numSubSteps;
    motionPositionsX(i,:) = mean(allTruePositionsX(range,:));
    motionPositionsY(i,:) = mean(allTruePositionsY(range,:));
end

% get the displacements for the blurred positions
motionDisplacementX = [diff(motionPositionsX)];
motionDisplacementY = [diff(motionPositionsY)];

% simulate noise 
noiseDisplacementX = randn(N,numSim)*(sigma);
noiseDisplacementY = randn(N,numSim)*(sigma);
noiseDisplacementX = noiseDisplacementX - repmat(mean(noiseDisplacementX),N,1);
noiseDisplacementY = noiseDisplacementY - repmat(mean(noiseDisplacementY),N,1);

% get the observed displacements
observedPositionsX = motionPositionsX + noiseDisplacementX;
observedPositionsY = motionPositionsY + noiseDisplacementY;

% store results as a struct
results.observedPositionsX = observedPositionsX;
results.observedPositionsY = observedPositionsY;
results.motionPositionsX = motionPositionsX;
results.motionPositionsY = motionPositionsY;
results.truePositionsX = truePositionsX;
results.truePositionsY = truePositionsY;
results.allTruePositionsX = allTruePositionsX;
results.allTruePositionsY = allTruePositionsY;








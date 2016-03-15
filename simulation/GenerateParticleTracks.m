function results = GenerateParticleTracks(Dindex,Pindex,Sindex,totSim,dt,params)

% parameters for exponentially distributed random track lengths
meanN = params.meanN;
minN = params.minN;
maxN = params.maxN;

k = 1;
X = {}; X_true = {}; deltaX = {}; deltaX_true = {}; 
for i = 1:length(Dindex)
    D = Dindex(i);
    locNoise = Sindex(i);
    numSim = round(totSim*Pindex(i));

    % get random track length
    N = round(RandomTrackLength(numSim,meanN,minN,maxN));
    for j = 1:length(N)
        
        results = SimulateIsotropicDiffusion(N(j),1,D,locNoise,dt);
        X{k} = [results.observedPositionsX results.observedPositionsY];
        deltaX{k} = diff(X{k});

        X_true{k} = [results.truePositionsX results.truePositionsY];
        deltaX_true{k} = diff(X_true{k});

        k = k + 1;
    end
    
end

results.X = X;
results.deltaX = deltaX;
results.X2 = X_true;
results.deltaX_true = deltaX_true;

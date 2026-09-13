function [] = ParLoop(workers)

currentPool = gcp('nocreate'); % Check if a parallel pool already exists

if isempty(currentPool)
    % If no pool exists, create one
    myCluster = parcluster('local');
    myCluster.NumWorkers = workers; % Set the desired number of workers
    parpool(myCluster, myCluster.NumWorkers);
else
    % A pool already exists; do nothing
    fprintf('A parallel pool is already active.\n');
end

end


% Time gpops
clc;
t0 = tic;
for nphases = 2 : 10
    fprintf('(nphase, T) = (%d, %1.1f)\n', nphases, 2.5);
    try
        TimeGpops_helper( nphases, 2.5 );
    end
end
t1 = toc(t0);

t0 = tic;
for nphases = 2 : 10
    fprintf('(nphase, T) = (%d, %d)\n', nphases, 6);
    try
        TimeGpops_helper( nphases, 6 );
    end
end
t2 = toc(t0);
function logP = logP_scan_topo(P, L, K, A, B, C, E, BA, Q)
%logP = logP_scan_topo(P, L, K, A, B, C, E, BA, Q)
%
%Scans the parameter space defined by the vectors K, B, C, E, BA, Q to
%find the log probability distribution for a topological model.
%
%K: Number of interacting nearest neighbours
%A: Inertia (set to 1)
%B: Alignment
%C: Attraction
%E: Angular noise
%BA: Blind Angle
%Q: Update rate

%Richard Mann (2011)

logP = zeros(length(K), length(B), length(C), length(E), length(BA), length(Q), 'single');


parfor k = 1:length(K)
    lp = zeros(length(B), length(C), length(E), length(BA), length(Q), 'single');
    Kk = K(k);
    for e = 1:length(E)
        Ee = E(e);
        for ba= 1:length(BA)
            BAba = BA(ba);
            for q = 1:length(Q)
                Qq = Q(q);
                lp(:, :, e, ba, q) = logP_sppABC_topo(P, L, Kk, A, B, C, Ee, BAba, Qq);
                
            end
        end
    end
    logP(k, :, :, :, :, :) = lp;
end

logP = squeeze(logP);

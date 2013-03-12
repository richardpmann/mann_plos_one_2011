function logP = logP_scan(P, L, R, A, B, C, E, BA, Q)
%logP = logP_scan(P, L, R, A, B, C, E, BA, Q)
%
%Scans the parameter space defined by the vectors R, B, C, E, BA, Q to
%find the log probability distribution for a geometrical model.
%
%R: interaction radius
%A: Inertia (set to 1)
%B: Alignment
%C: Attraction
%E: Angular noise
%BA: Blind Angle
%Q: Update rate

%Richard Mann (2011)
logP = zeros(length(R), length(B), length(C), length(E), length(BA), length(Q), 'single');


parfor r = 1:length(R) %Interaction radius
    lp = zeros(length(B), length(C), length(E), length(BA), length(Q), 'single');
    Rr = R(r);
    for e = 1:length(E)
        Ee = E(e);
        for ba= 1:length(BA)
            BAba = BA(ba);
            for q = 1:length(Q)
                Qq = Q(q);
                lp(:, :, e, ba, q) = logP_sppABC_euclid(P, L, Rr, A, B, C, Ee, BAba, Qq);
                
            end
        end
    end
    logP(r, :, :, :, :, :) = lp;
end

logP = squeeze(logP);

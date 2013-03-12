%Script to calculate a huge bunch of results needed for paper: Bayesian
%inference for identifying interaction rules in moving animal groups

clear all
tic

%set some parameters for the simulations
iterations = 5;
num_t_step = 10;
num_particles = 25;
Lsize = 10;
velocity = 1;

%Choose sections to do

flag1 = 1;
flag2 = 1;
flag3 = 1;
flag4 = 1;



R_true = 4;
B_true = 0.1;
C_true = 1;

E_true = 0.05*pi;
BA_true = pi/6;

Q_true = 1;


R = linspace(1, 5, 32); %make num R a multiple of 8 for parforing
K = 1:24; %make K a multiple of 8
B = linspace(0.0001, 3, 100);
C = linspace(0.0001, 3, 101);
E = pi*linspace(0.01, 0.1, 15);
BA = pi*linspace(0.1, pi, 21);


if flag1
    
    %Set aside some memory
    meanR = zeros(num_t_step, iterations);
    meanB = zeros(num_t_step, iterations);
    meanC = zeros(num_t_step, iterations);
    meanE = zeros(num_t_step, iterations);
    meanBA = zeros(num_t_step, iterations);
    Pentropy = zeros(num_t_step, iterations);
    BF = zeros(num_t_step, iterations);
    BFB = zeros(num_t_step, iterations);
    
    stdR = zeros(num_t_step, iterations);
    stdB = zeros(num_t_step, iterations);
    stdC = zeros(num_t_step, iterations);
    stdE = zeros(num_t_step, iterations);
    stdBA = zeros(num_t_step, iterations);
    
    meanR2 = zeros(num_t_step, iterations);
    meanB2 = zeros(num_t_step, iterations);
    meanC2 = zeros(num_t_step, iterations);
    meanE2 = zeros(num_t_step, iterations);
    meanBA2= zeros(num_t_step, iterations);
    Pentropy2 = zeros(num_t_step, iterations);
    BF2 = zeros(num_t_step, iterations);
    BFB2 = zeros(num_t_step, iterations);
    
    stdR2 = zeros(num_t_step, iterations);
    stdB2 = zeros(num_t_step, iterations);
    stdC2 = zeros(num_t_step, iterations);
    stdE2 = zeros(num_t_step, iterations);
    stdBA2 = zeros(num_t_step, iterations);
    
    logP = zeros(num_t_step, length(R), length(B), length(C), length(E), length(BA), 'single');
    logP2 = logP;
    logP_topo = logP;
    logP_topo2 = logP;
    
    for it = 1:iterations
        
        %simulate some data
        P = sppABC(Lsize, num_particles, R_true, velocity, 100, 1, B_true, C_true, E_true, BA_true, Q_true, 0);
        
        
        
        %perform inference over varying number of timesteps
        for i = 1:num_t_step
            
            logP(i, :, :, :, :, :) = logP_scan(P(:, :, i:i+1), Lsize, R, 1, B, C, E, BA, 1);%initial geo
            logP2(i, :, :, :, :, :) = logP_scan(P(:, :, 89+i:89+i+1), Lsize, R, 1, B, C, E, BA, 1);%steady-state geo
            logP_topo(i, :, :, :, :, :) = logP_scan_topo(P(:, :, i:i+1), Lsize, K, 1, B, C, E, BA, 1); %initial topo
            logP_topo2(i, :, :, :, :, :) = logP_scan_topo(P(:, :, 89+i:89+i+1), Lsize, K, 1, B, C, E, BA, 1);%stead-state topo
        end
        
        
        %Analysis of first 10 timesteps
        %Find marginal distributions for each parameter
        for i = 1:size(logP, 1)
            L = squeeze(lowexp(sum(logP(1:i, :, :, :, :, :), 1))); L = L/sum(L(:));
            
            L2 = squeeze(sum(sum(sum(sum(L, 2), 3), 4), 5));L2 = L2/sum(L2);
            % L2 = L2(:)./R(:); L2 = L2/sum(L2); %jeffreys prior
            [meanR(i, it), stdR(i, it)] = credible_interval(R, L2);
            
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 3), 4), 5));L2 = L2/sum(L2);
            % L2 = L2(:)./B(:); L2 = L2/sum(L2); %jeffreys prior
            [meanB(i, it), stdB(i, it)] = credible_interval(B, L2);
            
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 4), 5));L2 = L2/sum(L2);
            % L2 = L2(:)./C(:); L2 = L2/sum(L2); %jeffreys prior
            [meanC(i, it), stdC(i, it)] = credible_interval(C, L2);
            
            L2 = squeeze(sum(sum(sum(sum(L, 1),2), 3), 5));L2 = L2/sum(L2);
            [meanE(i, it), stdE(i, it)] = credible_interval(E, L2);
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 3), 4));L2 = L2/sum(L2);
            [meanBA(i, it), stdBA(i, it)] = credible_interval(BA, L2);
            
            Pentropy(i, it) = -1*sum(L(L(:)>0).*log(L(L(:)>0)));
            BF(i, it) = logmeanexp(sum(logP(1:i, :, :, :, :, :), 1)) - logmeanexp(sum(logP_topo(1:i, :, :, :, :, :), 1)); %Bayes factor geo v topo
            
            BFB(i, it) = logmeanexp(sum(logP(1:i, :, :, :, :, :), 1)) - logmeanexp(sum(logP(1:i, :, 1, :, :, :), 1)); %Bayes factor B v no B
            
        end
        
        %Analysis of final 10 timesteps
        %Find marginal distributions for each parameter
        for i = 1:size(logP2, 1)
            L = squeeze(lowexp(sum(logP2(1:i, :, :, :, :, :), 1))); L = L/sum(L(:));
            
            L2 = squeeze(sum(sum(sum(sum(L, 2), 3), 4), 5));L2 = L2/sum(L2);
            %  L2 = L2(:)./R(:); L2 = L2/sum(L2); %jeffreys prior
            [meanR2(i, it), stdR2(i, it)] = credible_interval(R, L2);
            
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 3), 4), 5));L2 = L2/sum(L2);
            %L2 = L2(:)./B(:); L2 = L2/sum(L2); %jeffreys prior
            [meanB2(i, it), stdB2(i, it)] = credible_interval(B, L2);
            
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 4), 5));L2 = L2/sum(L2);
            % L2 = L2(:)./C(:); L2 = L2/sum(L2); %jeffreys prior
            [meanC2(i, it), stdC2(i, it)] = credible_interval(C, L2);
            
            L2 = squeeze(sum(sum(sum(sum(L, 1),2), 3), 5));L2 = L2/sum(L2);
            [meanE2(i, it), stdE2(i, it)] = credible_interval(E, L2);
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 3), 4));L2 = L2/sum(L2);
            [meanBA2(i, it), stdBA2(i, it)] = credible_interval(BA, L2);
            
            
            Pentropy2(i, it) = -1*sum(L(L(:)>0).*log(L(L(:)>0)));
            BF2(i, it) = logmeanexp(sum(logP2(1:i, :, :, :, :, :), 1)) - logmeanexp(sum(logP_topo2(1:i, :, :, :, :, :), 1));
            BFB2(i, it) = logmeanexp(sum(logP2(1:i, :, :, :, :, :), 1)) - logmeanexp(sum(logP2(1:i, :, 1, :, :, :), 1));
        end
        
        
        toc
        it
    end
    
    clear logP* L L2
    disp('Section 1 complete')
    save HUGEBENCHMARK
    
end

if flag2
    
    %Tests on changing angular noise level
    
    %set aside some memory
    
    num_E_test = 20;
    E_test = pi*linspace(0.01, 1, num_E_test);
    
    meanR_E = zeros(num_E_test, iterations);
    meanB_E = zeros(num_E_test, iterations);
    meanC_E = zeros(num_E_test, iterations);
    meanE_E = zeros(num_E_test, iterations);
    meanBA_E = zeros(num_E_test, iterations);
    
    stdR_E = zeros(num_E_test, iterations);
    stdB_E = zeros(num_E_test, iterations);
    stdC_E = zeros(num_E_test, iterations);
    stdE_E = zeros(num_E_test, iterations);
    stdBA_E = zeros(num_E_test, iterations);
    Pentropy_E = zeros(num_E_test, iterations);
    
    
    for it = 1:iterations
        
        for i = 1:length(E_test) %vary angular noise
            %simulate data
            P = sppABC(Lsize, num_particles, R_true, velocity, 100, 1, B_true, C_true, E_test(i), BA_true, Q_true, 0);
            
            %perform inference on first 10 timesteps
            logP = logP_scan(P(:, :, 1:10), Lsize, R, 1, B, C, E_test, BA, 1);
            L = single(lowexp(logP));
            L = L/sum(L(:));
            clear logP
            
            %find marginal distributions
            L2 = squeeze(sum(sum(sum(sum(L, 2), 3), 4), 5));L2 = L2/sum(L2);
            meanR_E(i, it) = sum(L2(:).*R(:));
            stdR_E(i, it) = sqrt(sum(L2(:).*(R(:)-meanR_E(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 3), 4), 5));L2 = L2/sum(L2);
            meanB_E(i, it) = sum(L2(:).*B(:));
            stdB_E(i, it) = sqrt(sum(L2(:).*(B(:)-meanB_E(i)).^2));
            
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 4), 5));L2 = L2/sum(L2);
            meanC_E(i, it) = sum(L2(:).*C(:));
            stdC_E(i, it) = sqrt(sum(L2(:).*(C(:)-meanC_E(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 3), 5));L2 = L2/sum(L2);
            meanE_E(i, it) = sum(L2(:).*E_test(:));
            stdE_E(i, it) = sqrt(sum(L2(:).*(E_test(:)-meanE_E(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 3), 4));L2 = L2/sum(L2);
            meanBA_E(i, it) = sum(L2(:).*BA(:));
            stdBA_E(i, it) = sqrt(sum(L2(:).*(BA(:)-meanBA_E(i)).^2));
            
            Pentropy_E(i, it) = -1*sum(L(L(:)>0).*log(L(L(:)>0)));
            
        end
        toc
        it
    end
    
    clear logP L L2
    disp('Section 2 complete')
    save HUGEBENCHMARK
    
end


if flag3
    %tests on changing update rate, inference assuming Q=1
    
    num_Q_test = 10;
    Q = linspace(0.1, 1, num_Q_test);
    
    %set aside some memory
    meanR_Q = zeros(num_Q_test, iterations);
    meanB_Q = zeros(num_Q_test, iterations);
    meanC_Q = zeros(num_Q_test, iterations);
    meanE_Q = zeros(num_Q_test, iterations);
    meanBA_Q = zeros(num_Q_test, iterations);
    
    stdR_Q = zeros(num_Q_test, iterations);
    stdB_Q = zeros(num_Q_test, iterations);
    stdC_Q = zeros(num_Q_test, iterations);
    stdE_Q = zeros(num_Q_test, iterations);
    stdBA_Q = zeros(num_Q_test, iterations);
    Pentropy_Q = zeros(num_Q_test, iterations);
    
    
    for it = 1:iterations
        for i = 1:length(Q) %vary update rate
            %simulate data
            P = sppABC(Lsize, num_particles, R_true, velocity, 100, 1, B_true, C_true, E_true, BA_true, Q(i), 0);
            %perform inference on first 10 timesteps
            logP = logP_scan(P(:, :, 1:10), Lsize, R, 1, B, C, E, BA, 1);
            
            L = single(lowexp(logP));
            clear logP
            L = L/sum(L(:));
            
            %find marginal distributions
            L2 = squeeze(sum(sum(sum(sum(L, 2), 3), 4), 5));L2 = L2/sum(L2);
            meanR_Q(i, it) = sum(L2(:).*R(:));
            stdR_Q(i, it) = sqrt(sum(L2(:).*(R(:)-meanR_Q(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 3), 4), 5));L2 = L2/sum(L2);
            meanB_Q(i, it) = sum(L2(:).*B(:));
            stdB_Q(i, it) = sqrt(sum(L2(:).*(B(:)-meanB_Q(i)).^2));
            
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 4), 5));L2 = L2/sum(L2);
            meanC_Q(i, it) = sum(L2(:).*C(:));
            stdC_Q(i, it) = sqrt(sum(L2(:).*(C(:)-meanC_Q(i)).^2));
            
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 3), 5));L2 = L2/sum(L2);
            meanE_Q(i, it) = sum(L2(:).*E(:));
            stdE_Q(i, it) = sqrt(sum(L2(:).*(E(:)-meanE_Q(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(L, 1), 2), 3), 4));L2 = L2/sum(L2);
            meanBA_Q(i, it) = sum(L2(:).*BA(:));
            stdBA_Q(i, it) = sqrt(sum(L2(:).*(BA(:)-meanBA_Q(i)).^2));
            
            Pentropy_Q(i, it) = -1*sum(L(L(:)>0).*log(L(L(:)>0)));
            
        end
        toc
        it
    end
    disp('Section 3 complete')
    
    clear logP L L2
    
    save HUGEBENCHMARK
    
end



if flag4
    %tests on changing update rate, inference allows varying Q
    num_Q_test = 10;
    Q = linspace(0.1, 1, num_Q_test);
    iterations = 5;
    
    %set aside some memory
    meanR_Qfree = zeros(num_Q_test, iterations);
    meanB_Qfree = zeros(num_Q_test, iterations);
    meanC_Qfree = zeros(num_Q_test, iterations);
    meanE_Qfree = zeros(num_Q_test, iterations);
    meanBA_Qfree = zeros(num_Q_test, iterations);
    meanQ_Qfree = zeros(num_Q_test, iterations);
    
    stdR_Qfree = zeros(num_Q_test, iterations);
    stdB_Qfree = zeros(num_Q_test, iterations);
    stdC_Qfree = zeros(num_Q_test, iterations);
    stdE_Qfree = zeros(num_Q_test, iterations);
    stdBA_Qfree = zeros(num_Q_test, iterations);
    stdQ_Qfree = zeros(num_Q_test, iterations);
    Pentropy_Qfree = zeros(num_Q_test, iterations);
    
    for it = 1:iterations
        for i = 1:length(Q) %vary update rate
            
            %simulate some data
            P = sppABC(Lsize, num_particles, R_true, velocity, 100, 1, B_true, C_true, E_true, BA_true, Q(i), 0);
            
            %perform inference on first 10 timesteps
            logP = logP_scan(P(:, :, 1:10), Lsize, R, 1, B, C, E, BA, Q);
            
            L = single(lowexp(logP));
            clear logP
            L = L/sum(L(:));
            
            %find marginal distributions
            L2 = squeeze(sum(sum(sum(sum(sum(L, 2), 3), 4), 5), 6));L2 = L2/sum(L2);
            meanR_Qfree(i, it) = sum(L2(:).*R(:));
            stdR_Qfree(i, it) = sqrt(sum(L2(:).*(R(:)-meanR_Qfree(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(sum(L, 1), 3), 4), 5), 6));L2 = L2/sum(L2);
            meanB_Qfree(i, it) = sum(L2(:).*B(:));
            stdB_Qfree(i, it) = sqrt(sum(L2(:).*(B(:)-meanB_Qfree(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(sum(L, 1), 2), 4), 5), 6));L2 = L2/sum(L2);
            meanC_Qfree(i, it) = sum(L2(:).*C(:));
            stdC_Qfree(i, it) = sqrt(sum(L2(:).*(C(:)-meanC_Qfree(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(sum(L, 1), 2), 3), 5), 6));L2 = L2/sum(L2);
            meanE_Qfree(i, it) = sum(L2(:).*E(:));
            stdE_Qfree(i, it) = sqrt(sum(L2(:).*(E(:)-meanE_Qfree(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(sum(L, 1), 2), 3), 4), 6));L2 = L2/sum(L2);
            meanBA_Qfree(i, it) = sum(L2(:).*BA(:));
            stdBA_Qfree(i, it) = sqrt(sum(L2(:).*(BA(:)-meanBA_Qfree(i)).^2));
            
            L2 = squeeze(sum(sum(sum(sum(sum(L, 1), 2), 3), 4), 5));L2 = L2/sum(L2);
            meanQ_Qfree(i, it) = sum(L2(:).*Q(:));
            stdQ_Qfree(i, it) = sqrt(sum(L2(:).*(Q(:)-meanQ_Qfree(i)).^2));
            
            Pentropy_Qfree(i, it) = -1*sum(L(L(:)>0).*log(L(L(:)>0)));
            
        end
        toc
        it
    end
    disp('Section 4 complete')
    
    clear logP L L2
    
    save HUGEBENCHMARK
    
end

 
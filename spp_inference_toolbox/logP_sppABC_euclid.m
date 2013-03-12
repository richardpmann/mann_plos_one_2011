function logP = logP_sppABC_euclid(P,L,R,a,b,c,E, ba, p)
%logP = logP_sppABC_euclid(P,L,R,a,b,c,E, ba, p)
%
%Gives the log prob. of some data based on a geometrical SPP model
%
%L=side length of square
%N=number of particles
%R=interaction radius
%a,b,c weights
%E=error const
%ba= blind angle
%p = prob. of updating in time step.

%Richard Mann (2010)


N = size(P, 1); %Number of particles
t = size(P, 3); %Number of time steps


[B, C] = ndgrid(b, c); %mesh of B and C weights


logP = zeros(size(B)); %Prepare memory

for k=1:t-1 % For each time step
    
    for i=1:N % Go through every particle.
        
        %nh function calculates (for particle i)
        %1. its neighborhood (i.e. finds the particles in its neighborhood.)
        %2. the direction from it toward the center of mass of its neighborhood (CM)
        %3. the 'mean' direction of its neighbors (AL)
        
        CMAL=nhBlind3x(i,squeeze(P(:, :, k)),R,L, ba);
        
        CM=CMAL(1,:); % Direction toward the center of mass of the neighbors
        AL=CMAL(2,:); % The 'mean' direction of the neighbors
        n=CMAL(3,1); % Number of non-self neighbors
        
        %CM and AL are 0 if n = 0;
        
        
        %Calculate the predicted new direction if updating is done
        
        Dir_mean1=a*cos(P(i,3, k))+B.*AL(1)+C.*CM(1); %New direction of particle i
        Dir_mean2=a*sin(P(i,3, k))+B.*AL(2)+C.*CM(2);
        
        %The pdf is now a mixture of the updated direction and the
        %previous, weighted by the probability of updating
        
        lp = log(p) + log_circ_vmpdf(atan2(Dir_mean2, Dir_mean1), P(i, 3, k+1), 1/E^2);
        lq = log(1-p) + log_circ_vmpdf(P(i, 3, k), P(i, 3, k+1), 1/E^2); %prob for no action taken
        
        %%prob for no action taken
        m = max(lp, lq);
        lp_step = m + log(exp(lp-m) + exp(lq-m));
        
        logP = logP + lp_step; %add the log-prob from this time step and particle
        
        
        
    end
    
    
    
end

logP = squeeze(logP); %remove extraneous dimensions

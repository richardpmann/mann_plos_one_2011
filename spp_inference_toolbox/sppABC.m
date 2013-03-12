
function P = sppABC(L,N,R,velocity,t,a,b,c,E, ba, p, plot_flag)
%P = sppABC(L,N,R,velocity,t,a,b,c,E, ba, p, plot_flag)
%
%A generative SPP model based on a function from Daniel Strombom
%
%L=side length of square
%N=number of particles
%R=interaction radius
%velocity=speed per time step
%t=number of time steps.
%a,b,c weights
%E=error const
%ba= blind angle
%p = prob. of updating in time step.
%
%P = output positions and angular headings

%Richard Mann & Daniel Strombom (2010)

if ~exist('plot_flag')
    plot_flag =0;
end

% Create population record and initiate positions and directions
P=zeros(N,3, t+1);
P(:, 1:2, 1) = rand(N, 2)*L;
P(:, 3, 1) = rand(N, 1)*2*pi -pi;

for k=1:t % For each time step
    
    if plot_flag
        plot(P(:,1, k),P(:,2, k),'k.','markersize',10);
        hold on
        for r=1:N
            plot([P(r,1, k),P(r,1, k)+0.5*cos(P(r,3, k))],[P(r,2, k),P(r,2, k)+0.5*sin(P(r,3, k))],'r-');
            
        end
        
        hold off
        axis([0 L 0 L]);
        xlabel('X position')
        ylabel('Y position')
        axis manual
        drawnow
        
        
        
    end
    
    rand_choice = rand(N, 1) < p; %Choose particles to update this timestep
    rand_choice = find(rand_choice);
    P(:, 3, k+1) = P(:, 3, k)+circ_vmrnd(0, 1/E^2, N); %will replace if a decision is made
    NormDir = [cos(P(:, 3, k+1)), sin(P(:, 3, k+1))];
    if isempty(rand_choice)
        %nothing happens
    else
        
        for ri= 1:length(rand_choice) % Go through every particle.
            i = rand_choice(ri);
            %nh function calculates (for particle i)
            %1. its neighborhood (i.e. finds the particles in its neighborhood.)
            %2. the direction from it toward the center of mass of its neighborhood (CM)
            %3. the 'mean' direction of its neighbors (AL)
            
            %Choose either a metric or topological neighbourhood
            
            %metric
            CMAL=nhBlind3x(i,squeeze(P(:, :, k)),R,L, ba); % CMAL=[dir toward CM, 'mean' dir of Neighbors]
            
            %topological
            %CMAL=nhTopoBlind3x(i,squeeze(P(:, :, k)),R,L, ba); %R here is
            %the number of nearest neighbours
            
            CM=CMAL(1,:); % Direction toward the center of mass of the neighbors
            AL=CMAL(2,:); % The 'mean' direction of the neighbors
            n=CMAL(3,1); % Number of non-self neighbors
            
            %Calculate new headings and update
            if n==0
                Dir=a*[cos(P(i,3, k)),sin(P(i,3, k))]; %New direction of particle i
                NormDir(i, :)=Dir/sqrt(sum(Dir.^2));
            else
                Dir=a*[cos(P(i,3, k)),sin(P(i,3, k))]+b*AL+c*CM; %New direction of particle i
                NormDir(i, :)=Dir/sqrt(sum(Dir.^2)); %Normalized direction of particle i
            end
            
            P(i,3, k+1)=atan2(NormDir(i,2),NormDir(i,1)) + circ_vmrnd(0, 1/E^2); %New directional angle, inc. von-mises noise
            
            
        end
        
    end
    
    %Add small random variation to the speeds
    rv = 0.2*rand(N, 1);
    
    %Update new positions
    P(:,1, k+1)=mod(P(:,1, k)+(1+rv).*velocity.*NormDir(:,1),L); %New x-coordinate
    P(:,2, k+1)=mod(P(:,2, k)+(1+rv).*velocity.*NormDir(:,2),L); %New y-coordinate
    
    
end

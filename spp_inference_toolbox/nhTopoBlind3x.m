function [CMAL, NE]=nhTopoBlind3x(i,P,k, L,b)
%A topological neighbourhood function to determine the response to neighbours based on
%k nearest neighbour particles

%Richard Mann (2010)
ba = b/2; %for consistency with daniels

AG=zeros(k,3);
v=[-cos(P(i,3)),-sin(P(i,3))]; %opposite to particle i:s direction of travel

%First find distance to every other agent.

d = sqrt(sum(bsxfun(@minus, P(:, 1:2), P(i, 1:2)).^2, 2)); %normal distance
u=[P(:,1)-P(i,1),P(:,2)-P(i,2)];
theta=acos((v(1)*u(:,1)+v(2)*u(:, 2))./(sqrt(u(:,1).^2+u(:,2).^2))); 
d(abs(theta)<ba, 1) = NaN; 

[~, sortidx] = sort(d, 'ascend');

k = min(k, length(sortidx)-1); %if k is more than total number of other particles, limit it.cd

NNidx = sortidx(2:k+1); %the 1st point will be the animal itself


for count = 1:k
AG(count,:)=[P(NNidx(count),1),P(NNidx(count),2),P(NNidx(count),3)];
end


if size(AG,1)>0
    
    K=atan2(sum(sin(AG(:,3))),sum(cos(AG(:,3)))); %Mean angle of neighbors
    AL=[cos(K),sin(K)]; %'Mean' direction of neighbors
    
    CMx=mean(AG(:,1));
    CMy=mean(AG(:,2));
    
    if CMx==P(i,1) && CMy==P(i,2) %!!! I verkligheten b?r e: kunna bli, men antar avrundningsfel f?r att det blir s? el dylikt, och v?l:er d? forts?tt i tidigare riktning
        CMAL=[0,0;0,0;0,0];
    else
        CM=[CMx-P(i,1),CMy-P(i,2)]; %Direction toward center of mass of neighbors
        CM=(1/((CM(1,1)^2+CM(1,2)^2)^(1/2)))*CM; %Normed direction toward center of mass   
    
        CMAL=[CM;AL;size(AG,1),0];
    end
else
    CMAL=[0,0;0,0;0,0];
end

NE = NNidx;




function CMAL=nhBlind3x(i,P,R,L,b)
%A metric neighbourhood function to determine the response to neighbours based on
%an interactio  radius. 

%Daniel Strombom (2010) adapted by Richard Mann (2010)
ba=b/2;
N=size(P,1);
%NE=zeros(N,2);
AG=zeros(N,3);
count=0;
v=[-cos(P(i,3)),-sin(P(i,3))]; %opposite to particle i:s direction of travel


for j=1:size(P,1)
    if j~=i
        if (P(i,1)-P(j,1))^2+(P(i,2)-P(j,2))^2<R^2
            u=[P(j,1)-P(i,1),P(j,2)-P(i,2)];
            theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
            if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
                count=count+1;
                %NE(count,:)=[j,1];
                AG(count,:)=[P(j,1),P(j,2),P(j,3)];
            end
%         elseif (P(i,1)-(P(j,1)-L))^2+(P(i,2)-P(j,2))^2<R^2
%             u=[(P(j,1)-L)-P(i,1),P(j,2)-P(i,2)];
%             theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
%             if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
%                 count=count+1;
%                 %NE(count,:)=[j,2];
%                 AG(count,:)=[P(j,1)-L,P(j,2),P(j,3)];
%             end
%         elseif (P(i,1)-(P(j,1)-L))^2+(P(i,2)-(P(j,2)-L))^2<R^2
%             u=[(P(j,1)-L)-P(i,1),(P(j,2)-L)-P(i,2)];
%             theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
%             if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
%                 count=count+1;
%                 %NE(count,:)=[j,3];
%                 AG(count,:)=[P(j,1)-L,P(j,2)-L,P(j,3)];
%             end
%         elseif (P(i,1)-P(j,1))^2+(P(i,2)-(P(j,2)-L))^2<R^2
%             u=[P(j,1)-P(i,1),(P(j,2)-L)-P(i,2)];
%             theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
%             if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
%                 count=count+1;
%                 %NE(count,:)=[j,4];
%                 AG(count,:)=[P(j,1),P(j,2)-L,P(j,3)];
%             end
%         elseif (P(i,1)-(P(j,1)+L))^2+(P(i,2)-P(j,2))^2<R^2
%             u=[(P(j,1)+L)-P(i,1),P(j,2)-P(i,2)];
%             theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
%             if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
%                 count=count+1;
%                 %NE(count,:)=[j,5];
%                 AG(count,:)=[P(j,1)+L,P(j,2),P(j,3)];
%             end
%         elseif (P(i,1)-(P(j,1)+L))^2+(P(i,2)-(P(j,2)+L))^2<R^2
%             u=[(P(j,1)+L)-P(i,1),(P(j,2)+L)-P(i,2)];
%             theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
%             if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
%                 count=count+1;
%                 %NE(count,:)=[j,6];
%                 AG(count,:)=[P(j,1)+L,P(j,2)+L,P(j,3)];
%             end
%         elseif (P(i,1)-P(j,1))^2+(P(i,2)-(P(j,2)+L))^2<R^2
%             u=[P(j,1)-P(i,1),(P(j,2)+L)-P(i,2)];
%             theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
%             if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
%                 count=count+1;
%                 %NE(count,:)=[j,7];
%                 AG(count,:)=[P(j,1),P(j,2)+L,P(j,3)];
%             end
%         elseif (P(i,1)-(P(j,1)-L))^2+(P(i,2)-(P(j,2)+L))^2<R^2
%             u=[(P(j,1)-L)-P(i,1),(P(j,2)+L)-P(i,2)];
%             theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
%             if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
%                 count=count+1;
%                 %NE(count,:)=[j,8];
%                 AG(count,:)=[P(j,1)-L,P(j,2)+L,P(j,3)];
%             end
%         elseif (P(i,1)-(P(j,1)+L))^2+(P(i,2)-(P(j,2)-L))^2<R^2
%             u=[(P(j,1)+L)-P(i,1),(P(j,2)-L)-P(i,2)];
%             theta=acos((v(1)*u(1)+v(2)*u(2))/(sqrt(u(1)^2+u(2)^2))); %Angle between opposite i-direction and direction from i toward j.
%             if abs(theta)>ba %If particle j not in i:s blind angle add as neighbor
%                 count=count+1;
%                 %NE(count,:)=[j,9];
%                 AG(count,:)=[P(j,1)+L,P(j,2)-L,P(j,3)];
%             end
%         else
%             %NE=NE;
        end
    end
end

%NE=NE(1:count,:);
AG=AG(1:count,:);

if size(AG,1)>0
    
    K=atan2(sum(sin(AG(:,3))),sum(cos(AG(:,3)))); %Mean angle of neighbors
    AL=[cos(K),sin(K)]; %'Mean' direction of neighbors
    
    CMx=mean(AG(:,1));
    CMy=mean(AG(:,2));
    
    if CMx==P(i,1) && CMy==P(i,2) %!!! I verkligheten b?r ej kunna bli, men antar avrundningsfel f?r att det blir s? el dylikt, och v?ljer d? forts?tt i tidigare riktning
        CMAL=[0,0;0,0;0,0];
    else
        CM=[CMx-P(i,1),CMy-P(i,2)]; %Direction toward center of mass of neighbors
        CM=(1/((CM(1,1)^2+CM(1,2)^2)^(1/2)))*CM; %Normed direction toward center of mass   
    
        CMAL=[CM;AL;size(AG,1),0];
    end
else
    CMAL=[0,0;0,0;0,0];
end
%TESTS

%Neighborhood
%hold on
%plot(P(:,1),P(:,2),'.'); %All particles 
%plot(P(i,1),P(i,2),'rv'); %Particle i
%for k=1:size(NE,1)
 %   plot(P(NE(k),1),P(NE(k),2),'o'); %Neighbors of i marked with ring
 %   plot([P(NE(k),1),P(NE(k),1)+cos(P(NE(k),3))],[P(NE(k),2),P(NE(k),2)+sin(P(NE(k),3))],'k-'); %Direction of travel of neighbors
%end

%Center of mass (point)
%plot(mod(CMx+L,L),mod(CMy+L,L),'*'); %plotting center of mass

%Direction from individual i toward center of mass.
%plot(CM(1,1),CM(1,2),'x');

%Mean direction of neighbors
%plot(AL(1,1),AL(1,2),'+');


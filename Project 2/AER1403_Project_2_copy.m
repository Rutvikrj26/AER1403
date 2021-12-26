clear all; clc
%% Geometry & Constants
%Hole size
a = 0.2; b = 0.1;
%Cell Number
nx = 5;
ny = 5;
%Mesh Density Gradian
gradient = 1.1;
type = 1; %Mesh Type: 0-equally spaced, 1-high density at center

%Plate Properties
sigma = 200e6;
E = 210e9; %Young's Modulus [Pa]
v = 0.3; %Poisson Ratio
t = 0.001; %Thickness [mm]
D = (E/(1-v^2))*[1 v 0; v 1 0; 0 0 (1-v)/2]; %Matrix of Elastic Constants

%% Boundaries
%Ellipse shape
c_e = [a*cos(linspace(0,pi/2,ny+nx-1));b*sin(linspace(0,pi/2,ny+nx-1))];

%4 Boundaries
%c_1 = [zeros(1,ny);tan(linspace(0,pi/4,ny))];
c_2 = [tan(linspace(pi/4,0,nx));zeros(1,nx)+1];
c_3 = [zeros(1,ny)+1; tan(linspace(0,pi/4,ny))];
%c_4 = [tan(linspace(0,pi/4,nx)); zeros(1,nx)];

%% Generate Mesh
switch type
    case 0
        X=zeros(nx,ny);
        Y=zeros(nx,ny);
        for i = 1:ny
            X(i,:) = linspace(c_e(1,i),c_3(1,i),nx);
            Y(i,:) = linspace(c_e(2,i),c_3(2,i),nx);
        end
        for i = 1:nx
            X(i+ny-1,:) = linspace(c_e(1,i+ny-1),c_2(1,i),ny);
            Y(i+ny-1,:) = linspace(c_e(2,i+ny-1),c_2(2,i),ny);
        end
    case 1
    newX=zeros(nx,ny);
    newY=zeros(nx,ny);
    
% Generate Node Coordiantes
    R = nx+ny;
    for i = 1:ny
        newX(i,1) = c_e(1,i);
        newX(i,2) = newX(i,1)+(c_3(1,i)-c_e(1,i)) * ((1-gradient)/(1-gradient^(R-1)));
        for j = 3:R
            newX(i,j) = newX(i,j-1) + (newX(i,j-1)-newX(i,j-2))*gradient;
        end
        newX(i,R) = c_3(1,i);
    end
    for i = 1:nx
        newX(i+ny-1,1) = c_e(1,i+ny-1);
        newX(i+ny-1,2) = newX(i+ny-1,1)+(c_2(1,i)-c_e(1,i+ny-1)) * ((1-gradient)/(1-gradient^(R-1)));
        for j = 3:R
            newX(i+ny-1,j) = newX(i+ny-1,j-1) + (newX(i+ny-1,j-1)-newX(i+ny-1,j-2))*gradient;
        end
        newX(i+ny-1,R) = c_2(1,i);
    end
    for i = 1:nx
        newY(i,1) = c_e(2,i);
        newY(i,2) = newY(i,1)+(c_3(2,i)-c_e(2,i)) * ((1-gradient)/(1-gradient^(R-1)));
        for j = 3:R
            newY(i,j) = newY(i,j-1) + (newY(i,j-1)-newY(i,j-2))*gradient;
        end
        newY(i,R) = c_3(2,i);
    end
    for i = 1:ny
        newY(i+nx-1,1) = c_e(2,i+nx-1);
        newY(i+nx-1,2) = newY(i+nx-1,1)+(c_2(2,i)-c_e(2,i+nx-1)) * ((1-gradient)/(1-gradient^(R-1)));
        for j = 3:R
            newY(i+nx-1,j) = newY(i+nx-1,j-1) + (newY(i+nx-1,j-1)-newY(i+nx-1,j-2))*gradient;
        end
        newY(i+nx-1,R) = c_2(2,i);
    end
    X = newX;   Y = newY;
    clear newX newY;
end

%Ploting the Mesh
node_label = [1:length(X)]';

plot_1 = 0; 
if plot_1 == 1
    figure(1)
    hold on
    for i = 1:ny+nx-1
        plot(X(i,:),Y(i,:),'k')
        %text(X(i,:),Y(i,:),cellstr(num2str(node_label)),'VerticalAlignment','bottom','HorizontalAlignment','right')
        %node_label = node_label + 10;
    end
    plot(c_e(1,:), c_e(2,:),'k')
    plot(c_2(1,:), c_2(2,:),'k')
    plot(c_3(1,:), c_3(2,:),'k')
    plot(X,Y,'k')
    xlim([-0.2 1.2])
    ylim([-0.2 1.2])
    hold off
end

%% Nodes & Elements Numbering
figure(2)
hold on
for i = 1:ny+nx-1
    plot(X(i,:),Y(i,:),'k')
    text(X(i,:),Y(i,:),cellstr(num2str(node_label)),'VerticalAlignment','bottom','HorizontalAlignment','right')
    node_label = node_label + length(X);
end
plot(c_e(1,:), c_e(2,:),'k')
plot(c_2(1,:), c_2(2,:),'k')
plot(c_3(1,:), c_3(2,:),'k')
plot(X,Y,'k')
xlim([-0.2 1.2])
ylim([-0.2 1.2])
hold off

%Determine Node Number on top,btm,side
[M,N] = size(X);
LowerNode = 1:N;
UpperNode = linspace(ceil(M/2)*N,M*N,floor(N/2));
SideNode = LowerNode + (M-1)*N;

%Element Number
num_elemt = (M-1)*(N-1);
num_node = M*N;
order_elemt = zeros(num_elemt,4);
local_order = [1,2,2+N,1+N];

%Determine the node order and coordinates in each element.
n = 1;
for i = 1:M-1
    for j = 1:N-1
        order_elemt(n,:) = local_order;
        elemt_coordi_X(n,:) = [X(i,j),X(i,j+1),X(i+1,j+1),X(i+1,j)];
        elemt_coordi_Y(n,:) = [Y(i,j),Y(i,j+1),Y(i+1,j+1),Y(i+1,j)];
        local_order = local_order+1;
        n = n+1;
    end
    local_order = local_order+1;
end

%% Global Stiffness Matrix
Kg = zeros(2*num_node,2*num_node);
for i = 1:num_elemt
    L = gatherMat(num_node,order_elemt,i);
    Kg = Kg + L'*stiffnessMat(D,elemt_coordi_X(i,:),elemt_coordi_Y(i,:))*L*t;
end

%% Force Vecor
F = zeros(2*num_node,1);
% Unknown Vertical Forces(Y) Along Lower Nodes
for i = 1:length(LowerNode)
    F(2*LowerNode(i)) = nan;
end
% Known Vertical Forces Along Upper Nodes
for i = 1:length(UpperNode)-1
    F_node = (1/2)*t*abs(UpperNode(i)-UpperNode(i+1))*sigma;
    F(2*UpperNode(i),:) = F(2*UpperNode(i),:) + F_node;
    F(2*UpperNode(i+1),:) = F(2*UpperNode(i+1),:) + F_node;
end
% Unknown Horizontal Forces(X) Along Side Nodes
for i = 1:length(LowerNode)
    F(2*SideNode(i)-1) = nan;
end

%% Displacement Vector
D = nan*zeros(2*num_node,1);
for i = 1:length(D)
    if isnan(F(i))
        D(i) = 0;
    end
end

%% Solving
F_c = F;    K_c = Kg;   P = [];

% Partition the system from bottom up
for i = -length(D):-1
     if D(-i) == 0
        P = [P,-i];
        %Remove data from matrix
        F_c(-i) = [];
        K_c(-i,:) = []; K_c(:,-i) = [];
     end
end

D_c = K_c\F_c;
n = 0;
for i = 1:length(D)
    if isnan(D(i))
        n = n+1;
        D(i) = D_c(n);
    end
end

clear K_c D_c F_c;

F = Kg*D;
dx = D(1:2:end);    dy = D(2:2:end);

%%
n=1;
for i = 1:M
    for j =1:N
        new_X(i,j) = X(i,j)+dx(n);
        new_Y(i,j) = Y(i,j)+dy(n);
        n = n+1;
    end
end

figure(3)
hold on
plot(new_X,new_Y,'b')
plot(X,Y,'k')

plot_3 = 0; 
if plot_3 == 1
    figure(1)
    hold on
    for i = 1:ny+nx-1
        plot(X(i,:),Y(i,:),'k')
        plot(new_X(i,:),new_Y(i,:),'b')
    end
%     plot(c_e(1,:), c_e(2,:),'k')
%     plot(c_2(1,:), c_2(2,:),'k')
%     plot(c_3(1,:), c_3(2,:),'k')
    plot(X,Y,'k')
    plot(new_X,new_Y,'b')
    xlim([-0.2 1.2])
    ylim([-0.2 1.2])
    hold off
end
%% Deformed Mesh
% Undeformed Mesh Visualization
% scale = 60;
% visDeform = 0;
% if visDeform == 1
%     figure('Name','Deformed Mesh','Position', [100 100 600 600])
%     hold on
%     for i = 1:num_
%         plot(Nodes.X(Elements(i,[1:4,1])),Nodes.Y(Elements(i,[1:4,1])),'Color',[0 0 0 0.3])
%         plot((Nodes.X(Elements(i,[1:4,1]))+scale*Dx(Elements(i,[1:4,1]))'),(Nodes.Y(Elements(i,[1:4,1]))+scale*Dy(Elements(i,[1:4,1]))'),'r')
%     end
%     axis([-0.2 1.2 -0.2 1.2])
%     pbaspect([1,1,1])
%     daspect([1,1,1])
%     set(gca,'position',[0 0 1 1])
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(gca,'xcolor','w')
%     set(gca,'ycolor','w')
% end
%% Stiffness Matrix - Isoparametric Element
function K = stiffnessMat(D,X,Y)
C = [X',Y']; % Coordinate Matrix [m]
K = zeros(8,8);
for i = 1:2
    for j = 1:2
        eta = (2*i-3)/sqrt(3);
        xi = (2*j-3)/sqrt(3);
        J = (1/4)*[eta-1 1-eta 1+eta -eta-1; xi-1 -xi-1 1+xi 1-xi]*C;
        H = (1/4)*[eta-1 1-eta 1+eta -eta-1; xi-1 -xi-1 1+xi 1-xi];
        H = J\H;
        H = [H(1,1) 0 H(1,2) 0 H(1,3) 0 H(1,4) 0; 0 H(2,1) 0 H(2,2) 0 H(2,3) 0 H(2,4); H(2,1) H(1,1) H(2,2) H(1,2) H(2,3) H(1,3) H(2,4) H(1,4)];
        K = K + det(J)*H'*D*H;
    end
end
end

%% Gather Matrix Generation
function L = gatherMat(numNodes,Elements,N)
L = zeros(8,2*numNodes);
L(1:2,2*Elements(N,1)-1:2*Elements(N,1)) = eye(2);
L(3:4,2*Elements(N,2)-1:2*Elements(N,2)) = eye(2);
L(5:6,2*Elements(N,3)-1:2*Elements(N,3)) = eye(2);
L(7:8,2*Elements(N,4)-1:2*Elements(N,4)) = eye(2);
end

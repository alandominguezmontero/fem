clear;clc;clf;
%%                                       Definition of problem
n_el = 6;                                % number of elements
n_grid = n_el+1;                         % number of grids
%
A = 1000;                                % Rod Area [mm2]
E = 70000;                               % Young modulus [MPa]
L = 1000;                                % Length of the rod
P = -1e5 ;                               % Load value
%%
L = linspace(0,L,n_grid);                % x-coordinates of nodes
grids = [1:n_el; 
        2:n_el+1];                       % connectivities
A = A*ones(1,n_el);                      % area of each element
E = E*ones(1,n_el);                      % Young modulus of each element
%%                                       % BC, 1 = fixed, 0 = free	
bc = zeros(1,n_grid);
bc(1)   = 1;
bc(end) = 1;			                 	
%%                                       applied external forces on nodes
F=zeros(1,n_grid);
%F = P*[1/n_el*ones(1,n_grid)];          
F(floor((n_grid)/2)+1)=P;
%%                                       Stiffness matrix calculation
K = [zeros(n_grid,n_grid)];              
%
for i = 1:n_el                           
    l = L(grids(2,i)) - L(grids(1,i));   % compute element length
    k = E(i)*A(i)/l;                     % Rod stiffness (EA/L)
    k_el = [k, -k;
            -k, k];                         
                                         % global matrix assembly
    K(grids(1,i),grids(1,i)) = K(grids(1,i),grids(1,i)) + k_el(1,1);
    K(grids(1,i),grids(2,i)) = K(grids(1,i),grids(2,i)) + k_el(1,2);
    K(grids(2,i),grids(1,i)) = K(grids(2,i),grids(1,i)) + k_el(2,1);
    K(grids(2,i),grids(2,i)) = K(grids(2,i),grids(2,i)) + k_el(2,2);
end

% Modification of stiffness matrix with boundary conditions
for j = 1:n_grid
    if  bc(j) == 1
        K(j,j) = 1E+50;               
        F(j)   = 0;
    end
end
%
x = F/K;                                            % solution of [F]=[K][x] 
%%
subplot(221), plot(L,x,'-r*','LineWidth',2)
hold on
plot(linspace(0,L(end),n_grid),zeros(1,n_grid),'ko-')
title('Displacement Field');
ylabel('Displacement')
xlabel('x position')
grid on
% compute strain field

for i=1:n_el                                  
    l=L(grids(2,i))-L(grids(1,i));            % compute element length
    elong=x(grids(2,i))-x(grids(1,i));        % elongation of each element
    strain(2*i-1)=elong/l;                  
    xx(2*i-1)=L(grids(1,i));
    strain(2*i)=strain(2*i-1);
    xx(2*i)=L(grids(2,i));
end

subplot(222), plot(xx,strain,'-g','LineWidth',2) 
hold on
plot(linspace(0,L(end),n_grid),zeros(1,n_grid),'ko-')
title('Strain Field');
ylabel('Strain')
xlabel('x position')
grid on

% compute stress tensor

for i=1:n_el                                
    l=L(grids(2,i))-L(grids(1,i));           % compute element length
    elong=x(grids(2,i))-x(grids(1,i));       % elongation of each element
    stress(2*i-1)=E(i)*elong/l;              % to plot, calculate the "stress" at each node
    stress(2*i)=stress(2*i-1);               % stress in each element is uniform
    xx(2*i-1)=L(grids(1,i));
    xx(2*i)=L(grids(2,i));
end

subplot(223)
plot(xx,stress,'-b','LineWidth',2)  % plot the FEM stresses (constant stress element)
hold on
plot(linspace(0,L(end),n_grid),zeros(1,n_grid),'ko-')
title('Stress Field');
xlabel('x position')
ylabel('Stress')
grid on
% compute reaction forces
r = bc;

for k=1:n_grid                              
    if bc(k)==1
        if k==1
            r(k)=(x(k+1))*K(k,k+1);
        elseif k>1 & k<n_grid
            r(k)=(x(k+1)-x(k-1))*K(k,k+1);
        else k==n_grid
            r(k)=(x(k-1))*K(k-1,k);
        end           
    end
end 
subplot(224)
plot(L,r,'-k','LineWidth',2)  
hold on
plot(linspace(0,L(end),n_grid),zeros(1,n_grid),'ko-')
title('Reactions');
xlabel('x position')
ylabel('Forces')
grid on
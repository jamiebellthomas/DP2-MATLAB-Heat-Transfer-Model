% set problem constants
conductivity=100;
heat=5000; % heat applied at source
temp_amb=30; % temperature at external boundary
length_x=2.2;
length_y=2.2;
element_size=0.05; % will be rounded to ensure a whole number of elements
tol=1e-10; % to stop iteration

% set location of source - two values are [min, max]
% will be rounded to the nearest whole element
source_x=[1.0,1.2]; 
source_y=[1.0,1.2];

% calculate some values for the source
nx=round(length_x/element_size); % number of nodes must be an integer
ny=round(length_y/element_size);
sx=round(source_x/element_size); % nodes for source
sy=round(source_y/element_size);
sx=sx(1)+1:sx(2); % turn into ranges
sy=sy(1)+1:sy(2);
if ismember(1,sx) || ismember(nx,sx) || ismember(1,sy) || ismember(ny,sy)
    error('Source elements cannot be on the boundary.')
end

% fill matrices with initial values (ambient temperature)
T_new=zeros(nx,ny)+temp_amb;
T_old=T_new;

% set up constant for heat applied to each source node
c1=(heat/(length(sx)*length(sy)))/(4*conductivity);

% set up main iteration loop
errc=1; erro=errc; it=0;
while errc/erro>tol
    
    % calculate updated temperatures for internal nodes of problem matrix using iteration formula
    % because we never touch the outer nodes, we never need to explicitly enforce the boundary conditions!
    for ix=2:nx-1
        for iy=2:ny-1
            T_new(ix,iy)=0.25*T_old(ix+1,iy)+0.25*T_old(ix-1,iy)+0.25*T_old(ix,iy+1)+0.25*T_old(ix,iy-1);
        end
    end
    
    % apply heat load at the source
    T_new(sx,sy)=T_new(sx,sy)+c1;
    
    % calculate mean change across the matrix
    errc=mean2(abs(T_new-T_old));
    if it==1
        erro=errc;
    end
    
    % swap old and new values for next iteration and increment count
    T_old=T_new;
    it=it+1;
    % continue while loop (unless change ratio drops below tolerance)
end

% display result with contour plot
x=[1:size(T_new,2)]*element_size;
y=[1:size(T_new,1)]*element_size;
[C,j] = contourf(x,y,T_new,500);
set(j,'LineColor','none')
colorbar
axis('equal')

% calculate total heat transferred
top=sum(T_new(:,2)-T_new(:,1));
bottom=sum(T_new(:,ny-1)-T_new(:,ny));
left=sum(T_new(2,:)-T_new(1,:));
right=sum(T_new(nx-1,:)-T_new(nx,:));
heat_tot=(top+bottom+left+right)*conductivity
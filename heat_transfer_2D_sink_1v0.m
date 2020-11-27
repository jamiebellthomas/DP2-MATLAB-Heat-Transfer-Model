% set problem constants
conductivity=0.3;
convection_coeff=3;
temp_amb=30; % temperature at external boundary
temp_sink=5; % temperature at internal boundary (sink)
length_x=2.2;
length_y=2.2;
element_size=0.1; % will be rounded to ensure a whole number of elements
tol=1e-10; % to stop iteration

% set location of sink - two values are [min, max]
% will be rounded to the nearest whole element
source_x=[0.6,1.7];
source_y=[0.6,1.7];

% calculate some values for the sink
nx=round(length_x/element_size); % number of nodes must be an integer
ny=round(length_y/element_size);
sx=round(source_x/element_size); % nodes for sink
sy=round(source_y/element_size);
sx=sx(1)+1:sx(2); % turn into ranges
sy=sy(1)+1:sy(2);
%sx and sy are the unit steps of source's dimensions. So what this says is if the sink's
% boundaries hit any external boundries, it generates an error.
if ismember(1,sx) || ismember(nx,sx) || ismember(1,sy) || ismember(ny,sy)
    error('Sink elements cannot be on the boundary.')
end

% fill matrices with initial values (sink temperature)
T_new=zeros(nx,ny)+temp_sink;
T_old=T_new;

% set up constants for convection and conduction
%Based on the over all energy balance for slides and corners that are in
%the slides

%Nearly Edge Constant #1 
%(1/2*c2= Edge Constant #1) 
c2=(convection_coeff*element_size/conductivity)+2;
%Edge Constant / 2
c3=convection_coeff*element_size/conductivity;
%Corner Constant #1
c4=conductivity /(2*(conductivity+convection_coeff*element_size));
%Corner constant #2
c5=2*convection_coeff*element_size/conductivity;

% set up neighbour kernel for convolution
neighbours=[[0,.25,0];[.25,0,.25];[0,.25,0]];

% set up main iteration loop
errc=1; erro=errc; it=0;
while errc/erro>tol
    
    % calculate updated temperatures for internal nodes of problem matrix using convolution
    %Moves the kernal materix around the main matrix multiplying the
    %respect
    T_new(2:nx-1,2:ny-1)= conv2(T_old,neighbours,'valid');
    
    % calculate updated temperatures for edge nodes of problem matrix using convection condition
    T_new(2:nx-1,ny)=(1/(2*c2))*(2*T_old(2:nx-1,ny-1)+T_old(1:nx-2,ny)+T_old(3:nx,ny)+2*c3*temp_amb);
    T_new(2:nx-1,1)=(1/(2*c2))*(2*T_old(2:nx-1,2)+T_old(1:nx-2,1)+T_old(3:nx,1)+2*c3*temp_amb);
    T_new(nx,2:ny-1)=(1/(2*c2))*(2*T_old(nx-1,2:ny-1)+T_old(nx,1:ny-2)+T_old(nx,3:ny)+2*c3*temp_amb);
    T_new(1,2:ny-1)=(1/(2*c2))*(2*T_old(2,2:ny-1)+T_old(1,1:ny-2)+T_old(1,3:ny)+2*c3*temp_amb);
    
    % calculate updated temperatures for corner nodes of problem matrix using convection condition
    T_new(1,1)=c4*(T_old(2,1)+T_old(1,2)+c5*temp_amb);
    T_new(nx,ny)=c4*(T_old(nx-1,ny)+T_old(nx,ny-1)+c5*temp_amb);
    T_new(1,ny)=c4*(T_old(1,ny-1)+T_old(2,ny)+c5*temp_amb);
    T_new(nx,1)=c4*(T_old(nx-1,1)+T_old(nx,2)+c5*temp_amb);
    
    % enforce fixed temperature at the sink
    T_new(sx,sy)=temp_sink;
    
    % calculate mean change across the matrix
    T_sum=T_new-T_old;
    errc=mean2(abs(T_sum));
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
top=(T_new(sx,sy(1)-1)-T_new(sx,sy(1)));
bottom=(T_new(sx,sy(end)+1)-T_new(sx,sy(end)));
left=(T_new(sx(1)-1,sy)-T_new(sx(1),sy));
right=(T_new(sx(end)+1,sy)-T_new(sx(end),sy));
heat_tot=(top+bottom+left+right)*conductivity
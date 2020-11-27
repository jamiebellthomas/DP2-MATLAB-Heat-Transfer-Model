tic
% 3D Heat Sink
%Conduction Constant
k=0.031;
%Concevtion Coefficient
h=1.5;
%Ambient and sink temps
temp_amb=40;
temp_sink=0;
%Node dimensions
element_size=0.002;
%Factor that will determine when iterations stop
tol=1e-5;
%Dimensions 
length_x=0.565;
length_y=0.41;
length_z=0.138;
%Insulation Thickness
thickness=0.03;
%Sink Dimensions
source_x=[thickness, length_x-thickness];
source_y=[thickness, length_y-thickness];
source_z=[thickness, length_z-thickness];
%Number of nodes in each direction
nx=round(length_x/element_size); % number of nodes must be an integer
ny=round(length_y/element_size);
nz=round(length_z/element_size);

sx=round(source_x/element_size); % nodes for sink
sy=round(source_y/element_size);
sz=round(source_z/element_size);
sx=sx(1):sx(2); % turn into ranges
sy=sy(1):sy(2);
sz=sz(1):sz(2);
%Make sure sink does't cross external boundaries. 
if ismember(1,sx) || ismember(nx,sx) || ismember(1,sy) || ismember(ny,sy) || ismember(1,sz) || ismember(nz,sz)
    error('Sink elements cannot be on the boundary.')
end
%Starting matricies 
T_new=zeros(nx,ny,nz)+temp_sink;
T_old=T_new;
%Set up Kernal for conv2 function
kernal=[0,0,0;0,1/6,0;0,0,0];
kernal(:,:,2)=[0,1/6,0;1/6,0,1/6;0,1/6,0];
kernal(:,:,3)=[0,0,0;0,1/6,0;0,0,0];
%Define frequently used constants
c1=k*element_size;
c2=h*(element_size)^2;

c1_corner=(c1/(3*(c1+c2)));
c2_corner=((c2*temp_amb)/(c2+c1));
c1_edge=(c1/(4*c1+2*c2));
c2_edge=((c2*temp_amb)/(c2+2*c1));
c1_face=(c1/(5*c1+c2));
c2_face=((c2*temp_amb)/(c2+5*c1));

%Iteration counter
it=0;
mean_old=1;
mean_new=1;
error=1;

while error>tol
    %Internal Nodes
    T_new(2:nx-1,2:ny-1,2:nz-1)= convn(T_old,kernal,'valid');
    %Corners
    T_new(1,1,1)=c1_corner*(T_old(1,1,2)+T_old(1,2,1)+T_old(2,1,1))+c2_corner;
    T_new(nx,1,1)=c1_corner*(T_old(nx,1,2)+T_old(nx,2,1)+T_old(nx-1,1,1))+c2_corner;
    T_new(1,ny,1)=c1_corner*(T_old(1,ny,2)+T_old(1,ny-1,1)+T_old(2,ny,1))+c2_corner;
    T_new(1,1,nz)=c1_corner*(T_old(1,1,nz-1)+T_old(1,2,nz)+T_old(2,1,nz))+c2_corner;
    T_new(nx,ny,1)=c1_corner*(T_old(nx,ny,2)+T_old(nx,ny-1,1)+T_old(nx-1,ny,1))+c2_corner;
    T_new(nx,1,nz)=c1_corner*(T_old(nx,1,nz-1)+T_old(nx,2,nz)+T_old(nx-1,1,nz))+c2_corner;
    T_new(1,ny,nz)=c1_corner*(T_old(1,ny,nz-1)+T_old(1,ny-1,nz)+T_old(2,ny,nz))+c2_corner;
    T_new(nx,ny,nz)=c1_corner*(T_old(nx,ny,nz-1)+T_old(nx,ny-1,nz)+T_old(nx-1,ny,nz))+c2_corner;
    %Edges x axis
    T_new(2:nx-1,1,1)=c1_edge*(T_old(1:nx-2,1,1)+T_old(3:nx,1,1)+T_old(2:nx-1,2,1)+T_old(2:nx-1,1,2))+c2_edge;
    T_new(2:nx-1,ny,1)=c1_edge*(T_old(1:nx-2,ny,1)+T_old(3:nx,ny,1)+T_old(2:nx-1,ny-1,1)+T_old(2:nx-1,ny,2))+c2_edge;
    T_new(2:nx-1,1,nz)=c1_edge*(T_old(1:nx-2,1,nz)+T_old(3:nx,1,nz)+T_old(2:nx-1,2,nz)+T_old(2:nx-1,1,nz-1))+c2_edge;
    T_new(2:nx-1,ny,nz)=c1_edge*(T_old(1:nx-2,ny,nz)+T_old(3:nx,ny,nz)+T_old(2:nx-1,ny-1,nz)+T_old(2:nx-1,ny,nz-1))+c2_edge;
    %Edges y axis
    T_new(1,2:ny-1,1)=c1_edge*(T_old(1,1:ny-2,1)+T_old(1,3:ny,1)+T_old(2,2:ny-1,1)+T_old(1,2:ny-1,2))+c2_edge;
    T_new(nx,2:ny-1,1)=c1_edge*(T_old(nx,1:ny-2,1)+T_old(nx,3:ny,1)+T_old(nx-1,2:ny-1,1)+T_old(nx,2:ny-1,2))+c2_edge;
    T_new(1,2:ny-1,nz)=c1_edge*(T_old(1,1:ny-2,nz)+T_old(1,3:ny,nz)+T_old(2,2:ny-1,nz)+T_old(1,2:ny-1,nz-1))+c2_edge;
    T_new(nx,2:ny-1,nz)=c1_edge*(T_old(nx,1:ny-2,nz)+T_old(nx,3:ny,nz)+T_old(nx-1,2:ny-1,nz)+T_old(nx,2:ny-1,nz-1))+c2_edge;
    %Edges z axis
    T_new(1,1,2:nz-1)=c1_edge*(T_old(1,1,1:nz-2)+T_old(1,1,3:nz)+T_old(2,1,2:nz-1)+T_old(1,2,2:nz-1))+c2_edge;
    T_new(nx,1,2:nz-1)=c1_edge*(T_old(nx,1,1:nz-2)+T_old(nx,1,3:nz)+T_old(nx-1,1,2:nz-1)+T_old(nx,2,2:nz-1))+c2_edge;
    T_new(1,ny,2:nz-1)=c1_edge*(T_old(1,ny,1:nz-2)+T_old(1,ny,3:nz)+T_old(2,ny,2:nz-1)+T_old(1,ny-1,2:nz-1))+c2_edge;
    T_new(nx,ny,2:nz-1)=c1_edge*(T_old(nx,ny,1:nz-2)+T_old(nx,ny,3:nz)+T_old(nx-1,ny,2:nz-1)+T_old(nx,ny-1,2:nz-1))+c2_edge;    
    %Faces
    %xy plane faces
    T_new(2:nx-1,2:ny-1,1)=c1_face*(T_old(1:nx-2,2:ny-1,1)+T_old(3:nx,2:ny-1,1)+T_old(2:nx-1,1:ny-2,1)+T_old(2:nx-1,3:ny,1)+T_old(2:nx-1,2:ny-1,2))+c2_face; 
    T_new(2:nx-1,2:ny-1,nz)=c1_face*(T_old(1:nx-2,2:ny-1,nz)+T_old(3:nx,2:ny-1,nz)+T_old(2:nx-1,1:ny-2,nz)+T_old(2:nx-1,3:ny,nz)+T_old(2:nx-1,2:ny-1,nz-1))+c2_face; 
    %xz plane faces
    T_new(2:nx-1,1,2:nz-1)=c1_face*(T_old(1:nx-2,1,2:nz-1)+T_old(3:nx,1,2:nz-1)+T_old(2:nx-1,1,1:nz-2)+T_old(2:nx-1,1,3:nz)+T_old(2:nx-1,2,2:nz-1))+c2_face; 
    T_new(2:nx-1,ny,2:nz-1)=c1_face*(T_old(1:nx-2,ny,2:nz-1)+T_old(3:nx,ny,2:nz-1)+T_old(2:nx-1,ny,1:nz-2)+T_old(2:nx-1,ny,3:nz)+T_old(2:nx-1,ny-1,2:nz-1))+c2_face; 
    %yz plane faces
    T_new(1,2:ny-1,2:nz-1)=c1_face*(T_old(1,1:ny-2,2:nz-1)+T_old(1,3:ny,2:nz-1)+T_old(1,2:ny-1,1:nz-2)+T_old(1,2:ny-1,3:nz)+T_old(2,2:ny-1,2:nz-1))+c2_face; 
    T_new(nx,2:ny-1,2:nz-1)=c1_face*(T_old(nx,1:ny-2,2:nz-1)+T_old(nx,3:ny,2:nz-1)+T_old(nx,2:ny-1,1:nz-2)+T_old(nx,2:ny-1,3:nz)+T_old(nx-1,2:ny-1,2:nz-1))+c2_face; 
    
    %Reset Sink Temp
    T_new(sx,sy,sz)=temp_sink;

    %Set up errors:

    error=mean2(abs(T_new-T_old));
    
    T_old=T_new;
    
   
    

    it=it+1;
end

[A,B,C]=meshgrid(1:nx,1:ny,1:nz);
xslice = [ny];                            
yslice = [round(sx(round(length(sx)/2))) round((thickness/(2*length_x))*nx)];
zslice = [1 round(sz(round(length(sz)/2)))];
slice(T_new,xslice,yslice,zslice);
cb = colorbar;                
cb.Label.String = 'Temperature, C';
hold on

width=(sx(end)-sx(1)+1);
height=(sy(end)-sy(1)+1);
len=(sz(end)-sz(1)+1);
x1=[0 width width 0 0];
y1=[0 0 height height 0];
a1=[x1;x1];
b1=[y1;y1];
c1=[0 0 0 0 0;len len len len len];
a1=a1+sx(1);
b1=b1+sy(1);
c1=c1+sz(1);
mesh(b1,a1,c1,'FaceAlpha', 0.2);
axis equal;
hold on;

y2=[0 0 len len 0];
b2=[y2;y2];
c2=[0 0 0 0 0;height height height height height];


row_counter=1;
angle=(pi/2);
while row_counter<3;
    column_counter=1;
    while column_counter<6;
        b_new=(b2(row_counter,column_counter))*cos(angle)-(c2(row_counter,column_counter))*sin(angle);
        c_new=(b2(row_counter,column_counter))*sin(angle)+(c2(row_counter,column_counter))*cos(angle);
        b2(row_counter,column_counter)=b_new;
        c2(row_counter,column_counter)=c_new;
        column_counter=column_counter+1;
        disp(column_counter);
        
    end
    row_counter=row_counter+1;
end
b2=b2+height;
b2=b2+(sy(1));
c2=c2+(sy(1));
mesh(b2,a1,c2,'FaceAlpha', 0.2);

front=(T_new(sx,sy,sz(1))-T_new(sx,sy,(sz(1)-1)))*-1;
back=(T_new(sx,sy,sz(end))-T_new(sx,sy,sz(end)+1))*-1;
top=(T_new(sx,sy(1),sz)-T_new(sx,(sy(1)-1),sz))*-1;
bottom=(T_new(sx,sy(end),sz)-T_new(sx,(sy(end)+1),sz))*-1;
left=(T_new(sx(1),sy,sz)-T_new((sx(1)-1),sy,sz))*-1;
right=(T_new(sx(end),sy,sz)-T_new((sx(end)+1),sy,sz))*-1;







toc
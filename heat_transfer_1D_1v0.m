%
% set problem constants
%
conductivity=100;
length=4;
area=0.01;
%
% set boundary conditions
%
heat=200
%Heat applied in k
temp_bound=30;
%Heat at edges
spread=0.2;
%Spread is the width at which head is applied
location=1.6;
%Location where head is applied relative to LHS
%
% set initial values and solution parameters
%
nx=400;  % ( number of subdividing elements )
tol=1e-10;  %  ( tolerance parameter used to gauge when solution has converged )
it=0;  %  ( set up iteration count parameter )
errc=1;  %   ( set reference for initial loop )
erro=errc;
%
% fill matrices with initial values ( guesses )
%
T_new=zeros(nx)+temp_bound;
T_old=zeros(nx)+temp_bound;
%
% calculate element length
%
dx=length/nx;
%
% work out node numbers corresponding to position of heat load
%Numbeer of nodes in the loaded area.
no_nodes_load=(spread/length)*nx;
%heat per unit area that's applied
heat_unit_length=heat/spread;
% Tells where the centre of nodes are located
location_nodes=round((location/length)*nx);
%Lower limit of nodes
first_node_l=round(location_nodes-ceil((no_nodes_load/2)));
%Upper Limit nodes
last_node_l=round(location_nodes+ceil((no_nodes_load/2)));

%3 Lines above basically defines the node areas. 

%
% set up constant used in iteration calculation
%
c1=0.5*heat_unit_length*dx*dx/(conductivity*area);
   %{
   This constant is in the notes, it forms part of the equation of working
   out the temp of the current node 
   %}
%
% set up main iteration loop
%
while errc/erro>tol
   %
   % calculate updated temperatures using iteration formula
   %
   %This range is defined going from 2nd to penultimate (1st and last = 30)
   for ix=2:nx-1
    T_new(ix)=0.5*T_old(ix+1)+0.5*T_old(ix-1);

   %
   % apply heat load term on applicable loads
   %
        if (ix>first_node_l) && (ix<last_node_l)
        T_new(ix)=T_new(ix)+c1;
        end
   end
   %
   % calculate error values
   %
   errc=0;
   for ix=1:nx   
    errc=errc+abs(T_new(ix)-T_old(ix));   %  ( calculate error as sum of changes between iterations ) 
   end
    errc=errc/nx; 
   if it==1
    erro=errc;   %  ( set error reference to error value found in 1st iteration )
   end
   % swap old and new values ready for next iteration
   for ix=1:nx
    T_old(ix)=T_new(ix);
   end
   %
   % enforce boundary conditions
   %
   T_old(1)=temp_bound;
   T_old(nx)=temp_bound;
   %
   % increment iteration count
   %
   it=it+1;
   %
   % continue while loop (unless error ratio exceeds tolerance)
   %
end
%
% display result
%
plot(0:dx:length-dx,T_new);
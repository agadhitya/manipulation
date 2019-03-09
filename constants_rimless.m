
localDir = fileparts(mfilename('fullpath')) ;      
restoredefaultpath ;% clear paths before adding
addpath(fullfile(localDir, 'symbolic_functions')) ;

m      = [1 1]   ;  %mass matrix
l1      = 1     ;   %length of first link
l2      = 1     ;   % length of second link

grav  = 9.8;
    
I1 = ones(3);        %Inertia of the first link

I2 = ones(3);        %Interia of the second link

L_com=[0.5 1.5];     %COM length from origin

DH=[ 0 l1 0 0;...
     0 l2 0 0];      % DH parameters

n = length( DH(:,1) ) ;     % number of links


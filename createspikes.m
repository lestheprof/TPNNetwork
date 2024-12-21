% generate some spike input data
% format is <time neuron>
N_ap = 100 ;
N_ba = 120 ;

r1 = rand([1 N_ap]);
ap1(:,2) =  round(r1) + 1 ;
r2 = rand([1 N_ap]);
ap1(:,1) =  r2 ;
r1 = rand([1 N_ba]);
ba1(:,2) =  round(5 * r1) + 1 ;
r2 = rand([1 N_ba]);
ba1(:,1) =  r2 ;
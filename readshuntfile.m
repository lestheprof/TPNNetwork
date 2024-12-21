function shunts = readshuntfile(fname, simulation)
%readneuronfile reads file fname into a table and converts to a structure
%
% LSS 11 Dec 2024
% modified 15 Dec 2024, added noapicalshunts and nobasalshunts
% 
[n_lines, shunttable] = readnetwork(fname) ;
%
if simulation.N_TPNs > n_lines
    error("readshuntfile: simulation.N_TPNs > n_lines") ;
end
for tpnno = simulation.N_TPNs:-1:1 % for each TPN
    shunts(tpnno).apicalshuntduration = shunttable(tpnno, 1).apicalshuntduration;
    shunts(tpnno).basalshuntduration = shunttable(tpnno, 2).basalshuntduration;
    shunts(tpnno).noapicalshunts = shunttable(tpnno, 3).noapicalshunts ;
    shunts(tpnno).nobasalshunts = shunttable(tpnno, 4).nobasalshunts ;
    shunts(tpnno).apicalshuntweights = zeros([1, shunts(tpnno).noapicalshunts]) ; 
    shunts(tpnno).basalshuntweights =  zeros([1, shunts(tpnno).nobasalshunts]) ;
end
% Other members of shunts structure added later.
end

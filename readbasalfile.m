function basal = readbasalfile(fname, simulation)
%readbasalfile reads file fname into a table and converts to a structure
%
% LSS 11 Dec 2024
% updated LSS 15 Dec 2024, added n_basal_inputs
% 
[n_lines, basaltable] = readnetwork(fname) ;
%
if simulation.N_TPNs > n_lines
    error("readbasalfile: simulation.N_TPNs > n_lines") ;
end
for basno = simulation.N_TPNs:-1:1 % for each TPN
    basal(basno).n_basalinputs = basaltable(basno, 1).n_basalinputs;
    basal(basno).tau_basal = basaltable(basno, 2).tau_basal;
    basal(basno).C_basal = basaltable(basno, 3).C_basal ;
    basal(basno).R_basal = basaltable(basno, 4).R_basal ;
    basal(basno).R_synba_dendrite = basaltable(basno, 5).R_synba_dendrite ;
    basal(basno).C_basal_spine = basaltable(basno, 6).C_basal_spine ;
    basal(basno).R_basal_spine = basaltable(basno, 7).R_basal_spine ;
    basal(basno).R_synba_spine = basaltable(basno, 8).R_synba_spine ;
    basal(basno).synapsemultiplier = basaltable(basno, 9).synapsemultiplier ;
    basal(basno).basalsynapseweights = zeros([1 basal(basno).n_basalinputs]) ;
end
% Other members of basal structure added later.
end

function apical = readapicalfile(fname, simulation)
%readapicalfile reads file fname into a table and converts to a structure
%
% LSS 11 Dec 2024 started
% las modified 15 Dec2024, added n_apical_inputs
% 
[n_lines, apicaltable] = readnetwork(fname) ;
%
if simulation.N_TPNs > n_lines
    error("readapicalfile: simulation.N_TPNs > n_lines") ;
end
for apno = simulation.N_TPNs:-1:1 % for each TPN
    apical(apno).n_apicalinputs = apicaltable(apno, 1).n_apicalinputs;
    apical(apno).tau_apical = apicaltable(apno, 2).tau_apical;
    apical(apno).C_apical = apicaltable(apno, 3).C_apical ;
    apical(apno).R_apical = apicaltable(apno, 4).R_apical ;
    apical(apno).R_synap_dendrite = apicaltable(apno, 5).R_synap_dendrite ;
    apical(apno).C_apical_spine = apicaltable(apno, 6).C_apical_spine ;
    apical(apno).R_apical_spine = apicaltable(apno, 7).R_apical_spine ;
    apical(apno).R_synap_spine = apicaltable(apno, 8).R_synap_spine ;
    apical(apno).synapsemultiplier = apicaltable(apno, 9).synapsemultiplier ;
    apical(apno).apicalsynapseweights = zeros([1, apical(apno).n_apicalinputs]) ;
end
% Other members of apical structure added later.
end

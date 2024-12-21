function neuron = readneuronfile(fname, simulation)
%readneuronfile reads file fname into a table and converts to a structure
%
% LSS 11 Dec 2024
% 
[n_lines, neurontable] = readnetwork(fname) ;
%
if simulation.N_TPNs > n_lines
    error("readneuronfile: simulation.N_TPNs > n_lines") ;
end
for tpnno = simulation.N_TPNs:-1:1 % for each TPN
    neuron(tpnno).refractoryperiod = neurontable(tpnno, 1).refractoryperiod;
    neuron(tpnno).relrefperiod = neurontable(tpnno, 2).relrefperiod;
    neuron(tpnno).thresh_leap = neurontable(tpnno, 3).thresh_leap ;
    neuron(tpnno).thresh_decay = neurontable(tpnno, 4).thresh_decay ;
    neuron(tpnno).thresh_value = neurontable(tpnno, 5).thresh_value ;
    neuron(tpnno).maxnospikes = neurontable(tpnno, 6).maxnospikes ;
end
% Other members of neuron structure added later.
end

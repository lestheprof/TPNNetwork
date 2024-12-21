function IIneuron = readIIneuronfile(fname, simulation)
%readIIneuronfile reads IIneuronfile from file fname
%
% started LSS 15 Dec 2024
%
[n_lines, IItable] = readnetwork(fname) ;
%
if simulation.N_IIs > n_lines
    error("readIIneuronfile: simulation.N_IIs > n_lines") ;
end
for IIno = simulation.N_IIs:-1:1
    IIneuron(IIno).n_synapses = IItable(IIno,1).n_synapses ;
    IIneuron(IIno).tau = IItable(IIno,2).tau ;
    IIneuron(IIno).C = IItable(IIno,3).C ;
    IIneuron(IIno).R = IItable(IIno,4).R ;
    IIneuron(IIno).refractoryperiod = IItable(IIno,5).refractoryperiod ;
    IIneuron(IIno).relrefperiod = IItable(IIno,6).relrefperiod ;
    IIneuron(IIno).thresh_value = IItable(IIno,7).thresh_value ;
    IIneuron(IIno).thresh_leap = IItable(IIno,8).thresh_leap ;
    IIneuron(IIno).thresh_decay = IItable(IIno,9).thresh_decay ;
    IIneuron(IIno).maxnospikes = IItable(IIno,10).maxnospikes ;
    IIneuron(IIno).synapsemultiplier = IItable(IIno,11).synapsemultiplier ;
    IIneuron(IIno).weights = zeros([1 IIneuron(IIno).n_synapses]) ;
    IIneuron(IIno).spikes = zeros([1 IIneuron(IIno).maxnospikes]) ;
    
end
end
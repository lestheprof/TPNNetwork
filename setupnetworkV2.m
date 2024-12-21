function [simulation, neuron, basal,apical, shunts, apicalcurrent, basalcurrent, ...
    apicalactivation, basalactivation, ahactiv, IIneuron] = setupnetworkV2(simulation, neuron, basal,apical, shunts,IIneuron)
%% setupnetwork Sets up the network for simulation,
% V2 is used in the v4 of the simulation, and uses setupAlphaFunctionV2


%% (in structure simulation)
% duration: time in seconds for running the neuron
% timestep:  timestep length in seconds
% N_TPNs: number of TPN neurons
% N_IIs: number of inhibitory interneurons

%% (in structure neuron)
% parameters for spike generation
% thresh_value: initial threshold for dtetermining spike (applied to ahactiv
% value)
% refractoryperiod: refractory period during which no spikes will be
% generated.
% relrefperiod: relative refractory period after which threshold returns to
%  initial value
% thresh_leap : leap in threshold value after Refractory period
% thresh_decay : exponent of decay (exp(-thresh_decay * time))
% parameters for returning spike list
% neuronid: identity of this neuron (integer)
% maxnospikes: maximum number of spikes this neuron can return

%% (in structure basal (array of structures))
% nobasalinputs: number of basal inputs
% basalinputs: same format
% basalsynapsewights: weights on basal synapses, must be same length as nobasalinputs
% tau-basal: time constant for basal synapses (alpha function)
% concept is that dendrite has a local capacitance and membrane
% resistance. In addition there is a resistance describing the maximal
% conductance (=minimum resistance) of the dendrite that governs the amount
% of charge transferred. Keep one set of values for apical and one for
% basal.
% r_synba_dendrite: resistance of basal dendrite:
% c_basal: capacitance at basal dendrite
% r_basal: basal leakage resistance
% added for spines: currently same for all spines, but this is likely to
% be what gets modified for adaptation
% c_basal_spine: capacitance at basal spine
% r_basal_spine: basal lspineeakage resistance
% r_synba_spine: resistance of basal spine: weight
% (basalsynapseweights) adjusts this per synapse

%% (in structure apical (array of structures))
% noapicalinputs:  number of apical inputs
% apicalinputs: actual apical spike inputs in <time synapse_number> format
% apicalsynapseweights: weights on apical synapses, must be same length as noapicalinputs
% tau_apical: time constant for apical synapses (alpha function)
% concept is that dendrite has a local capacitance and membrane
% resistance. In addition there is a resistance describing the maximal
% conductance (=minimum resistance) of the dendrite that governs the amount
% of charge transferred. Keep one set of values for apical and one for
% basal.
% c_apical: capacitance at apical dendrite
% r_apical: apical leakage resistance
% r_synap_dendrite:
% added for spines: currently same for all spines, but this is likely to
% be what gets modified for adaptation
% c_apical_spine: capacitance at apical spinendrite: resistance (1/conductance) of apical dendrite.
% r_apical spine:
% r_synap_spine: resistance (1/conductance) of apical spine. weight
% (apicalsynapseweights) adjusts this per synapse% r_apical_spine: apical spine leakage resistance

%% (in structure shunts (array of strcutures))
% noapicalshunts: number of apical shunting synapses
% nobasalshunts: number of basal shunting synapses
% apicalshuntinputs: actual apical shunt inputs, in <time synapse_number> format
% basalshuntinputs: actual basal shunt inputs, in <time synapse_number> format
% apicalshuntweights: weights for apical shunt: 0 leq value leq 1
% basalshuntweights: weights for basal shunt: 0 leq value leq 1
% apicalshuntduration: time that apical shunt is active for (in seconds)
% basalshuntduration: time that basal shunt is active for (in seconds)

%% in structure IIneuron
% refractoryperiod:
% relrefperiod:
% thresh_leap
% thresh_decay
% thresh_value
% maxnospikes
% alpha_synapse: array from alpha function
% activation: array to hold activation
% II_th_inc: array to hold threshold increment
% II_threshold: array to hold actual threshold


% length in timesteps of simulation
simulation.simlength = ceil(simulation.duration/simulation.timestep) ;

% calculate alpha functions for apical and basal dendrites (tau need not be
% the same: we might later want a variety of these for excitatory and inhibitory synapses)
maxlengthaplha_apical = 1 ; % used later to ensure arrays are long enough
maxlengthalpha_basal = 1 ;
for tpnno = 1:simulation.N_TPNs
    apical(tpnno).alpha_apical = setupAlphaFunctionV2(simulation.timestep, apical(tpnno).tau_apical) ;
    if length(apical(tpnno).alpha_apical) > maxlengthaplha_apical
        maxlengthaplha_apical = length(apical(tpnno).alpha_apical) ;
    end
    basal(tpnno).alpha_basal = setupAlphaFunctionV2(simulation.timestep, basal(tpnno).tau_basal );
    if length(basal(tpnno).alpha_basal) > maxlengthalpha_basal
        maxlengthalpha_basal = length(basal(tpnno).alpha_basal) ;
    end
end

%% array declarations, 1st index is neuron no.
%% would these be better attached to the relevant structure?
% Specifically apicalpostsynapse and basalpostsynapse, as the others are 1 per neuron.
% array for apical synapse activations (extra length to stop late spikes
% falling off the end)
% the per-synapse post-synaptic values are used in this version
for tpnno = 1:simulation.N_TPNs
    apical(tpnno).apicalpostsynapse = zeros([apical(tpnno).n_apicalinputs (simulation.simlength  + 10 * maxlengthaplha_apical)]);
    % array for basal activations
    basal(tpnno).basalpostsynapse = zeros([basal(tpnno).n_basalinputs simulation.simlength + 10 * maxlengthalpha_basal]);
    % vector for apical current
    % vector for threshold: useful for seeing what has been going on
    % neuron(tpnno).TP_th_inc = calc_thresh_increment(neuron(tpnno).thresh_leap, neuron(tpnno).thresh_decay, ...
    %     neuron(tpnno).refractoryperiod, neuron(tpnno).relrefperiod, simulation.timestep);
    % neuron(tpnno).TP_threshold = ones([1, simulation.simlength + length(neuron(tpnno).TP_th_inc)]) * neuron(tpnno).thresh_value ;basalshuntweights(1, 1:2)
    neuron(tpnno).TP_threshold = ones([1, simulation.simlength + length(neuron(tpnno).thresh_increment)]) * neuron(tpnno).thresh_value ;

    % calculate apical leakiness per timestep:how much leaks away from
    % apicalactivation voltage in 1 timestep
    % C * dV/dt = I = V/R => dV/dt = V/(R*C)
    % so fraction that leaks away is timestep/(R*C)
    apical(tpnno).apicalfracleak = simulation.timestep / (apical(tpnno).R_apical * apical(tpnno).C_apical) ;
    basal(tpnno).basalfracleak = simulation.timestep / (basal(tpnno).R_basal * basal(tpnno).C_basal) ;
    % values below not currently used in TPNrunstep_v3, but should be in v4
    apical(tpnno).apicalspinefracleak = simulation.timestep/(apical(tpnno).R_apical_spine  * apical(tpnno).C_apical_spine) ;
    % apical(tpnno).apicalspinefracleak = 0 ;
    basal(tpnno).basalspinefracleak = simulation.timestep/(basal(tpnno).R_basal_spine  * basal(tpnno).C_basal_spine) ;
    % basal(tpnno).basalspinefracleak =  0 ;
end
apicalcurrent = zeros([simulation.N_TPNs simulation.simlength + 10 * maxlengthaplha_apical]) ;
% vector for basal current
basalcurrent = zeros([simulation.N_TPNs simulation.simlength + 10 * maxlengthalpha_basal]) ;
% vector for apical activation (charged by apicalcurrent, capacitor and
% parallel resistor)
apicalactivation = zeros([simulation.N_TPNs simulation.simlength + 10 * maxlengthaplha_apical]) ;
% vector for basal activation (charged by basalcurrent, capacitor and
% parallel resistor)f
basalactivation = zeros([simulation.N_TPNs simulation.simlength + 10 * maxlengthalpha_basal]) ;
% vector for axon hillock activation
ahactiv = zeros([simulation.N_TPNs simulation.simlength]) ;




% calc-thresh_increment has been calculated

% sort apical and basal inputs, additive and shunting, into time order
% both apical and basal inputs are a 2d array with
% N_spikes rows 2 columns <time neuron>
for tpnno = 1:simulation.N_TPNs

    if  (isfield(apical, 'apicalinputs') && ~isempty(apical(tpnno).apicalinputs) )% allow empty apical input
        apical(tpnno).apicalinputs = sortrows(apical(tpnno).apicalinputs, 1) ;
        % replace times with timestep numbers
        apical(tpnno).apicalinputs(:,1) = round(apical(tpnno).apicalinputs(:,1)/simulation.timestep);
        % check that the apical spike inputs are within range
        if (max(apical(tpnno).apicalinputs(:,2)) > apical(tpnno).n_apicalinputs)
            error("setupnetwork: apical input for TPN %d out of range", tpnno);
        end
    end
    if  (isfield(basal, 'basalinputs') && ~isempty(basal(tpnno).basalinputs)) % allow empty basal input
        basal(tpnno).basalinputs = sortrows(basal(tpnno).basalinputs, 1) ;
        % replace times with timestep numbers
        basal(tpnno).basalinputs(:,1) = round(basal(tpnno).basalinputs(:,1)/simulation.timestep);
        % check that the basal spike inputs are within range
        if (max(basal(tpnno).basalinputs(:,2)) > basal(tpnno).n_basalinputs)
            error("setupnetwork: basal input for TPN %d out of range", tpnno);
        end
    end
    if (isfield(shunts,'apicalshuntinputs') && ~isempty(shunts(tpnno).apicalshuntinputs))
        shunts(tpnno).apicalshuntinputs = sortrows(shunts(tpnno).apicalshuntinputs, 1);
        % replace times with timestep numbers
        shunts(tpnno).apicalshuntinputs(:,1) = round(shunts(tpnno).apicalshuntinputs(:,1)/simulation.timestep);
        % check that the apical spike inputs are within range
        if (max(shunts(tpnno).apicalshuntinputs(:,2)) > shunts(tpnno).noapicalshunts)
            error("setupnetwork: apical shunt input for TPN %d out of range", tpnno);
        end
    end
    if (isfield(shunts, 'basalshuntinputs') && ~isempty(shunts(tpnno).basalshuntinputs))
        shunts(tpnno).basalshuntinputs = sortrows(shunts(tpnno).basalshuntinputs, 1);
        % replace times with timestep numbers
        shunts(tpnno).basalshuntinputs(:,1) = round(shunts(tpnno).basalshuntinputs(:,1)/simulation.timestep);
        % check that the basal spike inputs are within range
        if (max(shunts(tpnno).basalshuntinputs(:,2)) > shunts(tpnno).nobasalshunts)
            error("setupnetwork: basal shunt input for TPN %d out of range", tpnno);
        end
    end
    % set up shunting synapses but only if there are shunting synapses
    if (shunts(tpnno).noapicalshunts > 0)
        shunts(tpnno).asduration = floor(shunts(tpnno).apicalshuntduration/simulation.timestep) ; % get apical shunt duration in timesteps
        if (shunts(tpnno).asduration < 1)
            shunts(tpnno).asduration = 1 ;
        end % 0 is possible, but needs reset to 1
        % calculate shunt value per timestep for each shunting synapse weight
        shunts(tpnno).ap_sh_wt_pt = zeros([1 length(shunts(tpnno).apicalshuntweights)]) ; % preallocate
        for sno = 1: length(shunts(tpnno).apicalshuntweights)
            shunts(tpnno).ap_sh_wt_pt(sno) = (1-shunts(tpnno).apicalshuntweights(sno))^(1/shunts(tpnno).asduration); % weight to apply per timestep
        end
    end
    if (shunts(tpnno).nobasalshunts  > 0)
        shunts(tpnno).bsduration = floor(shunts(tpnno).basalshuntduration/simulation.timestep) ; % get basal shunt duration in timesteps
        if (shunts(tpnno).bsduration < 1)
            shunts(tpnno).bsduration = 1 ;
        end % 0 is possible, but needs reset to 1
        shunts(tpnno).ba_sh_wt_pt = zeros([1 length(shunts(tpnno).basalshuntweights)]) ; % preallocate
        for sno = 1: length(shunts(tpnno).basalshuntweights)
            shunts(tpnno).ba_sh_wt_pt(sno) = (1 - shunts(tpnno).basalshuntweights(sno))^(1/shunts(tpnno).bsduration); % weight to apply per timestep
        end
    end

end
%  sort external II inpiuts, if any, and change times to timestep values
for IIno = 1:simulation.N_IIs
    if  (isfield(IIneuron, 'inputs') && ~isempty(IIneuron(IIno).inputs) )% allow empty II input
        IIneuron(IIno).inputs = sortrows(IIneuron(IIno).inputs, 1) ;
        % replace times with timestep numbers
        IIneuron(IIno).inputs(:,1) = round(IIneuron(IIno).inputs(:,1)/simulation.timestep);
        % check that the apical spike inputs are within range
        % if (max(IIneuron(IIno).inputs(:,2)) > IIneuron(IIno).n_inputs)
        %     error("setupnetwork: apical input for TPN %d out of range", IIno);
    end
end



%% setup II neuron information
for IIno = 1:simulation.N_IIs
    % set up activation
    % first calculate the max length of the alpha function
    IIneuron(IIno).alpha_synapse = setupAlphaFunctionV2(simulation.timestep, IIneuron(IIno).tau );
    IIneuron(IIno).activation = zeros ([1 simulation.simlength + 10 * length(IIneuron(IIno).alpha_synapse)]) ;
    % set up threhold array
    IIneuron(IIno).II_th_inc = calc_thresh_increment(IIneuron(IIno).thresh_leap, IIneuron(IIno).thresh_decay, ...
        IIneuron(IIno).refractoryperiod, IIneuron(IIno).relrefperiod, simulation.timestep);
    IIneuron(IIno).II_threshold = ones([1, simulation.simlength + length(IIneuron(IIno).II_th_inc)]) * IIneuron(IIno).thresh_value ;
    IIneuron(IIno).fracleak = simulation.timestep/(IIneuron(IIno).R * IIneuron(IIno).C) ;
end


end
function [basal, apical, shunts, IIsynapse] = setupweights(weightfile, basal, apical, shunts, IIsynapse)
% setupweights% sets up the weights reading from the table in weightfile


% read in table of weights from weightfile
[n_weights, weighttable] = readnetwork(weightfile) ;

for wtno = 1:n_weights % test each arc
    switch char(weighttable(wtno,1).neuron_type)
        case 'TPN' % two point neuron weights
            switch char(weighttable(wtno,3).syn_type)
                case 'A'
                    apical(weighttable(wtno,2).neuron_number).apicalsynapseweights(weighttable(wtno,4).syn_number) = weighttable(wtno,5).weight ;
                case 'B'
                    basal(weighttable(wtno,2).neuron_number).basalsynapseweights(weighttable(wtno,4).syn_number) = weighttable(wtno,5).weight ;
                case 'AS'
                    shunts(weighttable(wtno,2).neuron_number).apicalshuntweights(weighttable(wtno,4).syn_number) = weighttable(wtno,5).weight ;
                case 'BS'
                    shunts(weighttable(wtno,2).neuron_number).basalshuntweights(weighttable(wtno,4).syn_number) = weighttable(wtno,5).weight ;
                otherwise
                    error("setupweights: TPN synapse type %s is not valid", weighttable(wtno,3).syn_type) ;
            end
        case 'II' % inhibitory interneuron
            switch char(weighttable(wtno,3).syn_type)
                case 'S'
                    IIsynapse(weighttable(wtno,2).neuron_number).weights(weighttable(wtno,4).syn_number) = weighttable(wtno,5).weight ;
                otherwise
                    error("setupweights: II synapse type %s is not valid", weighttable(wtno,3).syn_type) ;
            end

        otherwise
            error("setupweights: Neuron type %s is not valid", weighttable(wtno,1).neuron_type)
    end

end
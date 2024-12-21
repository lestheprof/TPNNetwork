function [basalinputs, apicalinputs, apicalshuntinputs, basalshuntinputs, IIinputs] = readexternalinput(fname)
%readexternalinput reads external spiking input from file fname
%   B basal input, A apical input, AS apical shunt, BS basal shunt, II
%   inhibitory ionternmeuron
% LSS 17 Dec 2024 started
% 
% read file into table
[n_lines, inputtable] = readnetwork(fname) ;
% preallocate space for all the different inputs
basalinputs = zeros([n_lines 3]) ;
apicalinputs = zeros([n_lines 3]) ;
apicalshuntinputs = zeros([n_lines 3]) ;
basalshuntinputs = zeros([n_lines 3]) ;
IIinputs = zeros([n_lines 3]) ;
bno = 0;
ano = 0 ;
asno = 0 ;
bsno = 0 ;
iino = 0 ;
for inputno = 1:n_lines % for each input
    switch char(inputtable(inputno,1).inputtype)
        case 'A'
            ano = ano + 1;
            apicalinputs(ano, 1) = inputtable(inputno,2).Nno ;
            apicalinputs(ano, 2) = inputtable(inputno,3).time ;
            apicalinputs(ano, 3) = inputtable(inputno,4).synapseno ;
        case 'B'
            bno = bno + 1 ;
            basalinputs(bno, 1) = inputtable(inputno,2).Nno ;
            basalinputs(bno, 2) = inputtable(inputno,3).time ;
            basalinputs(bno, 3) = inputtable(inputno,4).synapseno ;
        case 'AS'
            asno = asno + 1 ;
            apicalshuntinputs(asno, 1) = inputtable(inputno,2).Nno ;
            apicalshuntinputs(asno, 2) = inputtable(inputno,3).time ;
            apicalshuntinputs(asno, 3) = inputtable(inputno,4).synapseno ;
        case 'BS'
            bsno = bsno + 1 ;
            basalshuntinputs(bsno, 1) = inputtable(inputno,2).Nno ;
            basalshuntinputs(bsno, 2) = inputtable(inputno,3).time ;
            basalshuntinputs(bsno, 3) = inputtable(inputno,4).synapseno ;
        case 'II'
            iino = iino+ 1 ;
            IIinputs(iino, 1) = inputtable(inputno,2).Nno ;
            IIinputs(iino, 2) = inputtable(inputno,3).time ;
            IIinputs(iino, 3) = inputtable(inputno,4).synapseno ;
        otherwise
            error("readexternalinput: invalid neuron type in external input") ;
    end
    % shorten arrays of input prior to returning them
    apicalinputs = apicalinputs(1:ano, :) ;
    basalinputs = basalinputs(1:bno, :) ;
    apicalshuntinputs = apicalshuntinputs(1:asno,:) ;
    basalshuntinputs = basalshuntinputs(1:bsno,:) ;
    IIinputs = IIinputs(1:iino,:) ;
end
end
function [ signal ] = cumulated_phase_rcomplement( sekvence )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sekvence=upper(sekvence); 
sekvence=seqrcomplement(sekvence);
a=1+1i;c=-1-1i;g=-1+1i;t=1-1i;
zapis=zeros(1,length(sekvence),'double');
zapis(sekvence=='A')=angle(a);
zapis(sekvence=='C')=angle(c);
zapis(sekvence=='G')=angle(g);
zapis(sekvence=='T')=angle(t);
signal=zapis; 

end


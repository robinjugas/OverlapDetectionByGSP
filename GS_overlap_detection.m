function [order] = GS_overlap_detection(fasta,threshold)
%% Readme

N=length(fasta);
for i=1:N %detectiong the length of shortest read
    contig_length(i)=length(fasta{i});
    shortest_length=min(contig_length);
    [pom,MAX]=max(contig_length);
end

seqs=upper(fasta);
% threshold=0.60; %covariance coefficient threshold
step=5; %length of covariance computing windows
steps_number=int32((shortest_length/4)/step); %choosen as 1/4 of read length

a=1+1i;c=-1-1i;g=-1+1i;t=1-1i; %numerical representation formula
for i=1:N
    read_length(i)=length(seqs{1,i});
    zapis=zeros(1,read_length(i),'double');
    zapis(seqs{i}=='A')=angle(a);
    zapis(seqs{i}=='C')=angle(c);
    zapis(seqs{i}=='G')=angle(g);
    zapis(seqs{i}=='T')=angle(t);
    SequencesNumS{i,1}=cumsum(zapis); %variable of numerical reads coming into overlap detection
%     SequencesNumS{i,1}=SequencesNumS{i,1}-mean(SequencesNumS{i,1});
end

%% Overlap detection
pool=[1:N]; %vector of remaining reads
order=[]; %vector of joined reads
unused=[]; %vector of disjoined reads
for j=1:N-1 %going through all contigs except the last
    if j==1 %start with the longest contig
        ACTUAL_contig=SequencesNumS{MAX,1}; %actual contig which is used for detection with others
        pool(pool==MAX)=[]; %delete actual read from pool
        order=[MAX]; %put the read in the result
    end  
    
    LAST=[]; %choosen read to be joined    
    covariance_matrix_x=zeros(N,steps_number); %covariance matrix N x steps_number for !!! 
    covariance_matrix_y=zeros(N,steps_number); %covariance matrix N x steps_number for !!!
    
    %poèítání korelaèní matice
    for p=1:length(pool) %computing covariance matrix step-by-step
        P=pool(p); %actual choosen read
        NEXT_contig=SequencesNumS{P,1}; %next contig, contig actually in computing
        
        %X actual sufix meeting another prefix
        for kx=1:steps_number
            if kx==1
                step=10;
                sig1x=ACTUAL_contig(end-kx*step+1:end);
                sig2x=NEXT_contig(1:kx*step);
            else
                step=5;
                sig1x=ACTUAL_contig(end-kx*step+1-5:end);
                sig2x=NEXT_contig(1:(kx*step)+5);
            end
            covariance_matrix_x(P,kx)=xcov(sig1x,sig2x,0,'coeff');
            
        end
        
        %Y actual prefix meeting another sufix
        for ky=1:steps_number
            if ky==1
                step=10;
                sig1y=ACTUAL_contig(1:ky*step);
                sig2y=NEXT_contig(end-ky*step+1:end);
            else
                step=5;
                sig1y=ACTUAL_contig(1:(ky*step)+5);
                sig2y=NEXT_contig(end-ky*step+1-5:end);
            end           
            covariance_matrix_y(P,ky)=xcov(sig1y,sig2y,0,'coeff');
        end
    end
    
    %Choosing max covariance coefficient over contigs
    matrixX=abs(covariance_matrix_x); %transform into absolute number
    matrixY=abs(covariance_matrix_y);
    [maxX]=max(matrixX(:)); %find maximum value in both matrices
    [maxY]=max(matrixY(:));
    [rx,cx] = find(matrixX == maxX,1); %find maximum value read in matrices rx-read,cx-position
    [ry,cy] = find(matrixY == maxY,1);
    
    
    if (maxX  > threshold) || (maxY > threshold) %covariance should be higher than threshold
        
        if maxX > maxY %sufix joined with some prefix
            LAST=rx; %write new last contig
            pool(pool==LAST)=[]; %delete from pool
            order=[order LAST]; %write new solution
            
            if covariance_matrix_x(rx,cx) < 0 %the LAST read is reverse complementary
                %transform into original sequence and join two reads at the maximum position
                pom_contig=[ACTUAL_contig(1:length(ACTUAL_contig)-(cx*step)-5) GS_cumulated_phase_rcomplement(seqs{LAST})];
                ACTUAL_contig=pom_contig; 
            else %is not reverse
                pom_contig=[ACTUAL_contig(1:length(ACTUAL_contig)-(cx*step)-5) SequencesNumS{LAST,1}];
                ACTUAL_contig=pom_contig;              
            end
            
        elseif maxX < maxY %prefix joined with some sufix
            LAST=ry; %write new last contig
            pool(pool==LAST)=[]; %delete from pool
            order=[LAST order]; %write new solution
            
            if covariance_matrix_y(ry,cy) < 0 %the LAST read is reverse complementary
                pom_contig=[GS_cumulated_phase_rcomplement(seqs{LAST}) ACTUAL_contig(1+5+cx*step:end)];
                ACTUAL_contig=pom_contig;
            else
                pom_contig=[SequencesNumS{LAST,1} ACTUAL_contig(1+5+cx*step:end)];
                ACTUAL_contig=pom_contig;
            end
        end
        
    else % no read joined, low maximum covariance, choose another random start read
%         if j==1
%             X=randi(length(pool));
%             pool(pool==X)=[];
%             ACTUAL_contig=SequencesNumS{X,1};
%             order=[order 0 X];
%         else
%             unused=pool;
%             return
%         end
        break
    end
end

end

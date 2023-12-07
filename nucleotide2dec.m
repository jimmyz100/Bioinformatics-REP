function decimal =  nucleotide2dec(Nucleotide_Seq)
% Function to convert nucleotide sequence to a base-10 number. The sequence
% is first converted into a base-4 number based on the following:
% G = 0; A = 1; T = 2; C = 3
% The base-4 number is then converted into a base-10 number
Base_4_String = num2str([]);
s = size(Nucleotide_Seq);
for i = 1:s(2) % Loop for conversion of nucleotide base to base-4 digit
    if Nucleotide_Seq(i) == 'G' 
        Base_4(i) = 0;
    elseif Nucleotide_Seq(i) == 'A'
        Base_4(i) = 1;
    elseif Nucleotide_Seq(i) == 'T'
        Base_4(i) = 2;
    elseif Nucleotide_Seq(i) == 'C'   
        Base_4(i) = 3;     
    else
        disp('Invalid base in file.')
    end
    b = num2str(Base_4(i));
    Base_4_String = strcat(Base_4_String,b);
end
decimal = base2dec(Base_4_String,4);
end
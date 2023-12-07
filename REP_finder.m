function REP_finder(Genome_Start)
% Genome Start: Number of the genome to start on based on the order of the
% list

% REP Expressions used in regexp (and inputs)
Min_Middle_Section = 0;
Max_Middle_Section = 32;
General_REP_Expression = ['(G[CT]C[CGT]GA).{' num2str(Min_Middle_Section) ',' num2str(Max_Middle_Section) '}(??@palindrome_flip($1))'];

GetChromoAccNumfromAssemSearch % Retrieves accession numbers

% Remove duplicates in the compiled accession number list from previous function and sorts in ascending order: 
for c = 1:length(List)
   list{c} = List{c,2};
end
Accession_List = sort(unique(list));

REP_Palindrome_List = []; % Creates list of all palindromes for development of key (later on)

% Make directories for organization of files
if 7 ~= exist('Genome Repeats data')
    mkdir('Genome Repeats data');
end
if 7 ~= exist('Genome Repeats data/REP-like')
    mkdir('Genome Repeats data/REP-like');
end

% Read existing file if a previous run is done
if exist('Genome Repeats data/REP-like/REP_One_Line_All_Genomes.csv','file') == 2
rep_one_line_table = readtable('Genome Repeats data/REP-like/REP_One_Line_All_Genomes.csv');
REP_One_Line = table2struct(rep_one_line_table);
end

switch nargin % Allows code to run regardless of presence of input
    case 1 % If genome_start was input
        if exist('Genome Repeats data/REP-like/REP_Counter_Palindrome_List.csv','file') == 2
        REP_Palindrome_List = table2array(readtable('Genome Repeats data/REP-like/REP_Counter_Palindrome_List.csv'));
        end
        if Genome_Start == 1 % Erases data from previous runs from recording table if beginning new run
            REP_One_Line = [];
        end
for a = Genome_Start:length(Accession_List) % Loops code for all accession numbers under "List"
    tic
warning('off', 'MATLAB:mode:EmptyInput') %suppresses unimportant warning 
downloadGenomeKEH(Accession_List{a}); %downloads the genome
fid=fileread(['NCBI data/Sequence/' Accession_List{a} '.sq']);
fid1=fopen(['NCBI data/Sequence/' Accession_List{a} '.sq']);
seq = textscan(fid1, '%s %*[^\n]','HeaderLines',1);
sequence=char(seq{1}); %saves the sequence without the header for later use
fclose('all');

% Name of Genome and Genome Length
Genome_name_Expression = [Accession_List{a} '(.\d\s(\w+\s\w+)|\s(\w+\s\w+)|\s(\[\w+\s\w+\])|\s(\[\w+\]\s\w+)|\s(\''\w+\s\w+\'')|\s(\w+\-\w+\s\w+)'];
Genome_name = regexp(fid,Genome_name_Expression,'tokens');
Genome_length = length(sequence);

% Clear REP data structure to prevent overlaps in future runs
REP = [];

% To test for general REP Sequences (General Expression)
[startIndex3, endIndex3] = regexp(sequence,General_REP_Expression);
[matches3,tokens3] = regexp(sequence,General_REP_Expression,'match','tokens');

% Determine distances between hits
if isempty(startIndex3) == 0
for m = 1:length(startIndex3)-1
    Distance(m) = startIndex3(m+1)-endIndex3(m);
end
REP_Seq_Mean_Distance = round(mean(Distance),2);
else
    REP_Seq_Mean_Distance = 0;
end

% Putting unique general hits into REP structure
if isempty(startIndex3) == 0
for i = 1:length(startIndex3)
    REP(i).Accession_Number = Accession_List{a};
    if isempty(Genome_name) == 0
    REP(i).Genome_Name = Genome_name{1};
    else
    REP(i).Genome_Name = 'Not Found on Program';
    end
    REP(i).StartIndex=startIndex3(i);
    REP(i).EndIndex=endIndex3(i);
    REP(i).Total_Seq=matches3(i);
    REP(i).Total_Length=endIndex3(i)-startIndex3(i)+1;
    REP(i).Palindrome_Seq=cell2mat(tokens3{i});
    REP(i).Palindrome_Description = 'General REP expression based on REP palindromes and 1-30 bp in between';
    REP(i).Palindrome_Length = length(REP(i).Palindrome_Seq);
end
end

% Producing REP Excel File
if isempty(startIndex3) == 0 % if there are general matches
    REP = nestedSortStruct(REP,'Palindrome_Length',1); % function to sort the palindrome length from smallest to biggest
    REP(1).Reg_Expression = General_REP_Expression;
    REP(1).Genome_Length = Genome_length;
    REP(1).Total_Hits = length(startIndex3);
    REP(1).Mean_Distance = REP_Seq_Mean_Distance;
    REP_table = struct2table(REP);    
% Determine average length of REP Hits
    REP_Total_Lengths = REP_table(:,6);
    REP_Total_Mean_Length = round(mean(table2array(REP_Total_Lengths)),2);
    REP(length(startIndex3)+1).Total_Length = REP_Total_Mean_Length;
else % if there are no matches, provides table displaying "No REP-like matches Found"
    REP(1).Seq_Found = 'No REP-like matches Found';
    REP(1).Reg_Expression = 'REP-like Expression based on whole consensus or based on consensus REP paindromes and 1-30 bp in between';
    REP(2).Reg_Expression = General_REP_Expression;
    REP(1).Genome_Length = Genome_length;
    REP_table = struct2table(REP);
    REP_Total_Mean_Length = 0;
end

% Calculations for one-line output
Expected_REP_number = round(Genome_length/(4^6)); % Expected random number of genomes by probability

% Determine total number of REP hits
Total_Hits = length(startIndex3);

% Determine number of repetitions of REP palindromic sequences
if isempty(startIndex3) == 0
    if length(startIndex3) ~= 1 % if there are multiple hits
REP_Palindrome_Seq = table2array(REP_table(:,7));
[REP_Unique_Palindrome_Seq, Useless_Index, Repetition_Index] = unique(REP_Palindrome_Seq);
Number_of_Repetitions = histcounts(Repetition_Index);
Repetitions_Seq = REP_Unique_Palindrome_Seq;
Repetitions_Number = Number_of_Repetitions';
REP_Palindrome_List = [REP_Palindrome_List Repetitions_Seq'];
    else % if there is only one hit
       Repetitions_Seq = table2array(REP_table(:,7));
       Repetitions_Number = 1;
       REP_Palindrome_List = [REP_Palindrome_List Repetitions_Seq];
    end
else
    Repetitions_Seq = 'None';
end

% Produce a one-line output for each genome
REP_One_Line(a).Accession_Number = Accession_List{a};
REP_One_Line(a).Genome_Name = Genome_name{1};
REP_One_Line(a).Genome_Length = Genome_length;
REP_One_Line(a).REP_Hits = Total_Hits;
REP_One_Line(a).Expected_Random_REP_Hits = Expected_REP_number;
REP_One_Line(a).REP_Mean_Length = REP_Total_Mean_Length;
REP_One_Line(a).REP_General_Mean_Distance = REP_Seq_Mean_Distance;
if isempty(startIndex3) == 0
    if length(startIndex3) == 1
        REP_One_Line(a).CodeFreq1 = [num2str(nucleotide2dec(Repetitions_Seq)) '-' num2str(Repetitions_Number)];
    else
    for f = 1:length(Repetitions_Seq)
        CodeFreq = ['CodeFreq' num2str(f)];
        REP_One_Line(a).(CodeFreq) = [num2str(nucleotide2dec(Repetitions_Seq{f})) '-' num2str(Repetitions_Number(f))];
    end
    end
else
        REP_One_Line(a).CodeFreq1 = Repetitions_Seq;
end
REP_table = struct2table(REP);
writetable(REP_table, ['Genome Repeats data/REP-like/' Accession_List{a} 'REP_matches.xlsx'])
writetable(REP_table, ['Genome Repeats data/REP-like/' Accession_List{a} 'REP_matches.csv'])
REP_counter = a; % Counter file to determine how many genomes processed
REP_Counter = array2table(REP_counter);
writetable(REP_Counter,'REP_Counter.csv')
REP_Counter_Palindrome_List = array2table(unique(REP_Palindrome_List));
if isempty(REP_Counter_Palindrome_List) == 0
writetable(REP_Counter_Palindrome_List,'Genome Repeats data/REP-like/REP_Counter_Palindrome_List.xlsx')
writetable(REP_Counter_Palindrome_List,'Genome Repeats data/REP-like/REP_Counter_Palindrome_List.csv')
end
REP_One_Line_Table = struct2table(REP_One_Line);
writetable(REP_One_Line_Table,'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.xlsx','Sheet','One-Line Data')
writetable(REP_One_Line_Table,'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.csv')
toc
end

% Determine number of genomes with hits vs no hits
Genome_List_Number = table2array(REP_One_Line_Table(:,1));

No_Hits_Number = 0;
No_Hits_Genomes = [];
Total_Genome = table2array(REP_One_Line_Table(:,4));
for m = 1:length(Total_Genome)
    if Total_Genome(m) == 0
        No_Hits_Number = No_Hits_Number + 1; % No_Hits_Number = Genomes with no hits
        No_Hits_Genomes = [No_Hits_Genomes m]; % List of rows where genome has no hits
    end
end
Hits_Number = length(Genome_List_Number) - No_Hits_Number; % Hits_Number = Genomes with hits
REP_One_Line(length(Genome_List_Number)+1).REP_Hits = Hits_Number;
REP_One_Line(length(Genome_List_Number)+2).REP_Hits = No_Hits_Number;

% Elliminate genomes without hits for subsequent mean calculations
REP_One_Line_Table_Calc = REP_One_Line_Table;
REP_One_Line_Table_Calc(No_Hits_Genomes,:) = [];

% Determine average REP_Mean_Length of the list of genomes with hits
Column_REP_Mean_Length = table2array(REP_One_Line_Table_Calc(:,6));
Average_REP_Mean_Length = mean(Column_REP_Mean_Length);
REP_One_Line(length(Genome_List_Number)+1).REP_Mean_Length = Average_REP_Mean_Length;

% Determine average REP_General_Mean_Distance of the list of genomes with
% hits
Column_REP_Mean_Distance = table2array(REP_One_Line_Table_Calc(:,7));
Average_REP_Mean_Distance = mean(Column_REP_Mean_Distance);
REP_One_Line(length(Genome_List_Number)+1).REP_General_Mean_Distance = Average_REP_Mean_Distance;
REP_One_Line_Table = struct2table(REP_One_Line);
writetable(REP_One_Line_Table,'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.xlsx','Sheet','One-Line Data')
writetable(REP_One_Line_Table,'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.csv')
case 0 % No input Genome_Start, so code starts from start of genome list
REP_One_Line = []; % Erases data from pervious runs to avoid confounding data
for a = 1:length(Accession_List) % Loops code for all accession numbers under "List"
    tic
warning('off', 'MATLAB:mode:EmptyInput') %suppresses unimportant warning 
downloadGenomeKEH(Accession_List{a}); %downloads the genome
fid=fileread(['NCBI data/Sequence/' Accession_List{a} '.sq']);
fid1=fopen(['NCBI data/Sequence/' Accession_List{a} '.sq']);
seq = textscan(fid1, '%s %*[^\n]','HeaderLines',1);
sequence=char(seq{1}); %saves the sequence without the header for later use
fclose('all');

% Name of Genome and Genome Length
Genome_name_Expression = [Accession_List{a} '(.\d\s(\w+\s\w+)|\s(\w+\s\w+)|\s(\[\w+\s\w+\])|\s(\[\w+\]\s\w+)|\s(\''\w+\s\w+\'')|\s(\w+\-\w+\s\w+)'];
Genome_name = regexp(fid,Genome_name_Expression,'tokens');
Genome_length = length(sequence);

% Clear REP data structure to prevent overlaps in future runs
REP = [];

% To test for general REP Sequences (General Expression)
[startIndex3, endIndex3] = regexp(sequence,General_REP_Expression);
[matches3,tokens3] = regexp(sequence,General_REP_Expression,'match','tokens');

% Determine distances between hits (general)
if isempty(startIndex3) == 0
for m = 1:length(startIndex3)-1
    Distance(m) = startIndex3(m+1)-endIndex3(m);
end
REP_Seq_Mean_Distance = round(mean(Distance),2);
else
    REP_Seq_Mean_Distance = 0;
end

% Putting unique general hits into REP structure
if isempty(startIndex3) == 0
for i = 1:length(startIndex3)
    REP(i).Accession_Number = Accession_List{a};
    if isempty(Genome_name) == 0
    REP(i).Genome_Name = Genome_name{1};
    else
    REP(i).Genome_Name = 'Not Found on Program';
    end
    REP(i).StartIndex=startIndex3(i);
    REP(i).EndIndex=endIndex3(i);
    REP(i).Total_Seq=matches3(i);
    REP(i).Total_Length=endIndex3(i)-startIndex3(i)+1;
    REP(i).Palindrome_Seq=cell2mat(tokens3{i});
    REP(i).Palindrome_Description = 'General REP expression based on REP palindromes and 1-30 bp in between';
    REP(i).Palindrome_Length = length(REP(i).Palindrome_Seq);
end
end

% Producing REP Excel File
if isempty(startIndex3) == 0 % if there are general matches
    REP = nestedSortStruct(REP,'Palindrome_Length',1); % function to sort the palindrome length from smallest to biggest
    REP(1).Reg_Expression = General_REP_Expression;
    REP(1).Genome_Length = Genome_length;
    REP(1).Total_Hits = length(startIndex3);
    REP(1).Mean_Distance = REP_Seq_Mean_Distance;
    REP_table = struct2table(REP);
    
% Determine average length of REP Hits
    REP_Total_Lengths = REP_table(:,6);
    REP_Total_Mean_Length = round(mean(table2array(REP_Total_Lengths)),2);
    REP(length(startIndex3)+1).Total_Length = REP_Total_Mean_Length;
else % if there are no matches, provides table displaying "No REP Sequences Found"
    REP(1).Seq_Found = 'No REP-like matches Found';
    REP(1).Reg_Expression = 'REP-like Expression based on whole consensus or based on consensus REP paindromes and 1-30 bp in between';
    REP(2).Reg_Expression = General_REP_Expression;
    REP(1).Genome_Length = Genome_length;
    REP_table = struct2table(REP);
    REP_Total_Mean_Length = 0;
end

% Calculations for one-line output
Expected_REP_number = round(Genome_length/(4^6)); % Expected random number of genomes by probability

% Determine total number of REP hits
Total_Hits = length(startIndex3);

% Determine number of repetitions of REP palindromic sequences
if isempty(startIndex3) == 0
    if length(startIndex3) ~= 1 % if there are multiple hits
REP_Palindrome_Seq = table2array(REP_table(:,7));
[REP_Unique_Palindrome_Seq, Useless_Index, Repetition_Index] = unique(REP_Palindrome_Seq);
Number_of_Repetitions = histcounts(Repetition_Index);
Repetitions_Seq = REP_Unique_Palindrome_Seq;
Repetitions_Number = Number_of_Repetitions';
REP_Palindrome_List = [REP_Palindrome_List Repetitions_Seq'];
    else % if there is only one hit
       Repetitions_Seq = table2array(REP_table(:,7));
       Repetitions_Number = 1;
       REP_Palindrome_List = [REP_Palindrome_List Repetitions_Seq];
    end
else
    Repetitions_Seq = 'None';
end

% Produce a one-line output for each genome
REP_One_Line(a).Accession_Number = Accession_List{a};
REP_One_Line(a).Genome_Name = Genome_name{1};
REP_One_Line(a).Genome_Length = Genome_length;
REP_One_Line(a).REP_Hits = Total_Hits;
REP_One_Line(a).Expected_Random_REP_Hits = Expected_REP_number;
REP_One_Line(a).REP_Mean_Length = REP_Total_Mean_Length;
REP_One_Line(a).REP_General_Mean_Distance = REP_Seq_Mean_Distance;
if isempty(startIndex3) == 0
    if length(startIndex3) == 1
        REP_One_Line(a).CodeFreq1 = [num2str(nucleotide2dec(Repetitions_Seq)) '-' num2str(Repetitions_Number)];
    else
    for f = 1:length(Repetitions_Seq)
        CodeFreq = ['CodeFreq' num2str(f)];
    REP_One_Line(a).(CodeFreq) = [num2str(nucleotide2dec(Repetitions_Seq{f})) '-' num2str(Repetitions_Number(f))];
    end
    end
else
        REP_One_Line(a).CodeFreq1 = Repetitions_Seq;
end
REP_table = struct2table(REP);
writetable(REP_table, ['Genome Repeats data/REP-like/' Accession_List{a} 'REP_matches.xlsx'])
writetable(REP_table, ['Genome Repeats data/REP-like/' Accession_List{a} 'REP_matches.csv'])
REP_counter = a;
REP_Counter = array2table(REP_counter);
writetable(REP_Counter,'REP_Counter.xlsx')
writetable(REP_Counter,'REP_Counter.csv')
REP_Counter_Palindrome_List = array2table(unique(REP_Palindrome_List));
if isempty(REP_Counter_Palindrome_List) == 0
writetable(REP_Counter_Palindrome_List,'Genome Repeats data/REP-like/REP_Counter_Palindrome_List.xlsx')
writetable(REP_Counter_Palindrome_List,'Genome Repeats data/REP-like/REP_Counter_Palindrome_List.csv')
end
REP_One_Line_Table = struct2table(REP_One_Line);
writetable(REP_One_Line_Table,'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.xlsx','Sheet','One-Line Data')
writetable(REP_One_Line_Table,'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.csv')
toc
end

% Determine number of genomes with hits vs no hits
No_Hits_Number = 0;
No_Hits_Genomes = [];
Total_Genome = table2array(REP_One_Line_Table(:,4));
for m = 1:length(Total_Genome)
    if Total_Genome(m) == 0
        No_Hits_Number = No_Hits_Number + 1; % No_Hits_Number = Genomes with no hits
        No_Hits_Genomes = [No_Hits_Genomes m]; % List of rows where genome has no hits
    end
end
Hits_Number = length(Total_Genome) - No_Hits_Number; % Hits_Number = Genomes with hits
REP_One_Line(length(Total_Genome)+1).REP_Hits = Hits_Number;
REP_One_Line(length(Total_Genome)+2).REP_Hits = No_Hits_Number;
% Clear extraneous numbers
REP_One_Line(length(Total_Genome)+1).REP_Hits = [];
REP_One_Line(length(Total_Genome)+2).REP_Hits = [];

% Elliminate genomes without hits for subsequent mean calculations
REP_One_Line_Table_Calc = REP_One_Line_Table;
REP_One_Line_Table_Calc(No_Hits_Genomes,:) = [];

% Determine average REP_Mean_Length of the list of genomes with hits
Column_REP_Mean_Length = table2array(REP_One_Line_Table_Calc(:,6));
Average_REP_Mean_Length = mean(Column_REP_Mean_Length);
REP_One_Line(length(Total_Genome)+1).REP_Mean_Length = Average_REP_Mean_Length;

% Determine average REP_General_Mean_Distance of the list of genomes with
% hits
Column_REP_Mean_Distance = table2array(REP_One_Line_Table_Calc(:,7));
Average_REP_Mean_Distance = mean(Column_REP_Mean_Distance);
REP_One_Line(length(Total_Genome)+1).REP_General_Mean_Distance = Average_REP_Mean_Distance;
REP_One_Line_Table = struct2table(REP_One_Line);
writetable(REP_One_Line_Table,'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.xlsx','Sheet','One-Line Data')
writetable(REP_One_Line_Table,'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.csv')
end

% Creates unique list of REP Palindromes for Code Key
REP_Unique_Palindrome_List = unique(REP_Palindrome_List);

% Write documentation with REP One Line Table
REP_One_Line_Documentation.REP_like_Expression_Used = General_REP_Expression;
REP_One_Line_Documentation.REP_Expression_Description = ['Expression of a REP-like sequence with defined 6-bp palindrome sections based on REP consensus sequences with ' num2str(Min_Middle_Section) '-' num2str(Max_Middle_Section) ' bp in between the palindrome repeats (Defined in Middle_Section parameters)'];
REP_One_Line_Documentation.Min_Middle_Section_bp = num2str(Min_Middle_Section);
REP_One_Line_Documentation.Max_Middle_Section_bp = num2str(Max_Middle_Section);
REP_One_Line_Documentation.Step_1 = 'Search NCBI assembly database for bacterial reference genomes and returns accession numbers using function GetChromoAccNumfromAssemSearch.m';
REP_One_Line_Documentation.Step_2 = 'Downloads genome sequence from NCBI database based on accession numbers using downloadGenomeKEH.m';
REP_One_Line_Documentation.Step_3 = 'Using MATLAB regexp function to locate indexes in genome with a sequence with 6-base REP-like (based on REP consensus) palindromes on both ends - expression found in "REP like Expression Used"';
REP_One_Line_Documentation.Step_4 = 'Descriptions of Parameters in Table explained:';
REP_One_Line_Documentation.Genome_Name = 'Name of Species from NCBI based on accession number';
REP_One_Line_Documentation.Genome_Length = 'Number of base pairs in whole genome';
REP_One_Line_Documentation.REP_Hits = 'Number of sequence matches based on the regexp function. Bottom two rows represent the number of genomes with REP-like hits, followed by number of genomes without REP-like hits';
REP_One_Line_Documentation.Expected_Random_REP_Hits = 'Expected Number of matches based on genome length divided by 4^6 (possible combinations of palindrome)';
REP_One_Line_Documentation.REP_Mean_Length = 'Average number of bp for all matches, Bottom Row displays the average of the mean lengths for the genome list';
REP_One_Line_Documentation.REP_General_Mean_Distance = 'Average distance in bp between each match, Bottom Row displays the average of the mean distances for the genome list';
REP_One_Line_Documentation.CodeFreq = 'Code-Frequency, Frequency of each 6-bp palindrome section sequence, Key for codes = palindrome sequences in next column. The key is based on converting the 6-bp sequence first into a 6 digit base-4 number, where G=0,A=1,T=2,C=3. Then, the base-4 number is converted into a base-10 number indicated by the "Code" here.';
for g = 1:length(REP_Unique_Palindrome_List)
REP_CodeKey = ['REP_CodeKey' num2str(g)];
REP_One_Line_Documentation.(REP_CodeKey) = [num2str(nucleotide2dec(REP_Unique_Palindrome_List{g})) '=' REP_Unique_Palindrome_List{g}];
end
REP_One_Line_Doc = struct2table(REP_One_Line_Documentation);
writetable(REP_One_Line_Doc, 'Genome Repeats data/REP-like/REP_One_Line_All_Genomes.xlsx','Sheet','Documentation')
writetable(REP_One_Line_Doc, 'Genome Repeats data/REP-like/REP_One_Line_All_Genomes_Doc.csv')
end
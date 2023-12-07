function Gene_Location(Genome_Start) 
% Genome_Start must be input. If start at the beginning, input 1 for
% Genome_Start.
warning('off', 'MATLAB:mode:EmptyInput') 

GetChromoAccNumfromAssemSearch % Retrieves accession numbers (in List)

% Remove duplicates in the compiled accession number list and sorts in ascending order: 
for c = 1:length(List)
   list{c} = List{c,2};
end
Accession_List = sort(unique(list));

% Makes directories to organize files
if 7 ~= exist('Genome Repeats data')
    mkdir('Genome Repeats data');
end
if 7 ~= exist('Genome Repeats data/General')
    mkdir('Genome Repeats data/General');
end

% Read existing file if a previous run is done
if exist('Genome Repeats data/General/General_One_Line_All_Genomes.xlsx','file') == 2
general_one_line_table = readtable('Genome Repeats data/General/General_One_Line_All_Genomes.xlsx');
General_One_Line = table2struct(general_one_line_table);
end

for a = Genome_Start:length(Accession_List)
    tic
downloadGenomeKEH(Accession_List{a}); 
[F,G] = analyze_feature_table(Accession_List{a});
for i = 1:length(G) 
    if G(i).End < G(i).Start % If end index is smaller than start index, inverts indices
        UG(i).Start = G(i).End;
        UG(i).End = G(i).Start;
    else
        UG(i).Start = G(i).Start;
        UG(i).End = G(i).End;
    end
    UG(i).Name = G(i).Name;
end

[part_sequence,startIndex,endIndex,General_table,Genome_name,General_Palindrome,Palindrome_Number,Min_Middle_Section,Max_Middle_Section] = general_repeat(Accession_List{a});
R = readtable(['Genome Repeats data/General/' Accession_List{a} 'gen_palindrome.xlsx']);
UGM = gene_locator(table2array(R(:,3)),table2array(R(:,4)),UG,Accession_List{a});
non_gene_entries = [];
for i = 1:length(UGM) % Elliminate non gene region sequences in gene-associated hits table
    if isempty(UGM(i).Name) == 1
        non_gene_entries = [non_gene_entries i];
    end
end
UGMref = UGM;
UGM_non_gene = UGMref(non_gene_entries);
UGM_non_gene = rmfield(UGM_non_gene,{'Name','Match_Attribute'});
UGMref(non_gene_entries) = [];
UGMref = nestedSortStruct(UGMref,'Start',1); % Sorts by ascending Start Index of hits
for c = 1:length(UGMref) % Finds lengths of hits
    UGMref(c).Total_Length=UGMref(c).End-UGMref(c).Start+1;
end
for c = 2:length(UGMref) % Finds distance between each hit
    UGMref_Distance(c) = UGMref(c).Start-UGMref(c-1).End;
end
UGMr_Seq_Mean_Distance = round(mean(UGMref_Distance),2);
UGMr = struct2table(UGMref);
UGMr_Total_Lengths = UGMr(:,6);
UGMr_Total_Mean_Length = round(mean(table2array(UGMr_Total_Lengths)),2);
UGMr = sortrows(UGMr,'Name','descend'); % Sorts by gene names
Gene_Total_Hits = height(UGMr);
writetable(UGMr,['Genome Repeats data/General/' Accession_List{a} 'GeneMatches.xlsx'])
if isempty(UGM_non_gene) == 0 % if there are hits with non gene matches
UGM_non_gene = nestedSortStruct(UGM_non_gene,'Start',1);
for c = 1:length(UGM_non_gene) % Finds lengths of hits
    UGM_non_gene(c).Total_Length=UGM_non_gene(c).End-UGM_non_gene(c).Start+1;
end
for c = 2:length(UGM_non_gene) % Finds distance between each hit
    UGM_non_gene_Distance(c) = UGM_non_gene(c).Start-UGM_non_gene(c-1).End;
end
UGM_non_gene_Seq_Mean_Distance = round(mean(UGM_non_gene_Distance),2);
UGMr_non_gene = struct2table(UGM_non_gene);
UGM_non_gene_Total_Lengths = UGMr_non_gene(:,4);
UGM_non_gene_Total_Mean_Length = round(mean(table2array(UGM_non_gene_Total_Lengths)),2);
Non_Gene_Total_Hits = height(UGMr_non_gene);
else 
    UGM_non_gene(1).Not_Found = 'None';
    UGMr_non_gene = struct2table(UGM_non_gene,'AsArray',true); 
    Non_Gene_Total_Hits = 0;
    UGM_non_gene_Total_Mean_Length = 0;
    UGM_non_gene_Seq_Mean_Distance = 0;
end
writetable(UGMr_non_gene,['Genome Repeats data/General/' Accession_List{a} 'NonGeneMatches.xlsx'])

% Calculations for One-Line output
Genome_length = length(part_sequence);

Expected_Hit_number = round(Genome_length/(4^(Palindrome_Number))); % Expected random number of genomes by probability

% Determine total number of general hits
if isempty(startIndex) == 0
Total_Hits = length(startIndex);
else
    Total_Hits = 0;
end

% Determine mean length of each hit
if isempty(startIndex) == 0
General_Total_Lengths = General_table(:,6);
General_Total_Mean_Length = round(mean(table2array(General_Total_Lengths)),2);
else
    General_Total_Mean_Length = 0;
end

% Determine distances between hits
if isempty(startIndex) == 0
for m = 1:length(startIndex)-1
    Distance(m) = startIndex(m+1)-endIndex(m);
end
General_Seq_Mean_Distance = round(mean(Distance),2);
else
    General_Seq_Mean_Distance = 0;
end

% Produce a one-line output for each genome
General_One_Line(a).Accession_Number = Accession_List{a};
General_One_Line(a).Genome_Name = Genome_name{1};
General_One_Line(a).Genome_Length = Genome_length;
General_One_Line(a).General_Hits = Total_Hits;
General_One_Line(a).General_Expected_Number_of_Hits = Expected_Hit_number;
General_One_Line(a).Gene_Hits = Gene_Total_Hits;
General_One_Line(a).Gene_Mean_Length = UGMr_Total_Mean_Length;
General_One_Line(a).Gene_Mean_Distance = UGMr_Seq_Mean_Distance;
General_One_Line(a).Non_Gene_Hits = Non_Gene_Total_Hits;
General_One_Line(a).Non_Gene_Mean_Length = UGM_non_gene_Total_Mean_Length;
General_One_Line(a).Non_Gene_Mean_Distance = UGM_non_gene_Seq_Mean_Distance;
General_One_Line(a).General_Mean_Length = General_Total_Mean_Length;
General_One_Line(a).General_Mean_Distance = General_Seq_Mean_Distance;
General_Counter = array2table(a);
writetable(General_Counter,'General_Counter.xlsx')
General_One_Line_Table = struct2table(General_One_Line);
writetable(General_One_Line_Table, 'Genome Repeats data/General/General_One_Line_All_Genomes.xlsx','Sheet','One-Line Data')
toc
end

% Calculations for Means for One_Line
Genome_List_Number = table2array(General_One_Line_Table(:,1));

% For gene hits
Column_Gene_Mean_Length = table2array(General_One_Line_Table(:,7));
Average_Gene_Mean_Length = mean(Column_Gene_Mean_Length);
General_One_Line(length(Genome_List_Number)+1).Gene_Mean_Length = Average_Gene_Mean_Length;

Column_Gene_Mean_Distance = table2array(General_One_Line_Table(:,8));
Average_Gene_Mean_Distance = mean(Column_Gene_Mean_Distance);
General_One_Line(length(Genome_List_Number)+1).Gene_Mean_Distance = Average_Gene_Mean_Distance;

% For non-gene hits
Column_Non_Gene_Mean_Length = table2array(General_One_Line_Table(:,10));
Average_Non_Gene_Mean_Length = mean(Column_Non_Gene_Mean_Length);
General_One_Line(length(Genome_List_Number)+1).Non_Gene_Mean_Length = Average_Non_Gene_Mean_Length;

Column_Non_Gene_Mean_Distance = table2array(General_One_Line_Table(:,11));
Average_Non_Gene_Mean_Distance = mean(Column_Non_Gene_Mean_Distance);
General_One_Line(length(Genome_List_Number)+1).Non_Gene_Mean_Distance = Average_Non_Gene_Mean_Distance;

% For all hits
Column_General_Mean_Length = table2array(General_One_Line_Table(:,12));
Average_General_Mean_Length = mean(Column_General_Mean_Length);
General_One_Line(length(Genome_List_Number)+1).General_Mean_Length = Average_General_Mean_Length;

Column_General_Mean_Distance = table2array(General_One_Line_Table(:,13));
Average_General_Mean_Distance = mean(Column_General_Mean_Distance);
General_One_Line(length(Genome_List_Number)+1).General_Mean_Distance = Average_General_Mean_Distance;
General_One_Line_Table = struct2table(General_One_Line);
writetable(General_One_Line_Table,'Genome Repeats data/General/General_One_Line_All_Genomes.xlsx','Sheet','One-Line Data')

% Write documentation with General One Line Table
General_One_Line_Documentation.General_Expression_Used = General_Palindrome;
General_One_Line_Documentation.General_Expression_Description = ['Expression of a palindrome sequence with any' num2str(Palindrome_Number) '-bp (number specified in Palindrome number) palindrome sections on both ends (inverted repeats) with  ' num2str(Min_Middle_Section) '-' num2str(Max_Middle_Section) ' bp in between the palindrome repeats (Defined in Middle_Section parameters)'];
General_One_Line_Documentation.Palindrome_Number = num2str(Palindrome_Number);
General_One_Line_Documentation.Min_Middle_Section_bp = num2str(Min_Middle_Section);
General_One_Line_Documentation.Max_Middle_Section_bp = num2str(Max_Middle_Section);
General_One_Line_Documentation.Step_1 = 'Search NCBI assembly database for bacterial reference genomes and returns accession numbers using function GetChromoAccNumfromAssemSearch.m';
General_One_Line_Documentation.Step_2 = 'Downloads genome sequence from NCBI database based on accession numbers using downloadGenomeKEH.m';
General_One_Line_Documentation.Step_3 = 'Using MATLAB regexp function to locate indexes in genome with a sequence with palindromes on both ends - expression found in "General Expression Used"';
General_One_Line_Documentation.Step_4 = 'Genes are identified using the analyze_feature_table.m to determine locations of annotated genes in genome';
General_One_Line_Documentation.Step_5 = 'Regexp match indexes are compared with those from annotated genes to determine which matches begin, end, or fully located in which genes - found in GeneMatches.xlsx';
General_One_Line_Documentation.Step_6 = 'Matches without a gene are separated and tabulated into NonGeneMatches.xlsx';
General_One_Line_Documentation.Step_7 = 'Descriptions of Parameters in Table explained:';
General_One_Line_Documentation.Genome_Name = 'Name of Species from NCBI based on accession number';
General_One_Line_Documentation.Genome_Length = 'Number of base pairs in whole genome';
General_One_Line_Documentation.General_Hits = 'Number of sequence matches based on the regexp function';
General_One_Line_Documentation.Gene_Hits = 'Number of sequence matches that are located within or partially within annotated genes - Specific details in GeneMatches.xlsx for each genome';
General_One_Line_Documentation.Non_Gene_Hits = 'Number of sequence matches that are not located within or partially within annotated genes - Specific details in NonGeneMatches.xlsx for each genome';
General_One_Line_Documentation.Gene_Mean_Length = 'Average number of bp for all "gene" matches - number in last row designates the average mean length for all genomes in this list';
General_One_Line_Documentation.Gene_Mean_Distance = 'Average distance in bp between each "gene" match - number in last row designates the average mean distance for all genomes in this list';
General_One_Line_Documentation.Non_Gene_Mean_Length = 'Average number of bp for all "non-gene" matches - number in last row designates the average mean length for all genomes in this list';
General_One_Line_Documentation.Non_Gene_Mean_Distance = 'Average distance in bp between each "non-gene" match - number in last row designates the average mean distance for all genomes in this list';
General_One_Line_Documentation.General_Mean_Length = 'Average number of bp for all matches - number in last row designates the average mean length for all genomes in this list';
General_One_Line_Documentation.General_Mean_Distance = 'Average distance in bp between each match - number in last row designates the average mean distance for all genomes in this list';
General_One_Line_Doc = struct2table(General_One_Line_Documentation);
writetable(General_One_Line_Doc, 'Genome Repeats data/General/General_One_Line_All_Genomes.xlsx','Sheet','Documentation')
end
function [part_sequence,startIndex,endIndex,General_table,Genome_Name,General_Palindrome,Palindrome_Number,Min_Middle_Section,Max_Middle_Section] = general_repeat(Accession_Number)

% Input Parameters for Palindrome Regular Expressions
Palindrome_Number = 8; % Number of bp in palindrome section
Min_Middle_Section = 10; % Minimum number of bp in middle section
Max_Middle_Section = 30; % Maximum number of bp in middle section

% Expressions for various palindrome sequences
General_Palindrome = ['(.{' num2str(Palindrome_Number) '}).{' num2str(Min_Middle_Section) ',' num2str(Max_Middle_Section) '}(??@palindrome_flip($1))']; % Finds general palindromes with defined palindrome number, middle section parameters

warning('off', 'MATLAB:mode:EmptyInput') %suppresses unimportant warning 
fid=fileread(['NCBI data/Sequence/' Accession_Number '.sq']);
fid1=fopen(['NCBI data/Sequence/' Accession_Number '.sq']);
seq = textscan(fid1, '%s %*[^\n]','HeaderLines',1);
sequence=char(seq{1}); %saves the sequence without the header for later use

% Name of Prokaryote
Genome_name_Expression = [Accession_Number '(.\d\s(\w+\s\w+)|\s(\w+\s\w+)|\s\[\w+\s\w+\])'];
Genome_name = regexp(fid,Genome_name_Expression,'tokens');
Genome_Name = Genome_name{1};
part_sequence = sequence;

% Clear General_Repeat structure to prevent future overlaps
General_Repeat = [];

% To test with general palindrome expression
[startIndex, endIndex] = regexp(part_sequence,General_Palindrome);
[matches,tokens] = regexp(part_sequence,General_Palindrome,'match','tokens');
for i = 1:length(startIndex)
    General_Repeat(i).Accession_Number = Accession_Number;
    General_Repeat(i).Genome_Name = Genome_Name;
    General_Repeat(i).StartIndex=startIndex(i);
    General_Repeat(i).EndIndex=endIndex(i);
    General_Repeat(i).Total_Seq=matches(i);
    General_Repeat(i).Total_Length=endIndex(i)-startIndex(i)+1;
    General_Repeat(i).Palindrome_Seq=tokens{i}{1};
    General_Repeat(i).Palindrome_Length = length(General_Repeat(i).Palindrome_Seq);
end
General_Repeat = nestedSortStruct(General_Repeat,'Palindrome_Length',1); % function to sort the palindrome length from smallest to biggest
General_Repeat(1).Reg_Expression = 'Palindromic Seq with 8 bp with 10-30 bp in between';
General_Repeat(2).Reg_Expression = General_Palindrome;
General_table = struct2table(General_Repeat);
writetable(General_table, ['Genome Repeats data/General/' Accession_Number 'gen_palindrome.xlsx'])
end
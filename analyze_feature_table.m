function [FeatureTable,Gene] = analyze_feature_table(Accession_Number)
feature_table=fileread(['NCBI data/Feature_Table/' Accession_Number '.ft']);
feat_expression = '(\d+)\t(\d+)\t(\w+)\n'; % General feature expression
descriptions = regexp(feature_table,feat_expression,'tokens');

% Expressions for Gene Features
Gene_Expression = '(\d+)\t(\d+)\t\w+\n\t\t\tgene\t(\w{3,4})'; % Primary label for gene
Gene_Expression1 = '(\d+)\t(\d+)\tgene\n(\d+)\t(\d+)\n\t\t\tgene\t(\w{3,4})'; % Label if gene has two regions
Gene_Expression2 = '(\d+)\t(\d+)\tgene\n(\d+)\t(\d+)\n(\d+)\t(\d+)\n\t\t\tgene\t(\w{3,4})'; % Label if gene has three regions
Gene_Descriptions = regexp(feature_table,Gene_Expression,'tokens');
Gene_Descriptions1 = regexp(feature_table,Gene_Expression1,'tokens');
Gene_Descriptions2 = regexp(feature_table,Gene_Expression2,'tokens');

% General Feature Table (denoted as 'F')
for i = 1:length(descriptions)
    FeatureTable(i).Start=str2double(descriptions{i}{1});
    FeatureTable(i).End=str2double(descriptions{i}{2});
    FeatureTable(i).Type=descriptions{i}{3};
end

% Formation of Gene Feature Table(denoted as 'G')
if isempty(Gene_Descriptions) == 0 % If gene names are in feature table
for i = 1:length(Gene_Descriptions)
Gene(i).Start = str2double(Gene_Descriptions{i}{1});
Gene(i).End = str2double(Gene_Descriptions{i}{2});
Gene(i).Name = Gene_Descriptions{i}{3};
end
for i = 1:length(Gene_Descriptions1)
    Gene(i+length(Gene_Descriptions)).Start = str2double(Gene_Descriptions1{i}{1});
    Gene(i+length(Gene_Descriptions)+length(Gene_Descriptions1)).Start = str2double(Gene_Descriptions1{i}{3});
    Gene(i+length(Gene_Descriptions)).End = str2double(Gene_Descriptions1{i}{2});
    Gene(i+length(Gene_Descriptions)+length(Gene_Descriptions1)).End = str2double(Gene_Descriptions1{i}{4});
    Gene(i+length(Gene_Descriptions)).Name = Gene_Descriptions1{i}{5};
    Gene(i+length(Gene_Descriptions)+length(Gene_Descriptions1)).Name = Gene_Descriptions1{i}{5};
end 
for i = 1:length(Gene_Descriptions2)
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)).Start = str2double(Gene_Descriptions2{i}{1});
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)+length(Gene_Descriptions2)).Start = str2double(Gene_Descriptions2{i}{3});
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)+2*length(Gene_Descriptions2)).Start = str2double(Gene_Descriptions2{i}{5});
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)).End = str2double(Gene_Descriptions2{i}{2});
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)+length(Gene_Descriptions2)).End = str2double(Gene_Descriptions2{i}{4});
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)+2*length(Gene_Descriptions2)).End = str2double(Gene_Descriptions2{i}{6});
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)).Name = Gene_Descriptions2{i}{7};
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)+length(Gene_Descriptions2)).Name = Gene_Descriptions2{i}{7};
    Gene(i+length(Gene_Descriptions)+2*length(Gene_Descriptions1)+2*length(Gene_Descriptions2)).Name = Gene_Descriptions2{i}{7};
end 
end
if isempty(Gene_Descriptions) == 1 % If no gene names in feature table
   General_Gene_Expression = '(\d+)\t(\d+)\tgene\n'; 
   General_Gene_Descriptions = regexp(feature_table,General_Gene_Expression,'tokens');
   for i = 1:length(General_Gene_Descriptions)
       Gene(i).Start = str2double(General_Gene_Descriptions{i}{1});
       Gene(i).End = str2double(General_Gene_Descriptions{i}{2});
       Gene(i).Name = 'Gene';
   end
end
end
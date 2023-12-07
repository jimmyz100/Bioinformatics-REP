function Gene_Name_Match = gene_locator(Start,End,Gene_Feature,Accession_Number) % Gene_Feature defined in analyze_feature_table function
% Start and End defined in excel table

% Set up of data structure
Gene_Name_Match(length(Start)).Accession_Number = [];
Gene_Name_Match(length(Start)).Start = [];
Gene_Name_Match(length(Start)).End = [];
Gene_Name_Match(length(Start)).Name = [];
Gene_Name_Match(length(Start)).Match_Attribute = [];

% Scan which indexes correspond with gene
for a = 1:length(Start)
    Gene_Name_Match(a).Accession_Number = Accession_Number;
    Gene_Name_Match(a).Start = Start(a);
    Gene_Name_Match(a).End = End(a);
    for i = 1:length(Gene_Feature)
        if Start(a) >= Gene_Feature(i).Start && End(a) <= Gene_Feature(i).End
            Gene_Name_Match(a).Name = Gene_Feature(i).Name; 
            Gene_Name_Match(a).Match_Attribute = 'Full';
        end
    end
end
for a = 1:length(Start)
    for i = 1:length(Gene_Feature)
    if isempty(Gene_Name_Match(a).Name) == 1
        if Start(a) <= Gene_Feature(i).End && End(a) >= Gene_Feature(i).End
            Gene_Name_Match(a).Name = Gene_Feature(i).Name;
            Gene_Name_Match(a).Match_Attribute = 'Begins';
        end
    end
    end
end
for a = 1:length(Start)
    for i = 1:length(Gene_Feature)
    if isempty(Gene_Name_Match(a).Name) == 1
        if Start(a) <= Gene_Feature(i).Start && End(a) >= Gene_Feature(i).Start
            Gene_Name_Match(a).Name = Gene_Feature(i).Name;
            Gene_Name_Match(a).Match_Attribute = 'Ends';
        end
    end
    end
end
end
function [AssemblyInfo] = getAssemblyInfo2( AssemblyAcc )
%parses assembly description and extracts information including
%    Plasmids: Count, length and Accession number
%    Chromosomes: Count, length and Accession number
%    RefSeq category (e.g. Reference genome)
%    Assembly level (e.g. complete genome)
%    Phylum

%Note: This only works for GCF prefix assembly accession numbers currently
try
   name=getAssNamefromAccNum(AssemblyAcc);   %gets assembly name from assembly accession
   space=regexp(name,' ');
   if length(space)>0    %replace spaces with underscores
       name=regexprep(name,' ','_');
   end
   slash=regexp(name,'/');
   if length(slash)>0
       name=regexprep(name,'/','_');
   end
   comma=regexp(name,',');   %delete comma from name
   if length(comma)>0
       name=regexprep(name,',','');
   end
   
   
   %fname=[AssemblyAcc '.assembly.txt'];
   URL='http://ftp.ncbi.nlm.nih.gov/genomes/all/';
   patF1='(GCF|GCA)_(\d{3})(\d{3})(\d{3})';
   tok=regexp(AssemblyAcc,patF1,'tokens');
   URL=[URL AssemblyAcc(1:3) '/' tok{1}{2} '/' tok{1}{3} '/' tok{1}{4} '/' AssemblyAcc '_' name '/' AssemblyAcc '_' name '_assembly_report.txt'];
   data=webread(URL);

   pat_level='Assembly level: (.*?)\n';   %grabs assembly level from assembly summary file 
   pat_p='Plasmid.{1,20}=\s+(\w+.\w+)\s+Primary Assembly\s+(\d+)\s+na';
   pat_c='Chromosome.{1,20}=\s+(\w+.\w+)\s+Primary Assembly\s+(\d+)\s+na';
   pat_refseqcat='RefSeq category: (.*?)\n';
   tok_p=regexp(data,pat_p,'tokens');    %Plasmids
   tok_c=regexp(data,pat_c,'tokens');   %Chromosomes
   tok_cat=regexp(data,pat_refseqcat,'tokens');%grabs RefSeq category from assembly summary file
   tok_level=regexp(data,pat_level,'tokens');  %grabs assembly level

   for i=1:length(tok_p)
       PlasmidAccession{i}=tok_p{i}{1};
       PlasmidLength{i}=tok_p{i}{2};
   end
   
   for i=1:length(tok_c)
       ChromoAccession{i}=tok_c{i}{1};
       ChromoLength{i}=tok_c{i}{2};
   end
   
   AssemblyInfo.AssemblyAccession=AssemblyAcc;
   AssemblyInfo.Name=name;
   
   if length(tok_p)==0
       AssemblyInfo.PlasmidAccession={};
       AssemblyInfo.PlasmidLength={'0'};
   else
       AssemblyInfo.PlasmidAccession=PlasmidAccession;
       AssemblyInfo.PlasmidLength=PlasmidLength;
   end
    
   if length(tok_c)==0
       AssemblyInfo.ChromoAccession={};
       AssemblyInfo.ChromoLength=0;
   else
       AssemblyInfo.ChromoAccession=ChromoAccession;
       AssemblyInfo.ChromoLength=ChromoLength;
   end
   
   AssemblyInfo.NumPlasmids=length(tok_p);
   AssemblyInfo.NumChromosomes=length(tok_c);
   if length(tok_cat)==0
      AssemblyInfo.RefSeqCategory=''; 
   else    
      AssemblyInfo.RefSeqCategory=tok_cat{1}{1};
   end
   
   if length(tok_level)==0
      AssemblyInfo.AssemblyLevel=''; 
   else    
      AssemblyInfo.AssemblyLevel=tok_level{1}{1};
   end
   
   AssemblyInfo.Phylum=getPhylum(AssemblyAcc);
   AssemblyInfo.error=0;
catch
    [AssemblyAcc  '    :DID NOT READ']
    AssemblyInfo.error=999;
end
end
%script that searches the assembly database using your search string
%example search string 'bacteria[orgn] AND "reference genome"[refseq category]'
%The example search string yields 120 assembly accession numbers, some of
%which have more than one chromosome
%
%Program returns List which contains a list of all of the chromosome
%accession numbers
AssAcc=getAccNumfromSearch('bacteria[orgn] AND "reference genome"[refseq category]');

clear List
N=AssAcc.Count;
iList=1;
iGCF=0;
for i=1:N
    %if strcmp(AssAcc.Accession{i},'GCF_900010755.1') || strcmp(AssAcc.Accession{i},'GCF_900065885.1')   %problem with comma in name, comma not in data file name
        
    %else
        if strcmp(AssAcc.Accession{i}(1:3),'GCF')
           iGCF=iGCF+1; 
           AssInfo=getAssemblyInfo2(AssAcc.Accession{i});
           if AssInfo.error~=999
               for j=1:AssInfo.NumChromosomes
                  List{iList,1}=AssInfo.Name;
                  List{iList,2}=AssInfo.ChromoAccession{j};
                  List{iList,3}=AssInfo.AssemblyLevel;
                  iList=iList+1;
               end
           end
        end
    %end
end
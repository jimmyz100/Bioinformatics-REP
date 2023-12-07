function [AccessionList]=getAccNumfromSearch(searchTerm)
%retieves accession numbers from Assembly database based on a NCBI search
%string. See example next which results in 4 accession numbers

%getAccNumfromSearch('acinetobacter AND "representative genome" AND "complete genome"')
options=weboptions('Timeout',15);
%first run search command and store output on NCBI server 
url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly';
url=[url '&term=' searchTerm '&usehistory=y'];
data=webread(url,options);   %reports back codes that allow access to stored information

%parse the needed information from the output
patCount='<eSearchResult><Count>(\d+)</Count>';
patWebEnv='<WebEnv>(.*)</WebEnv>';
patQueryKey='<QueryKey>(\d+)</QueryKey>';

tokCount=regexp(data,patCount,'tokens');
Count=tokCount{1}{1};
tokWebEnv=regexp(data,patWebEnv,'tokens');
WebEnv=tokWebEnv{1}{1};
tokQueryKey=regexp(data,patQueryKey,'tokens');
QueryKey=tokQueryKey{1}{1};

%Now set up a loop to access the data summary for each of the ID's obtained
%above. The loop is used to minimize the number of calls to esummary. This
%is important only when the search results in a large number of hits. This 
%loop structure is recommeded by NCBI 
patAcc = '<AssemblyAccession>(GC[AF]_\d+.\d+)</AssemblyAccession>';
patName = '<AssemblyName>(.+?)</AssemblyName>';
patTitle='<Organism>(.*?)</Organism>';
retdel=500;   %arbitrary block size
ict1=1;
ict2=1;
ict3=1;
options=weboptions('Timeout',15);
for i=1:retdel:str2double(Count)
    url2='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=';
    url2 = [url2 'assembly&WebEnv=' WebEnv '&query_key=' QueryKey];
    url2 = [url2 '&retstart=' num2str(i) '&retmax=' num2str(retdel)];
    data2=webread(url2,options);  %Reports back retdel
    
    match=regexp(data2,patAcc,'tokens');
    for j=1:length(match)
        Acc{ict1}=match{j}{1};
        ict1=ict1+1;
    end
    
    match2=regexp(data2,patTitle,'tokens');
    for j=1:length(match2)
        Title{ict2}=match2{j}{1};
        ict2=ict2+1;
    end
    
    match3=regexp(data2,patName,'tokens');
    for j=1:length(match3)
        Name{ict3}=match3{j}{1};
        ict3=ict3+1;
    end
    
end
 if ict1~=ict2
     AccessionList.Count=-999;
 end
 AccessionList.Count=ict1-1;
 AccessionList.Accession=Acc;
 AccessionList.Title=Title;
 AccessionList.Name=Name;
 
end


function AssemblyName=getAssNamefromAccNum(Accession)
    %retieves AssemblyName from assembly accession number

    %getAssNamefromAccNum('GCF_001932555.1')

    options=weboptions('Timeout',15);
    %first run esearch command with accession number to get id
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly';
    url=[url '&term=' Accession ];
    data=webread(url,options);   %reports back assembly id

    %parse the needed information from the output
    patId='<Id>(\d+)</Id>';
    tok=regexp(data,patId,'tokens');
    Id=tok{1}{1};

    url2='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=';
    url2 = [url2 'assembly&id=' Id];
    data2=webread(url2,options);  %Reports back the assembly summary
    
    patName='<AssemblyName>(.+?)</AssemblyName>';
    tok=regexp(data2,patName,'tokens');
    AssemblyName=tok{1}{1};     %assembly name

end


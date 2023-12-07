function downloadGenomeKEH(Accession)
    %got rid of 'wt' on write commands --> no extra CR in files
    %Check to see if database folder exists already, create if it does not
    if 7 ~= exist('NCBI data')
        mkdir('NCBI data');
    end
    if 7 ~= exist('NCBI data/Sequence')
        mkdir('NCBI data/Sequence');
    end
    if 7 ~= exist('NCBI data/Feature_Table')
        mkdir('NCBI data/Feature_Table');
    end


    %Check if feature table has been downloaded previously
    if exist(['NCBI data/Feature_Table/' Accession '.ft']) ~= 2
        URL=['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=' Accession '&rettype=ft'];
        feature=urlread(URL);
        fid=fopen(['NCBI data/Feature_Table/' Accession '.ft'],'w');
        fprintf(fid, '%s', feature);
        fclose(fid);
    end

    %Check if sequence has been downloaded previously
    if exist(['NCBI data/Sequence/' Accession '.sq']) ~= 2
        URL=['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=' Accession '&rettype=fasta'];
        Seq=fastaread(URL);
        fid=fopen(['NCBI data/Sequence/' Accession '.sq'],'w');
        fprintf(fid, ['>' Seq.Header '\n']);
        fprintf(fid, [Seq.Sequence '\n']);
        fclose(fid);
    end
end
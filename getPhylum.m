function t = getPhylum(AssemblyAccession)
%Given Assembly Accession #, report back phylum
doc2 = urlread(['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&amp;term=' AssemblyAccession]);
tok = regexp(doc2, '<Id>([0-9]+)</Id>','tokens');

doc2 = urlread(['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=' tok{1}{1}]);
tok = regexp(doc2, '<Taxid>([0-9]+)</Taxid>','tokens');

doc = urlread(['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=' tok{1}{1}]);
pat='<ScientificName>(\w*)</ScientificName>\n[ ]*<Rank>phylum';
tok = regexp(doc, pat,'tokens');

if ~isempty(tok)
   t = tok{1}{1};
else
    t='';
end
end
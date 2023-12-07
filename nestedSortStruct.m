function [sortedStruct index] = nestedSortStruct(aStruct, fieldNamesCell, directions)

if ~isstruct(aStruct)
    error('first input supplied is not a struct.')
end % if

if sum(size(aStruct)>1)>1 % if more than one non-singleton dimension
    error('I don''t want to sort your multidimensional struct array.')
end % if

if ~iscell(fieldNamesCell)
    if isfield(aStruct, fieldNamesCell) % if fieldNamesCell is a simple string of a valid fieldname
        [sortedStruct index] = sortStruct(aStruct, fieldNamesCell);
        return
    else
        error('second input supplied is not a cell array or simple string of a fieldname.')
    end % if isfield
end % if ~iscell

if ~isfield(aStruct, fieldNamesCell)
    for ii=find(~isfield(aStruct, fieldNamesCell))
        fprintf('%s is not a fieldname in the struct.\n', fieldNamesCell{ii})
    end % for
    error('at least one entry in fieldNamesCell is not a fieldname in the struct.')
end % if

fieldFlag = 0;
for ii=1:length(fieldNamesCell)
    fieldEntry = aStruct(1).(fieldNamesCell{ii});
    if ~( ((isnumeric(fieldEntry) || islogical(fieldEntry)) && numel(fieldEntry)==1) || ischar(fieldEntry) )
        fprintf('%s is not a valid fieldname by which to sort.\n', fieldNamesCell{ii})
        fieldFlag = 1;
    end % if
end % for ii

if fieldFlag
    error('at least one fieldname is not a valid one by which to sort.')
end

if nargin < 3 % if directions doesn't exist
    directions = ones(1, length(fieldNamesCell));
else % check directions if it does exist
    if ~(isnumeric(directions) && all(ismember(directions, [-1 1])))
        error('directions, if given, must be a single number or a vector with 1 (ascending) and -1 (descending).')
    end % if ~(...
    
    if numel(directions)==1
        directions = directions * ones(1, length(fieldNamesCell)); % create vector from single element
    elseif length(fieldNamesCell)~=length(directions)
        error('fieldNamesCell and directions vector are different lengths.')
    end % if numel...
end % if exist...

[dummy fieldNamesIdx] = ismember(fieldNamesCell, fieldnames(aStruct));

aCell = squeeze(struct2cell(aStruct))';

[sortedCell index] = sortrows(aCell, fieldNamesIdx .* directions);

sortedStruct = aStruct(index); % apply the index to the struct array
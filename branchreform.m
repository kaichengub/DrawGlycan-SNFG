function reformed = branchreform(allsgp)
% BRANCHREFORM: discriptions
%
% Syntax:
% N/A
%
% Input:
% N/A
%
% Output:
% N/A
%
% Note:
% N/A
%
% Example:
% N/A
%
% Children function: 
% N/A
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2016, Research Foundation for State University of New York. All rights reserved
%

reformed = cell(size(allsgp));
for i = 1:numel(allsgp)
    thisgly = allsgp{i};
    thisgly = strrep(thisgly,'[','');
    thisgly = strrep(thisgly,']','');
     thisglynocurl = thisgly;
    [optvalstart,optvalend] = regexp(thisgly,'["''].*?["'']','start','end');
    for j = 1:length(optvalstart)
        thisglynocurl(optvalstart(j):optvalend(j)) = ['"',repmat('?',1,optvalend(j)-optvalstart(j)-1),'"'];
    end
    indtemp = strfind(thisglynocurl,'{');
    levelindex = zeros(2,length(thisgly));
    indtemp2 = strfind(thisglynocurl,'}');
    if ~ismember(1,indtemp)
        error('Did you forget something?');
    end
    levelindex(1,indtemp) = 1;
    levelindex(1,indtemp2) = -1;
    levelindex(2,1) = 1;
    for j = 2:size(levelindex,2)
        levelindex(2,j) = levelindex(2,j-1) + levelindex(1,j);
    end
    wholestr = zeros(1,length(thisgly));
    wholestr(indtemp) = 1;
    wholestr(indtemp2) = 1;
    allbond = cell(sum(wholestr)/2,1);
    writeind = 1;
    tempstr = '';
    readind = 1;
    while readind <= length(wholestr)
        if wholestr(readind) ~= 1  % gather character to form the monosac. string
            tempstr = [tempstr,thisgly(readind)];
        else
            if ~isempty(tempstr)
                bondinfo = regexp(tempstr,'[a-z?][\d?]+-[\d?]+','match');
                allbond{writeind} = bondinfo{1};
                tempstr = '';
                writeind = writeind + 1;
            end
        end
        readind = readind + 1;
    end

    %% calculate distance of each monosac
    letterindex = zeros(1,length(thisgly));
    letterindex(regexp(thisglynocurl,'[^{}]')) = 1;
    distance = letterindex.*levelindex(2,:);
    distance = distance(indtemp+1);  % all monosac's
    if length(indtemp) > 1
        for j = 2:length(indtemp)
            letterindex(indtemp(j):end) = letterindex(indtemp(j):end)/(j-1)*j;
        end
    end
    bondmap = zeros(length(distance));
    readind = 1;
    while readind < length(distance)
        if (distance(readind + 1) - distance(readind)) == 1
            bondmap(readind,readind + 1) = 1;  % consecutive numbers indicate bond
            readind = readind + 1;
        elseif (distance(readind + 1) - distance(readind)) < 1  % if chain is broken, go back to find its fork point
            thisind = distance(readind + 1);  % where it's broken
            itsforkpt = find(distance(1:readind) == thisind - 1,1,'last');  % where is the fork point
            bondmap(itsforkpt,readind + 1) = 1;  % mark this bond
            readind = readind + 1;  % keep going on
        end
    end
    isfork = find(sum(bondmap,2) > 1);
    if any(isfork)
        for j = length(isfork):-1:1
            forkchildren = find(bondmap(isfork(j),:));
            branches = cell(size(forkchildren));
            bonds = cell(length(forkchildren),2);
            branchpos = zeros(length(forkchildren),2);
            for k = 1:length(forkchildren)
                branches{k} = glytreetracker(bondmap,forkchildren(k),[],'down');
                branchpos(k,1) = min(find(letterindex == forkchildren(k))) - 1;
                branchpos(k,2) = find(levelindex(2,branchpos(k,1):end) == levelindex(2,branchpos(k,1))-1,1) + branchpos(k,1)-1;
                tempbond = allbond{forkchildren(k)};
                carbon = regexp(tempbond,'[a-z?]','match');
                nrnum = regexp(tempbond,'[0-9?]+','match');
                bonds{k,1} = carbon{1};
                bonds{k,2} = nrnum{end};
            end
            branchseq = 1:size(bonds,1);
            [~,ind1] = sort(bonds(:,1));
            bonds = bonds(ind1,:);
            branchseq = branchseq(ind1);
            [~,ind2] = sort(bonds(:,2));
            branchseq = branchseq(ind2);
            thisglyhead = thisgly(1:min(min(branchpos))-1);
            thisglytail = thisgly(max(max(branchpos))+1:end);
            thisglymiddle = thisgly(min(min(branchpos)):max(max(branchpos)));
            letterindhead = letterindex(1:min(min(branchpos))-1);
            letterindtail = letterindex(max(max(branchpos))+1:end);
            letterindmiddle = letterindex(min(min(branchpos)):max(max(branchpos)));
            lvlindhead = levelindex(:,1:min(min(branchpos))-1);
            lvlindtail = levelindex(:,max(max(branchpos))+1:end);
            lvlindmiddle = levelindex(:,min(min(branchpos)):max(max(branchpos)));
            minmsind = cellfun(@min,branches);
            maxmsind = cellfun(@max,branches);
            bondmaphead = bondmap(1:min(minmsind)-1,:);
            bondmaptail = bondmap(max(maxmsind)+1:end,:);
            branchpos = branchpos-min(min(branchpos))+1;
            newglymiddle = [];
            newletterindmiddle = [];
            newlvlindmiddle = [];
            newbondmapmiddle = [];
            for k = 1:length(branchseq)
                newglymiddle = [newglymiddle,thisglymiddle(branchpos(branchseq(k),1):branchpos(branchseq(k),2))];
                newletterindmiddle = [newletterindmiddle,letterindmiddle(branchpos(branchseq(k),1):branchpos(branchseq(k),2))];
                newlvlindmiddle = [newlvlindmiddle,lvlindmiddle(:,branchpos(branchseq(k),1):branchpos(branchseq(k),2))];
                newbondmapmiddle = [newbondmapmiddle;bondmap(branches{branchseq(k)},:)];
            end
            thisgly = [thisglyhead,newglymiddle,thisglytail];
            letterindex = [letterindhead,newletterindmiddle,letterindtail];
            levelindex = [lvlindhead,newlvlindmiddle,lvlindtail];
            bondmap = [bondmaphead;newbondmapmiddle;bondmaptail];
        end
    end
    reformed{i} = thisgly;
end
end
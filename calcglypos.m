function plotinfoout = calcglypos(sgpseq,options,specialoptions)
% function [msposout,bondmapout,allmonosacout,directionseqout] = calcglypos(allsgp,options)
% CALCGLYPOS: calculate monosaccharide positions of glycan(s)
%
% Syntax:
% [msposout,bondmapout,allmonosacout,alllinkout,directionseqout] = calcglypos(allsgp,options)
%
% Input:
% allsgp: n x 1 cell array of strings, the glycans to be analyzed.
% options: structure, contains options for drawing glycans. See DRAWGLYCAN
% for details.
%
% Output:
% msposout: n x 1 cell array of m x 2 numerical matrices, monosaccharide
% positions in the glycan. n equals to the number of glycans in  'allsgp',
% m equals to the number of monosaccharides in the corresponding glycans.
% bondmapout: n x 1 cell array of m x m numerical matrices, the linkage
% information of monosaccharides in glycan.
% allmonosacout: n x 1 cell array of m x 1 cell array of strings, the name
% of the monosaccharides in glycan.
% alllinkout: n x 1 cell array of m x 1 cell array of strings, the
% glycosidic bond information of glycans.
% directionseqout: n x 1 cell array of m x 1 numerical matrices, the
% orientation of monosaccharides.
%
% Note:
% Solo perpendicular monosaccharides will be
% plotted perpendicularly, e.g. the top vertex of Fuc will point to the
% left if the orientation of the glycan is 'up', when Fuc is the only
% monosaccharide in the glycan.
%
% Example:
% N/A. Run examples in DRAWGLYCAN and set breakpoints.
%
% Children function:
% DRAWRAWTREE, BRANCHEQUALIZER, SMARTPMDISPENSER, BONEADDPM, ADDSTUBPM

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2016, Research Foundation for State University of New York. All rights reserved
%

%% OPTIONAL: reform the glycan structure based on linkage information, sort the branches by linkage number
%% TO USE: de-comment the following command:
% allsgp = branchreform(allsgp);

alllinkout = regexp(sgpseq,DrawGlycanPara.regexp_monosaclinkage,'match');
for i = 1:length(alllinkout)
    alllinkout{i} = alllinkout{i}(2:end-1);
end
%% Curly-bracket handling: remove monosac. outside curly brac.
%% Start calculation

%% sep. monosac. and bond
thisgly = sgpseq;
thisgly = strrep(thisgly,'[','');
thisgly = strrep(thisgly,']','');
indtemp = strfind(thisgly,'{');
levelindex = zeros(2,length(thisgly));
indtemp2 = strfind(thisgly,'}');
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
readind = 1;
allms = cell(sum(wholestr)/2,1);
writeind = 1;
tempstr = '';
while readind <= length(wholestr)
    if wholestr(readind) ~= 1  % gather character to form the monosac. string
        tempstr = [tempstr,thisgly(readind)];
    else
        if ~isempty(tempstr)
            bondinfo = strfind(tempstr,')');
            if ~isempty(bondinfo)  % using parenthesis
                allms{writeind} = tempstr(max(bondinfo) + 1:end);
            else  % No parenthesis means there's only monosac.
                allms{writeind} = tempstr;
            end
            tempstr = '';
            writeind = writeind + 1;
        end
    end
    readind = readind + 1;
end
%% calculate distance of each monosac
letterindex = zeros(1,length(thisgly));
letterindex(regexp(thisgly,'[^{}]')) = 1;
distance = letterindex.*levelindex(2,:);
distance = distance(indtemp+1);  % all monosac's


ispm = find(ismember(lower(allms),lower(options.perpendicularmonosac)));
%% customized perpendicular monosac.
if ~isempty(specialoptions)
    spoptfld = fieldnames(specialoptions);
    if ismember('P',upper(spoptfld))
        specialpm = specialoptions.(spoptfld{ismember(upper(spoptfld),'P')});
        specialpm = cell2mat(specialpm(:,2));
    else
        specialpm = [];
    end
else
    specialpm = [];
end

ispm = [ispm;specialpm(:)];
directionseq = zeros(size(distance));


%% Decide working mode
if isempty(ispm)
    drawglycanmode = 1;  % bone only
else
    if ismember(1,ispm)
        drawglycanmode = 2;  % PMS only
    else
        drawglycanmode = 3;  % bone+PMS
    end
end


%% build adjacency matrix "bondmap"
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


%% Handling 3 different situations
switch drawglycanmode
    case 1
        mspos = drawrawtree(bondmap,distance,[]);
        mspos = branchequalizer(mspos,bondmap);
        
    case 2
        mspos = drawrawtree(bondmap,distance,[]);
        mspos = branchequalizer(mspos,bondmap);
        
    case 3
        bonebondmap = bondmap;  % bone part
        bonedistance = distance;
        pmchildlist = {};  % Each element in 'pmschildlist' is a sub-tree connected to the main structure,
        % the first monosac in each element is 'thispms'
        while ~isempty(ispm)  % retrieve the monosac index of PMS containing sub-trees
            thispm = ispm(find(ispm,1,'first'));
            tobetracked = thispm;
            temppmchild = [];
            while ~isempty(tobetracked)
                if any(bondmap(tobetracked(1),:))
                    tobetracked = [tobetracked find(bondmap(tobetracked(1),:))];
                end
                temppmchild = [temppmchild;tobetracked(1)];
                tobetracked(1) = [];
            end
            pmchildlist = [pmchildlist;temppmchild];
            ispm = setdiff(ispm,temppmchild);
        end
        parentms = zeros(size(pmchildlist));
        allpmchild = [];  % this is used later to draw the bone structure of glycan
        for j = 1:length(pmchildlist)
            parentms(j) = find(bondmap(:,pmchildlist{j}(1)));
            allpmchild = [allpmchild;pmchildlist{j}];
        end
        [parentmsnum,~,parentmsind] = unique(parentms);
        pmgroup = cell(max(parentmsind),2);  % The number of elements in 'pmsgroup' equals the number of
        for j = 1:size(pmgroup,1)
            pmgroup{j,1} = parentmsnum(j);  % 1st column contains the parent monosac
            pmgroup{j,2} = pmchildlist(parentmsind == j);  % 2nd column contains children PMS
        end
        bonebondmap(:,allpmchild) = 0;
        bonebondmap(allpmchild,:) = 0;
        bonedistance(allpmchild) = 0;
        tempbonebondmap = bonebondmap;
        tempbonebondmap(allpmchild,:) = [];
        tempbonebondmap(:,allpmchild) = [];
        mspos = drawrawtree(bondmap,distance,allpmchild);  % Here mspos contains all pos info of bone monosac
        bonemspos = mspos;
        bonemspos(allpmchild,:) = [];
        bonemspos = branchequalizer(bonemspos,tempbonebondmap);
        mspos(setdiff(1:length(distance),allpmchild),:) = bonemspos;
        %% start dealing with PM
        [pmposstorage,pmindstorage,directionseq] = smartPMdispenser(pmgroup,bondmap,distance,directionseq);
        allfork = find(sum(bonebondmap,2) > 1);
        if ~isempty(allfork)
            [mspos,safezone,PMroot] = boneaddPM(allfork,mspos,bonebondmap,bondmap,pmgroup,pmindstorage,pmposstorage,allpmchild);
            stub = min(find(sum(bonebondmap,2) > 1,1,'first')) - 1;  % if there is a fork
            mspos = addstubPM(mspos,bondmap,stub,bonebondmap,PMroot,pmposstorage,pmindstorage,safezone,directionseq);
        else
            safezone = find(bonedistance);
            stub = safezone(end);
            PMroot = cell2mat(pmgroup(:,1));
            [mspos,directionseq] = addstubPM(mspos,bondmap,stub,bonebondmap,PMroot,pmposstorage,pmindstorage,safezone,directionseq);
        end
end
mspos(:,2) = -mspos(:,2);  % glycan has been drawn upside down, flip it up.
mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);  % normalize all coordinates


%% Option: orientation
if strcmpi(options.orientation,'left')
    mspos = -mspos;
    mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);
elseif strcmpi(options.orientation,'up')
    tmp = mspos(:,1);
    mspos(:,1) = -mspos(:,2);
    mspos(:,2) = tmp;
    mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);
elseif strcmpi(options.orientation,'down')
    tmp = mspos(:,1);
    mspos(:,1) = mspos(:,2);
    mspos(:,2) = -tmp;
    mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);
end
plotinfoout.mspos = mspos;
plotinfoout.bondmap = bondmap;
plotinfoout.alllinkout = alllinkout(:);
plotinfoout.allms = allms;
plotinfoout.directionseq = directionseq(:);
%% pass info from specialoptions to plotinfoout
if ~isempty(specialoptions)
    spfieldnames = fieldnames(specialoptions);
    bonddesname = spfieldnames(ismember(spfieldnames,[DrawGlycanPara.intglybondinfo,DrawGlycanPara.intglymodinfo]));
    bonddescription = {};
    for i = 1:length(bonddesname)
        tempbonddes = specialoptions.(bonddesname{i});
        bonddescription = [bonddescription;tempbonddes(:,2),repmat(bonddesname(i),size(tempbonddes,1),1),tempbonddes(:,1)];
    end
    plotinfoout.bonddescription = bonddescription;
    spfieldnames = setdiff(spfieldnames,bonddesname);
    for i = 1:length(spfieldnames)
        plotinfoout.(spfieldnames{i}) = specialoptions.(spfieldnames{i});
    end
end
end


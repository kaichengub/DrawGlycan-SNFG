function [pmposstorage,pmindstorage,directionseq] = smartPMdispenser(pmsgroup,bondmap,distance,directionseq)
% SMARTPMDISPENSER: smart penpendicular monosac. structure dispenser. Helps
% arrange perpendicular glycan sub-structures within the bone structure
% 
% Syntax:
% [pmsposstorage,pmsindstorage,directionseq] = smartPMdispenser(pmsgroup,bondmap,distance,directionseq)
% 
% Input:
% pmsgroup: n x 2 cell array, 1st column contains the serial number of the
% monosac. that have PM attached to it, 2nd column contains the serial
% number of the corresponding PM's.
% bondmap: m x m numerical array, the linkage map of the glycan.
% distance: 1 x m numerical array, the distance of the monosac. in the
% glycan.
% directionseq: 1 x m numerical array, the orientation of the monosac. in
% the glycan.
% 
% Output:
% pmposstorage: n x 2 cell array of m x 2 numerical arrays, the positioning
% info. of the PM's. The 1st column % contains the position info for the PM's
% to be drawn upwards, 2nd column contains those to be drawn downwards.
% pmindstorage: n x 2 cell array of m x 1 numerical arrays, the serial
% numbers of the PM's to be drawn upwards or downwards. 1st column contains
% upwards ones, 2nd column contains downwards ones.
% directionseq: 1 x m numerical array, the orientation of the monosac.
% after upwards/downwards adjustment.
% 
% Note:
% If more than one PM is directly attached to the same monosac. in the bone
% structure, SMARTPMDISPENSER will automatically dispense the PM subtrees
% evenly around the root monosac.: the first points down, second points up,
% third points down again.
%
% Example:
% N/A. Set breakpoints in main program to see how this function works.
%
% Children function: 
% DRAWRAWTREE, BRANCHEQUALIZER
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2016, Research Foundation for State University of New York. All rights reserved
%


pmposstorage = cell(size(pmsgroup));
pmindstorage = cell(size(pmsgroup));
for i = 1:size(pmsgroup,1)  % one subtree at a time
    tmpsubtreeroot = pmsgroup{i,1};  % temp sub tree root
    if numel(pmsgroup{i,2}) ~= 1  % if there are more than 1 arm, need to consider distribution of arms
        uparms = 1:floor(numel(pmsgroup{i,2})/2);
        downarms = floor(numel(pmsgroup{i,2})/2)+1:numel(pmsgroup{i,2});
    else
        uparms = [];
        downarms = 1;
    end
    if ~isempty(uparms)
        alluparms = pmsgroup{i,2}(uparms);
%         perpenind = intersect(find(directionseq),cell2mat(alluparms));
%         directionseq(perpenind) = 180;
        tmparmcont = [];
        for j = 1:numel(alluparms)
            tmparmcont = [tmparmcont;alluparms{j}];
        end
        tmparmbondmap = bondmap([tmpsubtreeroot;tmparmcont],[tmpsubtreeroot;tmparmcont]);
        tmparmdistance = distance([tmpsubtreeroot;tmparmcont]);
        uparmmspos = drawrawtree(tmparmbondmap,tmparmdistance,[]);
        uparmmspos(:,1) = uparmmspos(:,1) - uparmmspos(1,1);
        uparmmspos(:,2) = uparmmspos(:,2) - uparmmspos(1,2);
        uparmmspos = branchequalizer(uparmmspos,tmparmbondmap);
        tmprow = -uparmmspos(:,1);
        uparmmspos(:,1) = uparmmspos(:,2);
        uparmmspos(:,2) = tmprow;
        uparmmspos = uparmmspos(2:end,:);
        pmposstorage{i,1} = uparmmspos;
        temppmsind = pmsgroup{i,2}(uparms);
        pmsind = [];
        for j = 1:numel(temppmsind)
            pmsind = [pmsind;temppmsind{j}];
        end
        pmindstorage{i,1} = pmsind;
    end
    if ~isempty(downarms)
        alldownarms = pmsgroup{i,2}(downarms);
        tmparmcont = [];
        for j = 1:numel(alldownarms)
            tmparmcont = [tmparmcont;alldownarms{j}];
        end
        tmparmbondmap = bondmap([tmpsubtreeroot;tmparmcont],[tmpsubtreeroot;tmparmcont]);
        tmparmdistance = distance([tmpsubtreeroot;tmparmcont]);
        downarmmspos = drawrawtree(tmparmbondmap,tmparmdistance,[]);
        downarmmspos(:,1) = downarmmspos(:,1) - downarmmspos(1,1);
        downarmmspos(:,2) = downarmmspos(:,2) - downarmmspos(1,2);
        downarmmspos = branchequalizer(downarmmspos,tmparmbondmap);
        tmprow = downarmmspos(:,1);
        downarmmspos(:,1) = downarmmspos(:,2);
        downarmmspos(:,2) = tmprow;
        downarmmspos = downarmmspos(2:end,:);
        pmposstorage{i,2} = downarmmspos;  % turn PMS subtree 90 deg
        temppmsind = pmsgroup{i,2}(downarms);
        pmsind = [];
        for j = 1:numel(temppmsind)
            pmsind = [pmsind;temppmsind{j}];
        end
        pmindstorage{i,2} = pmsind;
    end
end
end
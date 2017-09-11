function paintglycan(plotinfoout,options)
% PAINTGLYCAN: Draws monosaccharides using symbols, glycosidic bonds using lines,
% and also depicts glycan fragmentations
%
% Syntax:
% paintglycan(msposout,bondmapout,alllinkout,bondbreak,allmonosacout,directionseqout,options)
%
% Input:
% msposout: n x 1 cell array of m x 2 numerical array, the position info of
% glycans.
% bondmapout: n x 1 cell array of m x m numerical array, the linkage info
% of monosac. in glycans.
% alllinkout: n x 1 cell array of m x 1 cell array of strings, the
% glycosidic bond info of glycans.
% bondbreak: n x 1 cell array of 1 x 3 cell array. 1st elelment is
% fragmentation position, 2nd is type, 3rd is content.
% allmonosacout: n x 1 cell array of m x 1 cell array of strings, the
% monosac. in each glycan.
% directionseqout: n x 1 cell array of 1 x m numerical array, the
% orientation of monosac. in each glycan.
% options: structure, see help document in DRAWGLYPEP for detail.
%
% Output:
% Figure, a new figure must be created with 'hold' is set to 'on' prior to
% execution of the program.
%
% Note:
% N/A.
%
% Example:
% N/A. Set breakpoints in program.
%
% Children function:
% GLYVIS
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2016, Research Foundation for State University of New York. All rights reserved
%

msposout = {plotinfoout.mspos};
bondmapout = {plotinfoout.bondmap};
alllinkout = {plotinfoout.alllinkout};
allmonosacout = {plotinfoout.allms};
directionseqout = {plotinfoout.directionseq};
if isfield(plotinfoout,'CHAR')
    identity = {plotinfoout.CHAR};
else
    identity = cell(size(plotinfoout));
end
if isfield(plotinfoout,'bonddescription')
    bondescription = {plotinfoout.bonddescription};
else
    bondescription = cell(size(plotinfoout));
end

if strcmpi(options.workingmode,'g')
    switch lower(options.orientation)
        case 'up'
            previousmspos = -[max(msposout{1}(:,1)) + options.structspacing,0];
        case 'down'
            previousmspos = -[max(msposout{1}(:,1)) + options.structspacing,0];
        case 'left'
            previousmspos = -[0,max(msposout{1}(:,2)) + options.structspacing];
        case 'right'
            previousmspos = -[0,max(msposout{1}(:,2)) + options.structspacing];
    end
%     previousmspos = [0,max(msposout{1}(:,2)) + options.structspacing];
    availableMSlist = properties(MonoColor);
    for i = 1:numel(msposout)
        mscolor = cell(size(allmonosacout{i}));
        msshape = mscolor;
        %% special case: deoxyhexose 6dAlt & 6dTal
        for j = 1:numel(allmonosacout{i})
            if ismember(allmonosacout{i}{j},availableMSlist)
                mscolor{j} = MonoColor.(allmonosacout{i}{j});
                msshape{j} = MonoShape.(allmonosacout{i}{j});
            elseif strcmp(allmonosacout{i}{j},'6dTal')
                mscolor{j} = 'Light Blue';
                msshape{j} = 'Filled Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dAlt')
                mscolor{j} = 'Pink';
                msshape{j} = 'Filled Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dGul')
                mscolor{j} = 'Orange';
                msshape{j} = 'Filled Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dAltNAc')
                mscolor{j} = 'Pink';
                msshape{j} = 'Divided Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dTalNAc')
                mscolor{j} = 'Light Blue';
                msshape{j} = 'Divided Triangle';
            elseif strcmp(allmonosacout{i}{j},'4eLeg')
                mscolor{j} = 'Light Blue';
                msshape{j} = 'Flat Diamond';
            else
                mscolor{j} = 'White';
                msshape{j} = 'Flat Hexagon';
            end
        end
        mspara.shape = msshape;
        mspara.color = mscolor;
        mspara.size = ones(size(mscolor)) * options.monosacsize;
        mspara.periwidth = ones(size(mscolor)) * options.msperiwidth;
        mspara.name = allmonosacout{i};
        switch lower(options.orientation)
            case 'up'
                previousmspos = previousmspos + [max(msposout{i}(:,1)) + options.structspacing,0];
            case 'down'
                previousmspos = previousmspos + [max(msposout{i}(:,1)) + options.structspacing,0];
            case 'left'
                previousmspos = previousmspos + [0,max(msposout{i}(:,2)) + options.structspacing];
            case 'right'
                previousmspos = previousmspos + [0,max(msposout{i}(:,2)) + options.structspacing];
        end
        mspos = msposout{i} + repmat(previousmspos,size(msposout{i},1),1);
        glyvis(mspos,bondmapout{i},alllinkout{i},bondescription{i},mspara,directionseqout{i},identity{i},options)
        switch lower(options.orientation)
            case 'up'
                previousmspos = [max(mspos(:,1)),0];
            case 'down'
               previousmspos = [max(mspos(:,1)),0];
            case 'left'
                previousmspos = [0,max(mspos(:,2))];
            case 'right'
                previousmspos = [0,max(mspos(:,2))];
        end
    end
elseif strcmpi(options.workingmode,'gp')
    availableMSlist = properties(MonoColor);
    for i = 1:numel(msposout)
        mscolor = cell(size(allmonosacout{i}));
        msshape = mscolor;
        for j = 1:numel(allmonosacout{i})
            if ismember(allmonosacout{i}{j},availableMSlist)
                mscolor{j} = MonoColor.(allmonosacout{i}{j});
                msshape{j} = MonoShape.(allmonosacout{i}{j});
            elseif strcmp(allmonosacout{i}{j},'6dTal')
                mscolor{j} = 'Light Blue';
                msshape{j} = 'Filled Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dAlt')
                mscolor{j} = 'Pink';
                msshape{j} = 'Filled Triangle';
            else
                mscolor{j} = 'White';
                msshape{j} = 'Flat Hexagon';
            end
        end
        mspara.shape = msshape;
        mspara.color = mscolor;
        mspara.size = ones(size(mscolor)) * options.monosacsize;
        mspara.periwidth = ones(size(mscolor)) * options.msperiwidth;
        mspara.name = allmonosacout{i};
        mspos = msposout{i};
        glyvis(mspos,bondmapout{i},alllinkout{i},bondescription{i},mspara,directionseqout{i},identity{i},options)
    end
end
end
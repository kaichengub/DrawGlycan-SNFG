function drawglycan(allsgp,varargin)
% DRAWGLYCAN: draw glycan or glycopeptide. This is the main program
%
% Syntax:
% drawglycan(allsgp,options)
% drawglycan(allsgp,optionname1,optionvalue1,...)
%
% Input:
% allsgp: string, glycan and glycopeptide string in IUPAC format.
% options: struct, user defined options listed below in Table format.
% optionname and optionvalue: user defined option name and value pairs.
%
% Output:
% N/A
%
% Available options (case sensitive):
% ________________________________________________________________________________________
% |  option name              purpose                   option value                     |
% |--------------------------------------------------------------------------------------|
% |  orientation           Orientation of the        String:'up','down','left','right'.  |
% |                        glycan structure          Default is 'right'.                 |
% |  monosacsize           Monosac. size             Number in range: [0,1].             |
%                                                    Default=0.5.                        |
% |  msperiwidth           Monosac. border           Any integer. Default=1.             |
% |                        thickness                                                     |
% |  bondwidth             Molecular bond            Any integer. Default=2.             |
% |                        thickness                                                     |
% |  perpendicularmonosac  List of perpendicular     1 x n cell array of strings.        |
% |                        monosacs.                 Default= {'Fuc','Xyl'}.             |
% |  showlink              Diplay linkage info.      String, 'yes' or 'no'.              |
% |                                                  Default = 'yes'.                    |
% |  fontsize              Linkage data font size.   Any integer. Default = 12.          |
% |  linkinfodist          Linkage text position.    Any number. Default = 0.7.          |
% |  linkinfotheta         Text orientation with     Any number in range = [-90,90].     |
% |                        respect to linkage.       Default = 30 (degrees).             |
% |  bondbreaksiglength    Length of glycosidic      Any number. Default = 0.5.          |
% |                        bond fragmentation line.                                      |
% |  structspacing         Spacing between glycans   Any number. Default = 1.            |
% |                        listed in cell array.                                         |
% |  aaspacing             Spacing between amino     Any number. Default = 0.75.         |
% |                        acids.                                                        |
% |  fileout               File saving specs.        String, output figure file          |
% |                                                  name. Default = empty.              |
% |  visible               Display figure.           String,'on' or 'off',Default='on'.  |
% ----------------------------------------------------------------------------------------
% Example:
% drawglycan('A[GalNAc(a1-3)](-N "B2" -C  "X4")AB[Xyl(b1-3)GalNAc(a1-3)]B(-C "X2.5")C[Fuc(a1-3)]CDD')
%
% Children function:
% USG, CALCGLYPOS, ESTIFIGSIZE, PAINTGLYCAN, GETFRAGINFO, REARRANGEGLYPEP,
% PLOTPEPBREAK
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2016, Research Foundation for State University of New York. All rights reserved
%
% v1.1 update note: draw one structure at a time, the function of drawing
% multiple structure from a cell input has been cancelled.
defoptname = {'orientation','monosacsize','msperiwidth',...
    'bondwidth','perpendicularmonosac',...
    'linkinfodist','linkinfotheta','bondbreaksiglength',...
    'showlink','fontsize','workingmode',...
    'structspacing','aaspacing','fileout',...
    'visible','inputformat','displaybrackettext'};

defoptvalue = {'left',.5,1,...
    2,{'Fuc','Xyl'},...
    .7,30,.5,...
    'yes',12,'default',...
    1,.75,'',...
    'on','iupac','no'};

defopt = usg('gen',defoptname,defoptvalue);
if ~isempty(varargin)
    if isstruct(varargin{1})
        options = usg('mod',defopt,varargin{1});
    elseif ischar(varargin{1})
        options = usg('mod',defopt,varargin(1:2:end),varargin(2:2:end));
    end
else
    options = defopt;
end

if ismember(upper(options.orientation),{'UP','RIGHT'})
    options.linkinfotheta = abs(options.linkinfotheta);
elseif ismember(upper(options.orientation),{'DOWN','LEFT'})
    options.linkinfotheta = -abs(options.linkinfotheta);
end
switch lower(options.inputformat)
    case 'glycam'
        if ~any(strfind(allsgp,'('))
            allsgp = glycam2iupac(allsgp);  % not tested
        end
        options.inputformat = 'iupac';
end
if strcmpi(options.inputformat,'iupac')
    if strcmpi(options.workingmode,'DEFAULT') || strcmpi(options.workingmode,'GP')
        [glyseq,pepseq,glypos] = distggp(allsgp);  % glyseq is cell, pepseq is char
        if isempty(pepseq)
            options.workingmode = 'g';
        else
            options.workingmode = 'gp';
        end
    elseif strcmpi(options.workingmode,'G')
        glyseq = allsgp;
        pepseq = '';
        glypos = [];
    end
end
switch upper(options.workingmode)
    case 'G'
        if strcmpi(options.inputformat,'iupac')
            glyseq = strtrim(glyseq);
            [cleanseq,specialoptions] = getglydressing(glyseq);
            glysgp = IUPAC2Sgp(cleanseq,2);
        elseif strcmpi(options.inputformat,'sgp')
            glysgp = allsgp;
            %% need a module to recognize sgp based special options
            specialoptions = [];
        end
        if ~isempty(specialoptions)  % get option position in sgp seq, because it needs conversion
            spoptfld = fieldnames(specialoptions);  % special option fields
            numMono = length(strfind(glysgp,'{'));
            for i = 1:length(spoptfld)
                tempspopt = specialoptions.(spoptfld{i});
                for j = 1:size(tempspopt,1)
                    tempspopt{j,2} = bsxfun(@minus,numMono + 1,tempspopt{j,2});  % convert option associated monosac location from IUPAC to SGP
                end
                specialoptions.(spoptfld{i}) = tempspopt;
            end
        end
        plotinfoout = calcglypos(glysgp,options,specialoptions);  % only draw 1 glycan structure, for multiple structures, call func. multiple times
        % specialoptions are merged into plotinfoout
        if ~isempty(specialoptions)
            plotinfoout = calcattachmentpos(plotinfoout,options);
        end
        figure, hold on, set(gcf,'visible',options.visible)
        paintglycan(plotinfoout,options)
        paintattachment(plotinfoout,options)
        estimatefigsize(plotinfoout,'',options);
    case 'GP'  % default format: IUPAC
        options.orientation = 'up';
        glyseq = cellfun(@strtrim,glyseq,'uniformoutput',false);
        [cleanseq,specialoptions] = cellfun(@getglydressing,glyseq,'uniformoutput',false);
        glysgp = cellfun(@IUPAC2Sgp,cleanseq,num2cell(repmat(2,numel(glyseq),1)),'uniformoutput',false);
        plotinfoout = cell(size(glysgp));
        for i = 1:length(glysgp)
            if ~isempty(specialoptions{i})  % get option position in sgp seq, because it needs conversion
                spoptfld = fieldnames(specialoptions{i});  % special option fields
                numMono = length(strfind(glysgp{i},'{'));
                for j = 1:length(spoptfld)
                    tempspopt = specialoptions{i}.(spoptfld{j});
                    for k = 1:size(tempspopt,1)
                        tempspopt{k,2} = bsxfun(@minus,numMono + 1,tempspopt{k,2});  % convert option associated monosac location from IUPAC to SGP
                    end
                    specialoptions{i}.(spoptfld{j}) = tempspopt;
                end
            end
            plotinfoout{i} = calcglypos(glysgp{i},options,specialoptions{i});  % only draw 1 glycan structure, for multiple structures, call func. multiple times
            % specialoptions are merged into plotinfoout
            if ~isempty(specialoptions{i})
                plotinfoout{i} = calcattachmentpos(plotinfoout{i},options);
            end
        end
        [pepbreak,cleanpepseq] = getseqoptions(pepseq,'p');
        [plotinfoout,peppos] = rearrangeglypep(cleanpepseq,glypos,plotinfoout,options);
        figure, hold on, set(gcf,'visible',options.visible)
        for i = 1:length(cleanpepseq)
            text(peppos(i),-1,char(double(cleanpepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
                'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
        end
        for i = 1:length(plotinfoout)
            paintglycan(plotinfoout{i},options)
            paintattachment(plotinfoout{i},options)
        end
        estimatefigsize(plotinfoout,peppos,options)
    case 'P'
        
end
set(gca,'position',[0,0,1,1],'visible','off')
axis equal
if strcmpi(options.workingmode,'gp')
    plotpepbreak(peppos,pepbreak,glypos,options)  % this is the last step because frag marker length is related to figure size
end

if ~isempty(options.fileout)
    saveas(gcf,options.fileout)
end
if strcmpi(options.visible,'off')
    close gcf
end
end
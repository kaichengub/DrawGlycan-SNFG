function estimatefigsize(plotinfoout,peppos,options)
% ESTIMATEFIGSIZE: Designate figure size according to the size of glycans
% and peptide backbone.
%
% Syntax:
% estimatefigsize(plotinfoout,peppos,options)
%
% Input:
% plotinfoout: m x n cell array, contains the plot information used in
% painting program, 1st column contains the positioning info for each
% monosaccharide in each glycan.
% peppos: 1 x n numerical array, the position of each amino acid in the
% peptide backbone. The numbers are the x position of amino acid. They will
% be put at [x,-1] with "verticalalignment" set to "top".
% options: struct, options.
%
%
% Output:
% horirec: coordinates for the left and right boundaries of the figure
% maxvert: figure height
% Figure size will be adjusted based on the estimation program
%
% Note:
% Figure size estimation is based on 70 pixels per unit. Program will
% read proposed glycan structure size and peptide backbone
% length, then reset figure size.
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

% xlim = [0,0];
% ylim = [0,0];
% if isempty(peppos)
%     mode = 'g';
% else
%     mode = options.workingmode;
% end
pixelsperunit = DrawGlycanPara.pixelsperunit;
currentfig = gcf;
switch lower(options.workingmode)
    case 'g'
        mainmspos = plotinfoout.mspos;
        if isfield(plotinfoout,'CBITEM')
            cbitemmspos = plotinfoout.CBITEM.cbplotinfosto;
        else
            cbitemmspos = [];
        end
        if isfield(plotinfoout,'adducttextpos')
            adducttextpos = plotinfoout.adducttextpos;
        else
            adducttextpos = [];
        end
        letterwidth = options.fontsize*1.9 * 20/8/pixelsperunit/12;
        xlim = [min(mainmspos(:,1)),max(mainmspos(:,1))];  % initial guess
        ylim = [min(mainmspos(:,2)),max(mainmspos(:,2))];
        switch lower(options.orientation)
            case 'up'
                if ~isempty(cbitemmspos)
                    structwidth = 0;
                    structheight = 0;
                    for i = 1:length(cbitemmspos)
                        structwidth = structwidth + max(cbitemmspos(i).mspos(:,1)) - min(cbitemmspos(i).mspos(:,1)) + 1;
                        structheight = max([structheight,max(cbitemmspos(i).mspos(:,2)-cbitemmspos(i).mspos(1,2))]);
                    end
                    cbxlim = [min(cbitemmspos(1).mspos(:,1)),min(cbitemmspos(1).mspos(:,1)) + structwidth - 1];
                    ylim(2) = ylim(2) + structheight + 1;
                    xlim = [min([xlim(1),cbxlim(1)]),max([xlim(2),cbxlim(2)])];
                end
                
                
            case 'down'
                if ~isempty(cbitemmspos)
                    structwidth = 0;
                    structheight = 0;
                    for i = 1:length(cbitemmspos)
                        structwidth = structwidth + max(cbitemmspos(i).mspos(:,1)) - min(cbitemmspos(i).mspos(:,1)) + 1;
                        structheight = max([structheight,-min(cbitemmspos(i).mspos(:,2)-cbitemmspos(i).mspos(1,2))]);
                    end
                    cbxlim = [min(cbitemmspos(1).mspos(:,1)),min(cbitemmspos(1).mspos(:,1)) + structwidth - 1];
                    ylim(1) = ylim(1) - structheight - 1;
                    xlim = [min([xlim(1),cbxlim(1)]),max([xlim(2),cbxlim(2)])];
                end
                
                
            case 'left'
                if ~isempty(cbitemmspos)
                    structwidth = 0;
                    structheight = 0;
                    for i = 1:length(cbitemmspos)
                        structwidth = max([structwidth,-min(cbitemmspos(i).mspos(:,1)-cbitemmspos(i).mspos(1,1))]);
                        structheight = structheight + max(cbitemmspos(i).mspos(:,2)) - min(cbitemmspos(i).mspos(:,2)) + 1;
                    end
                    cbylim = [min(cbitemmspos(1).mspos(:,2)),min(cbitemmspos(1).mspos(:,2)) + structheight - 1];
                    xlim(1) = xlim(1) - structwidth - 1;
                    ylim = [min([ylim(1),cbylim(1)]),max([ylim(2),cbylim(2)])];
                end
                
            case 'right'
                if ~isempty(cbitemmspos)
                    structwidth = 0;
                    structheight = 0;
                    for i = 1:length(cbitemmspos)
                        structwidth = max([structwidth,max(cbitemmspos(i).mspos(:,1)-cbitemmspos(i).mspos(1,1))]);
                        structheight = structheight + max(cbitemmspos(i).mspos(:,2)) - min(cbitemmspos(i).mspos(:,2)) + 1;
                    end
                    cbylim = [min(cbitemmspos(1).mspos(:,2)),min(cbitemmspos(1).mspos(:,2)) + structheight - 1];
                    xlim(2) = xlim(2) + structwidth + 1;
                    ylim = [min([ylim(1),cbylim(1)]),max([ylim(2),cbylim(2)])];
                end
                
        end
        
        if ~isempty(adducttextpos)
            adduct = plotinfoout.ADDUCT;
            adducttext = '';
            for i = 1:size(adduct,1)
                adducttext = [adducttext,adducttext,' '];
            end
            %                 adducttext = adducttext(1:end-1);
            xlim(2) = adducttextpos(1) + (length(adducttext) - 1) * letterwidth;
        end
        set(gca,'XLim',xlim + [-1,1],'YLim',ylim + [-1,1]);
        set(currentfig,'position',[currentfig.Position(1:2)/4,(diff(xlim) + 2)*pixelsperunit,(diff(ylim)+2)*pixelsperunit]);
        
        
    case 'gp'
        LR = [min(peppos),max(peppos)];
        maxvert = 0;
        for i = 1:length(plotinfoout)
            thismspos = plotinfoout{i}.mspos;
            maxvert = max([maxvert,max(thismspos(:,2))]);
            thismsLR = [min(thismspos(:,1)),max(thismspos(:,1))];
            %         thismsright = min(thismspos(:,1));
            fldnames = lower(fieldnames(plotinfoout{i}));
            if ismember('adducttextpos',fldnames)
                adduct = plotinfoout{i}.ADDUCT;
                adducttext = '';
                for j = 1:size(adduct,1)
                    adducttext = [adducttext,adduct{j,1},' '];
                end
                adductLR = [min(plotinfoout{i}.bracketspos(:,1)),plotinfoout{1}.adducttextpos(1) + (length(adducttext)-1) * DrawGlycanPara.letterwidthpxl/DrawGlycanPara.pixelsperunit];
                %             adductright = plotinfoout{1}.adducttextpos(1) + (length(adducttext)-1) * 22.85/10;
            else
                adductLR = [];
            end
            if ismember('cbitem',fldnames)
                cbLR = [0,0];
            cbplotinfos = plotinfoout{i}.CBITEM.cbplotinfosto;
            for j = 1:size(cbplotinfos,1)
                cbitempos = plotinfoout{i}.CBITEM.cbplotinfosto(j).mspos;
                cbLR = [min([cbLR(1),min(cbitempos(:,1))]),max([cbLR(2),max(cbitempos(:,1))])];
%                 plotinfoout{i}.CBITEM.cbplotinfosto(j).mspos = plotinfoout{i}.CBITEM.cbplotinfosto(j).mspos + repmat([rootpos(i),0],size(plotinfoout{i}.CBITEM.cbplotinfosto(j).mspos,1),1);
            maxvert = max([maxvert,max(cbitempos(:,2))]);
            end
                
            else
                cbLR = [];
            end
            tempallLR = [thismsLR;adductLR;cbLR];
            LR = [min([LR(1),min(tempallLR(:,1))]),max([LR(2),max(tempallLR(:,2))])];
        end

%         horirec = [min(plotinfoout{1,1}(:,1) - plotinfoout{1,1}(1,1)),min(plotinfoout{1,1}(:,1) - plotinfoout{1,1}(1,1))];
        % horirec: 1st element records the coordinate of the glycan to the
        % furthest left, 2nd records the furthest right
%         maxvert = 0;
%         for i = 1:size(plotinfoout,1)
%             horirec = horirec + [0,max(plotinfoout{i,1}(:,1)) - min(plotinfoout{i,1}(:,1)) + options.structspacing];
%             % It is possible that in REARRANGEGLYPEP the monosac position was
%             % changed, need to compensate for that
%             if max(plotinfoout{i,1}(:,2)) >= maxvert
%                 maxvert = max(plotinfoout{i,1}(:,2));
%             end
%         end
        if maxvert == 0
            maxvert = 1;  % compensation for glycans with max. dist. = 0, one monosac or one monosac with fucose attached.
        end
        endcomp = diff(peppos);  % compensate space for both ends of peptide backbone
        if ~isempty(endcomp)  % need to consider this only when num. of AA > 1
            if LR(1) > min(peppos) - endcomp(1)
                LR(1) = min(peppos) - endcomp(1);  % if peptide backbone is long
            end
            if LR(2) < max(peppos) + endcomp(end)
                LR(2) = max(peppos) + endcomp(end);
            end
        end
        set(gca,'XLim',LR+[-1,1])
        set(gca,'YLim',[-2,maxvert+1])
        set(gcf,'position',[currentfig.Position(1:2)/4,(sum(abs(LR))+2)*pixelsperunit,(maxvert+3)*pixelsperunit])
end
end

function plotpepbreak(peppos,pepbreak,glycosylpos,options)
% PLOTPEPBREAK: plot peptide backbone fragmentation in glycopeptide figure.
% 
% Syntax:
% plotpepbreak(peppos,pepbreak,glycosylpos,options)
% 
% Input:
% peppos: 1 x n numerical array, the position of the amino acids in the
% backbone.
% pepbreak: 1 x 3 cell array, 1st element contains m x 1 numerical array,
% the position of the fragmentation mark, 2nd contains m x 1 cell array of
% strings, the type of the marks, 3rd contains m x 1 cell array of strings,
% the content of the marks.
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


if isempty(pepbreak)
    for i = 1:length(glycosylpos)
        plot([peppos(glycosylpos(i)),peppos(glycosylpos(i))],[-1,-options.monosacsize/2],'k','linewidth',options.bondwidth)
    end
else
    pepbreaktyp = fieldnames(pepbreak);
    pepbreakpos = [];
    pepbreakcont = {};
    pepbreaktype = {};
    for i = 1:length(pepbreaktyp)
        if ismember(pepbreaktyp{i},{'N','C'})
            pepbreakinfo = pepbreak.(pepbreaktyp{i});
            pepbreakpos = [pepbreakpos;cell2mat(pepbreakinfo(:,2))];
            pepbreaktype = [pepbreaktype;repmat(pepbreaktyp(i),size(pepbreakinfo,1),1)];
            pepbreakcont = [pepbreakcont;pepbreakinfo(:,1)];
        end
    end
%     pepbreakpos = pepbreak{1};
    [unibreakpos,~,unibreakposind] = unique(pepbreakpos);
%     pepbreaktyp = pepbreak{2};
%     pepbreakcont = pepbreak{3};
    currentfig = gcf;
    height = currentfig.Position(4);
    currentax = gca;
    ylim = currentax.YLim;
    pixelperunit = height/diff(ylim);
    letterheight = options.fontsize*1.5 * 23/pixelperunit/12;
    if any(ismember('N',pepbreaktyp))
        for i = 1:length(glycosylpos)
            plot([peppos(glycosylpos(i)),peppos(glycosylpos(i))],[-1+letterheight/(1.5/12*options.fontsize),-options.monosacsize/2],'k','linewidth',options.bondwidth)
        end  % leave room for peptide fragmentation text by shortening bond length connecting AA and monosac.
    elseif any(ismember('C',pepbreaktyp))
        for i = 1:length(glycosylpos)
            plot([peppos(glycosylpos(i)),peppos(glycosylpos(i))],[-1,-options.monosacsize/2],'k','linewidth',options.bondwidth)
        end
    end
    for i = 1:length(unibreakpos)
        if i == length(peppos)
            posx = peppos(i) + options.aaspacing/2;
        else
            posx = mean([peppos(unibreakpos(i)),peppos(unibreakpos(i)+1)]);
        end
        posy = [-1,-1-letterheight];
        plot([posx,posx],posy,'k')
        localpepbreaktyp = pepbreaktype(unibreakposind == i);
        localpepbreakcont = pepbreakcont(unibreakposind == i);
        for j = 1:numel(localpepbreaktyp)
            if strcmpi(localpepbreaktyp{j},'N')
                posx2 = posx - options.bondbreaksiglength/3;
                plot([posx,posx2],[posy(2),posy(2)],'k')
                text(posx,posy(2),char(double(localpepbreakcont{j})),'fontsize',options.fontsize,'VerticalAlignment','top','HorizontalAlignment','right')
            elseif strcmpi(localpepbreaktyp{j},'C')
                posx2 = posx + options.bondbreaksiglength/3;
                plot([posx,posx2],[posy(1),posy(1)],'k')
                text(posx,posy(1),char(double(localpepbreakcont{j})),'fontsize',options.fontsize,'VerticalAlignment','bottom')
            end
        end
    end
end
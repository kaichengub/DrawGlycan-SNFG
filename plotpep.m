function plotpep(pepseq,peppos,options)
% PLOTPEP: display peptide backbone at specified position
%
% Syntax:
% plotpep(pepseq,peppos,options)
% 
% Input:
% pepseq: string, the amino acid sequence of the peptide backbone.
% peppos: 1 x n numerical array, the x-value of the position of the amino
% acids.
% options: structure, display options.
% 
% Output:
% N/A
% 
% Note:
% Peptides are plotted at y = -1, aligned at the top.
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


for i = 1:length(peppos)
    text(peppos(i),-1,char(double(pepseq(i))),'FontName','Helvetica','VerticalAlignment','top',...
        'HorizontalAlignment','center','fontsize',options.fontsize*1.5);
end
end
function [gly,pep,glypos] = distggp(IUPACinput)
% DISTGGP: Parse the input IUPAC sequence to provide glycan and peptide
% string, identify the input string, distinguishing between glycan and
% glycopeptide.
%
% Syntax:
% [gly,pep,glypos] = distggp(IUPACinput)
%
% Input:
% IUPACinput: string, the IUPAC sequence of the
% glycan/glycopeptide to be identified.
%
% Output:
% gly: n x 1 cell array of strings, the glycan sequence of
% glycan/glycopeptide.
% pep: string, the peptide sequence of glycopeptide, this output is blank
% when inputs are glycan only.
% glypos: n x 1 matrix, position where glycan is attached to peptide backbone, 
% counting from N-terminus.
%
% Note:
% Input must be in a string.
% Input must be consists of pure glycan(s) or contains only one
% glycopeptide. Mixture of glycan and glycopeptide is prohibited.
%
% Example:
% [gly,pep,glypos] = distggp('FKT[(??-?)GalNAc[(??-?)Fuc](??-?)Gal]GTK')
% 
% gly = 
% 
%     '(??-?)GalNAc[(??-?)Fuc](??-?)Gal'
% 
% 
% pep =
% 
% FKTGTK
% 
% 
% glypos =
% 
%      3
% 
%
% [gly,pep,glypos] = distggp('(??-?)GalNAc[(??-?)Fuc](??-?)Gal')
% 
% gly =
% 
% (??-?)GalNAc[(??-?)Fuc](??-?)Gal
% 
% 
% pep =
% 
%      ''
% 
% 
% 
% glypos =
% 
%      0
% 
% 
% Children function:
% N/A

% ver 1.1, added output "glypos"

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2017, Research Foundation for State University of New York. All rights reserved
%
thissgp = IUPACinput;
thissgp = strtrim(thissgp);
[comstart,comend] = regexp(thissgp,DrawGlycanPara.regexp_optionvalue,'start','end');
tempthissgp = thissgp;
%% block user-customized info using whitespaces
for i = 1:length(comstart)
    tempthissgp(comstart(i):comend(i)) = ' ';
end
%% identify glycan from glypeptide by counting brackets
opensqbr = strfind(tempthissgp,'[');
closesqbr = strfind(tempthissgp,']');
[sqbrpos,ind] = sort([opensqbr closesqbr]);
sqbrcounter = [ones(size(opensqbr)),-ones(size(closesqbr))];
sqbrcounter = sqbrcounter(ind);
sqbrcounterind = arrayfun(@(x) sum(sqbrcounter(1:x)),1:length(sqbrcounter));
sqbrecorder = [0 sqbrcounterind] == 0;
sqbrecorder = [sqbrpos(sqbrecorder(1:end-1));sqbrpos(sqbrcounterind == 0)];
glypos = zeros(size(sqbrecorder,2),1);
%% memorizing glycan position
pepseq = thissgp;  % now we don't know if it's a pep or gly, call it pepseq for now
glyseq = cell(size(sqbrecorder,2),1);
for i = size(sqbrecorder,2):-1:1
    pepseq(sqbrecorder(1,i):sqbrecorder(2,i)) = '';
    glyseq{i} = thissgp(sqbrecorder(1,i)+1:sqbrecorder(2,i)-1);
end
%% peptide options: e.g. fragmentation
[pepoptstart,pepoptend] = regexp(pepseq,'\((.*?)\)','start','end');
pepseq2 = pepseq;  % Need this duplicate for removing lowercase letters in pep. fragmentation info.
for i = length(pepoptstart):-1:1
    pepseq2(pepoptstart(i):pepoptend(i)) = '';
end
if any(isstrprop(pepseq2,'lower'))  % There is at least one lowercase letter, this is a glycan. Correct input is very important for proper functioning
    gly = IUPACinput;
    pep = '';
else  % input is glycopeptide
    gly = glyseq;
    pep = strtrim(pepseq);
end
if ~isempty(pep)
    %% Get glycosylation position after eliminating peptide fragment info
    glyseqcompen = 0;
    for i = 1:size(sqbrecorder,2)
        glypos(i) = sqbrecorder(1,i) - 1 - glyseqcompen;
        for j = length(pepoptstart):-1:1
            if any(glypos >= pepoptend(j))
                glypos(glypos >= pepoptend(j)) = glypos(glypos >= pepoptend(j)) - (pepoptend(j) - pepoptstart(j) + 1);
            end
        end
        glyseqcompen = glyseqcompen + sqbrecorder(2,i) - sqbrecorder(1,i) + 1;
    end
end
end
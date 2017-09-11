classdef DrawGlycanPara
    
    properties (Constant)
        %% regular expression
        regexp_monosac_SGP = '(?<={)[^{}]+[{}]*?';
        regexp_monosac_IUPAC = '(?<=[\[\]\)])[^\[\]]+?\)';
        regexp_monosaclinkage = '\((.*?)\)';
        regexp_option = '\-[A-z]+';
        regexp_optionvalue = '["''�].*?["''�]';
        regexp_cbitem = '{.+?}';
        regexp_MonoPattern_Split = '[A-Z_a-z_0-9]*\((.*?)\)'
        regexp_BracketPattern = '[\[\]]';
        regexp_AA = '[A-Z][a-z]*';
        %% variables by categories
        intglymodinfo = {'U','D','P'};
        intglybondinfo = {'R','NR','BOLD','ZIG','DASH','WAVY','WEDGE'};
        extglymodinfo = {'RS','RE','ADDUCT','ANOMER','CURLY'};
        glyidentityinfo = {'CHAR'};
        aamodinfo = {'N','C'};
        %% constants
        tipradius = 0.25;
        pixelsperunit = 70;
        letterwidthpxl = 22.85/10;
        letterheightpxl = 22.85/24;
        
    end
end
function newLinkInfo = linkrepair(LinkInfo)
LinkInfo = LinkInfo{1}(2:end-1);
if isempty(LinkInfo)  % no link/frag info
    newLinkInfo = '(??-?';
    fraginfo = '';
else  % have link/frag info
    fragstart = regexp(LinkInfo,'-R|-NR|-U|-D','start');
    if isempty(fragstart)
        fraginfo = '';
        fragstart = length(LinkInfo) + 1;
    else
        fraginfo = [' ',LinkInfo(fragstart:end)];
    end
    glybond = LinkInfo(1:min(fragstart)-1);
    glybond = strrep(glybond,',','-');
    if isempty(glybond)
        newLinkInfo = '(??-?';
    elseif length(glybond) == 1  % bond info is question mark/1 letter/1 digit
        glybondascii = double(LinkInfo(1:min(fragstart)-1));
        if glybondascii == 63 %% question mark
            newLinkInfo = '(??-?';
        elseif glybondascii >47 && glybondascii < 58  % 0-9
            newLinkInfo = ['(??-',LinkInfo(1:min(fragstart)-1)];
        elseif (glybondascii >64 && glybondascii < 91) || (glybondascii >96 && glybondascii < 123)  % a-z and A-Z
            newLinkInfo = ['(',LinkInfo(1:min(fragstart)-1),'?-?'];
        else  % unrecognizable
            newLinkInfo = '(??-?';
        end
    else  % bond info is longer
        carbon = regexp(glybond,'[a-z]','once');
        qmark = regexp(glybond,'?','once');
        pairnum = regexp(glybond,'[\d/\\]+\-[?\d/\\]+','once','match');
        [anynum,anynumstart] = regexp(glybond,'[\d/\\]+','once','match','start');
        if isempty(carbon)  % long bond info, no letter, question mark/number only
            if isempty(qmark)  % long bond info, no letter, no question mark
                if isempty(pairnum)  %  long bond info, no letter, no question mark, number not in pair
                    newLinkInfo = ['(??-',glybond];
                else  %  long bond info, no letter, no question mark, number in pair
                    newLinkInfo = ['(?',glybond];
                end
            elseif strcmp(glybond,'??-?')
                newLinkInfo = '(??-?';
            elseif ~isempty(anynum) && ~isempty(qmark)  % long bond info, no letter, question mark/number both exist
                if isempty(pairnum)  % long bond info, no letter, question mark exist, number not in pair
                    if qmark < anynumstart  % (?#)
                        newLinkInfo = ['(??-',anynum];
                    else  % (#?)
                        newLinkInfo = ['(?',anynum,'-?'];
                    end
                else % long bond info, no letter, question mark exist, number in pair
                    newLinkInfo = ['(?',pairnum];
                end
            else  % long bond info, unrecognizable
                newLinkInfo = '(??-?';
            end
        else  % long bond info, have carbon letter
            if isempty(qmark)  % long bond info, have carbon letter, no question mark
                if isempty(pairnum)  % long bond info, have carbon letter, no question mark, number not in pair
                    newLinkInfo = ['(',glybond(carbon),'?-',anynum];
                else  % long bond info, have carbon letter, no question mark, number in pair
                    newLinkInfo = ['(',glybond(carbon),pairnum];
                end
            else  % long bond info, have carbon letter, have question mark
                newLinkInfo = ['(',glybond(carbon),'?-?'];
            end
        end
    end
end
newLinkInfo = [newLinkInfo,fraginfo];
end
function paintattachment(plotinfoout,options)
%% adduct, anomer curlybracket,repeat
mspos = plotinfoout.mspos;
bondmap = plotinfoout.bondmap;
tr = DrawGlycanPara.tipradius;


%% 1.cbitem
%% analysis of curly bracet add-on structure
if isfield(plotinfoout,'CBITEM')
    cbitem = plotinfoout.CBITEM;
    waypoints = cbitem.bracketwaypoints;
    plot(waypoints(2:3,1),waypoints(2:3,2),'k')
    plot(waypoints(6:7,1),waypoints(6:7,2),'k')
    mainmspos = plotinfoout.mspos;
    switch lower(options.orientation)
        case 'up'
            drawcurve(pi/2,pi,waypoints(1,1),waypoints(1,2),tr)  % I
            drawcurve(pi*3/2,pi*2,waypoints(4,1),waypoints(4,2),tr)  % II
            drawcurve(pi,pi*3/2,waypoints(5,1),waypoints(5,2),tr)  % III
            drawcurve(0,pi/2,waypoints(8,1),waypoints(8,2),tr)  % IV
        case 'down'
            drawcurve(pi,pi*3/2,waypoints(1,1),waypoints(1,2),tr)  % I
            drawcurve(0,pi/2,waypoints(4,1),waypoints(4,2),tr)  % II
            drawcurve(pi/2,pi,waypoints(5,1),waypoints(5,2),tr)  % III
            drawcurve(pi*3/2,pi*2,waypoints(8,1),waypoints(8,2),tr)  % IV
        case 'left'
            drawcurve(pi,pi*3/2,waypoints(1,1),waypoints(1,2),tr)  % I
            drawcurve(0,pi/2,waypoints(4,1),waypoints(4,2),tr)  % II
            drawcurve(pi*3/2,pi*2,waypoints(5,1),waypoints(5,2),tr)  % III
            drawcurve(pi/2,pi,waypoints(8,1),waypoints(8,2),tr)  % IV
        case 'right'
            drawcurve(pi*3/2,pi*2,waypoints(1,1),waypoints(1,2),tr)  % I
            drawcurve(pi/2,pi,waypoints(4,1),waypoints(4,2),tr)  % II
            drawcurve(pi,pi*3/2,waypoints(5,1),waypoints(5,2),tr)  % III
            drawcurve(0,pi/2,waypoints(8,1),waypoints(8,2),tr)  % IV
    end
    cboption = options;
    cboption.workingmode = 'g';
    paintglycan(cbitem.cbplotinfosto,cboption)
end

%% 2. repeat
if isfield(plotinfoout,'RS')
    repeatstart = plotinfoout.RS;
    repeatend = plotinfoout.RE;
    if strcmpi(options.orientation,'UP')
        textpos = 'end';
    elseif strcmpi(options.orientation,'DOWN')
        textpos = 'start';
        inter = repeatend;
        repeatend = repeatstart;
        repeatstart = inter;
    elseif strcmpi(options.orientation,'LEFT')
        textpos = 'start';
    elseif strcmpi(options.orientation,'RIGHT')
        textpos = 'end';
        inter = repeatend;
        repeatend = repeatstart;
        repeatstart = inter;
    end
    for i = 1:size(repeatend,1)
        thismonosacpos = mspos(repeatend{i,2},:);
        prevmonosacind = find(bondmap(:,repeatend{i,2}),1,'first');
        if ~isempty(prevmonosacind)
            prevmonosacpos = mspos(prevmonosacind,:);
        else  % only at the beginning of the structure, it's the nature of tree structure
            switch lower(options.orientation)
                    case 'up'
                        prevmonosacpos = thismonosacpos - [0,1];
                    case 'down'
                        prevmonosacpos = thismonosacpos + [0,1];
                    case 'left'
                        prevmonosacpos = thismonosacpos + [1,0];
                    case 'right'
                        prevmonosacpos = thismonosacpos - [1,0];
            end
        end
        if thismonosacpos(1)-prevmonosacpos(1) ~= 0
            alpha = atan((thismonosacpos(2)-prevmonosacpos(2))/(thismonosacpos(1)-prevmonosacpos(1)));
        else
            alpha = pi/2;
        end
        xdif = thismonosacpos(1)-prevmonosacpos(1);
        ydif = thismonosacpos(2)-prevmonosacpos(2);
        msdistance = sqrt(xdif^2+ydif^2);
        if xdif > 0 && ydif <= 0
            alpha = alpha + 2*pi;
        elseif xdif <= 0 && ydif < 0
            alpha = alpha + pi;
        elseif xdif < 0 && ydif >= 0
            alpha = alpha + pi;
        end  % fix alpha value
        centx = prevmonosacpos(1)*.3 + thismonosacpos(1)*.7;
        centy = prevmonosacpos(2)*.3 + thismonosacpos(2)*.7;
        linestart = [centx + options.bondbreaksiglength/1.6 * cos(alpha + pi/2),centy + options.bondbreaksiglength/1.6 * sin(alpha + pi/2)];
        lineend = [centx - options.bondbreaksiglength/1.6 * cos(alpha + pi/2),centy - options.bondbreaksiglength/1.6 * sin(alpha + pi/2)];
        plot([linestart(1),lineend(1)],[linestart(2),lineend(2)],'k','linewidth',options.bondwidth*.75);
        bracketlinefix = [options.bondbreaksiglength/ 4 * cos( 0 /180*pi+alpha),options.bondbreaksiglength/ 4 * sin( 0 /180*pi+alpha)];
        bracketlineend1 = linestart + bracketlinefix;
        bracketlineend2 = lineend + bracketlinefix;
        plot([bracketlineend1(1),linestart(1)],[bracketlineend1(2),linestart(2)],'k','linewidth',options.bondwidth*.75);
        plot([bracketlineend2(1),lineend(1)],[bracketlineend2(2),lineend(2)],'k','linewidth',options.bondwidth*.75);
        if ~isempty(repeatend{i,1})
            switch textpos
                case 'start'
                    brackettextpos = [linestart(1) + (linestart(1) - thismonosacpos(1))* 0.4,linestart(2) + (linestart(2) - thismonosacpos(2))* 0.4];
                    text(brackettextpos(1),brackettextpos(2),repeatend{i,1},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',options.fontsize * 1.25);
                case 'end'
                    brackettextpos = [lineend(1) + (lineend(1) - thismonosacpos(1))* 0.4,lineend(2) + (lineend(2) - thismonosacpos(2))* 0.4];
                    text(brackettextpos(1),brackettextpos(2),repeatend{i,1},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',options.fontsize * 1.25);
            end
        end
    end
    
    
    for i = 1:size(repeatstart,1)
        thismonosacpos = mspos(repeatstart{i,2},:);
        nextmonosacind = find(bondmap(repeatstart{i,2},:),1,'first');
        if ~isempty(nextmonosacind)
            nextmonosacpos = mspos(nextmonosacind,:);
        else  % end of antenna
            prevmonosacind = find(bondmap(:,repeatstart{i,2}),1,'first');
            if ~isempty(prevmonosacind)
                prevmonosacpos = mspos(prevmonosacind,:);
            else  % only at the beginning of the structure, it's the nature of tree structure
                switch lower(options.orientation)
                    case 'up'
                        prevmonosacpos = thismonosacpos - [0,1];
                    case 'down'
                        prevmonosacpos = thismonosacpos + [0,1];
                    case 'left'
                        prevmonosacpos = thismonosacpos + [1,0];
                    case 'right'
                        prevmonosacpos = thismonosacpos - [1,0];
                end
            end
            nextmonosacpos = 2*thismonosacpos - prevmonosacpos;
        end
        if thismonosacpos(1)-nextmonosacpos(1) ~= 0
            alpha = atan((nextmonosacpos(2)-thismonosacpos(2))/(nextmonosacpos(1)-thismonosacpos(1)));
        else
            alpha = pi/2;
        end
        xdif = nextmonosacpos(1) - thismonosacpos(1);
        ydif = nextmonosacpos(2) - thismonosacpos(2);
        msdistance = sqrt(xdif^2+ydif^2);
        if xdif > 0 && ydif <= 0
            alpha = alpha + 2*pi;
        elseif xdif <= 0 && ydif < 0
            alpha = alpha + pi;
        elseif xdif < 0 && ydif >= 0
            alpha = alpha + pi;
        end  % fix alpha value
        centx = nextmonosacpos(1)*.3 + thismonosacpos(1)*.7;
        centy = nextmonosacpos(2)*.3 + thismonosacpos(2)*.7;
        linestart = [centx + options.bondbreaksiglength/1.6 * cos(alpha + pi/2),centy + options.bondbreaksiglength/1.6 * sin(alpha + pi/2)];
        lineend = [centx - options.bondbreaksiglength/1.6 * cos(alpha + pi/2),centy - options.bondbreaksiglength/1.6 * sin(alpha + pi/2)];
        plot([linestart(1),lineend(1)],[linestart(2),lineend(2)],'k','linewidth',options.bondwidth*.75);
        bracketlinefix = [options.bondbreaksiglength/ 4 * cos( 0 /180*pi+alpha),options.bondbreaksiglength/ 4 * sin( 0 /180*pi+alpha)];
        brackettextfix = [options.bondbreaksiglength/ 4 * cos( 90 /180*pi+alpha),options.bondbreaksiglength/ 4 * sin( 90 /180*pi+alpha)];
        bracketlineend1 = linestart - bracketlinefix;
        bracketlineend2 = lineend - bracketlinefix;
        plot([bracketlineend1(1),linestart(1)],[bracketlineend1(2),linestart(2)],'k','linewidth',options.bondwidth*.75);
        plot([bracketlineend2(1),lineend(1)],[bracketlineend2(2),lineend(2)],'k','linewidth',options.bondwidth*.75);
        if ~isempty(repeatstart{i,1})
            switch textpos
                case 'start'
                    brackettextpos = [linestart(1) + (linestart(1) - thismonosacpos(1))* 0.4,linestart(2) + (linestart(2) - thismonosacpos(2))* 0.4];
                    text(brackettextpos(1),brackettextpos(2),repeatstart{i,1},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',options.fontsize * 1.25);
                case 'end'
                    brackettextpos = [lineend(1) + (lineend(1) - thismonosacpos(1))* 0.4,lineend(2) + (lineend(2) - thismonosacpos(2))* 0.4];
                    %                     text(lineend(1)-1.0*brackettextfix(1),lineend(2)-1.0*bracketlinefix(2),repeatend{i,1},'HorizontalAlignment','center');
                    text(brackettextpos(1),brackettextpos(2),repeatstart{i,1},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',options.fontsize * 1.25);
            end
        end
    end
    
end

%% 3.bracket
if isfield(plotinfoout,'bracketspos')
    bracketspos = plotinfoout.bracketspos;
    adduct = plotinfoout.ADDUCT;
    adducttextpos = plotinfoout.adducttextpos;
    plot(bracketspos(1:4,1),bracketspos(1:4,2),'k');
    plot(bracketspos(5:8,1),bracketspos(5:8,2),'k');
    adducttext = '';
    for i = 1:size(adduct,1)
        adducttext = [adducttext,adduct{i,1},' '];
    end
    adducttext = adducttext(1:end - 1);
    text(adducttextpos(1),adducttextpos(2),adducttext,'fontsize',options.fontsize * 1.5);
end
end

function drawcurve(a,b,h,k,r)
% Plot a circular arc as a pie wedge.
% a is start of arc in radians,
% b is end of arc in radians,
% (h,k) is the center of the circle.
% r is the radius.
% Try this:   plot_arc(pi/4,3*pi/4,9,-4,3)
% Author:  Matt Fig
t = linspace(a,b);
x = r*cos(t) + h;
y = r*sin(t) + k;
plot(x,y,'k');
end
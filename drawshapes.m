function drawshapes(shape,origin,diameter,fill,linewidth,rotation,name,identity,options)
% DRAWSHAPES: draw pre-defined shape with user-defined parameters
%
% Syntax:
% drawshapes(shape,origin,diameter,fill,linewidth)
%
% Input:
% shape: string, available options are: 'Filled Circle', 'Filled Square',
% 'Crossed Square', 'Divided Diamond', 'Filled Triangle', 'Divided
% Triangle', 'Flat Rectangle', 'Filled Star', 'Filled Diamond', 'Flat
% Hexagon', 'Pentagon'.
% origin: 1 x 2 vector, the center point of the shape.
% diameter: number, the diameter of the shape.
% fill: string, the filling color of the shape. Available opitons are:
% 'White', 'Blue', 'Green', 'Yellow', 'Light blue', 'Pink', 'Purple',
% 'Brown', 'Orange', 'Red'.
% linewidth, number, the line width of the shape perimeter
%
% Output:
% N/A
%
% Note:
% N/A
%
% Example:
% figure, hold on, drawshapes('flat rectangle',[0,0],3,'light blue',3,45), set(gca,'xlim',[-3,3]), axis equal
%
% Children function:
% N/A
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2016, Research Foundation for State University of New York. All rights reserved
%


% colorname = {'white','blue','green','yellow','light blue','pink','purple','brown','orange','red'};
% colorval = [1 1 1;0 .564 .737;0 .651 0.318;1 .831 0;.56 .8 .914;0.964 0.619 0.631;.647 .262 0.6;.631 .478 .302;.956 .471 .125;.929 .11 .141];
% fill = colorval(ismember(lower(colorname),lower(fill)),:);
fill = MonoColor.(fill);
color = [0 0 0];
theta = rotation/180*pi;
switch identity
    case 'SHAPE'
        if strcmpi(shape,'filled circle')
            pos = [origin(1) - .5*diameter,origin(2)-.5*diameter,diameter,diameter];
            rectangle('Position',pos,'Curvature',[1 1],'facecolor',fill,'edgecolor',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Filled Square')
            r = sqrt(2)/2*diameter;
            vertexes = [origin(1)+r*cos(1/4*pi-theta),origin(2)+r*sin(1/4*pi-theta);...
                origin(1)+r*cos(1/4*pi+theta),origin(2)-r*sin(1/4*pi+theta);...
                origin(1)-r*cos(1/4*pi-theta),origin(2)-r*sin(1/4*pi-theta);...
                origin(1)-r*cos(1/4*pi+theta),origin(2)+r*sin(1/4*pi+theta);...
                origin(1)+r*cos(1/4*pi-theta),origin(2)+r*sin(1/4*pi-theta)];
            patch(vertexes(1:4,1),vertexes(1:4,2),fill,'edgecolor',[1 1 1])
            plot(vertexes(:,1),vertexes(:,2),'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Crossed Square')
            r = sqrt(2)/2*diameter;
            vertexes = [origin(1)-r*cos(1/4*pi+theta),origin(2)+r*sin(1/4*pi+theta);...
                origin(1)+r*cos(1/4*pi-theta),origin(2)+r*sin(1/4*pi-theta);...
                origin(1)+r*cos(1/4*pi+theta),origin(2)-r*sin(1/4*pi+theta);...
                origin(1)-r*cos(1/4*pi-theta),origin(2)-r*sin(1/4*pi-theta);...
                origin(1)-r*cos(1/4*pi+theta),origin(2)+r*sin(1/4*pi+theta)];
            patch(vertexes(1:3,1),vertexes(1:3,2),fill,'edgecolor',[1 1 1])
            patch(vertexes(3:5,1),vertexes(3:5,2),'white','edgecolor',[1 1 1])
            plot([vertexes(:,1);vertexes(3,1)],[vertexes(:,2);vertexes(3,2)],'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Flat Rectangle')
            r = sqrt(5)/4*diameter;
            alpha = asin(1/sqrt(5));
            vertexes = [origin(1)+r*cos(alpha-theta),origin(2)+r*sin(alpha-theta);...
                origin(1)+r*cos(alpha+theta),origin(2)-r*sin(alpha+theta);...
                origin(1)-r*cos(alpha-theta),origin(2)-r*sin(alpha-theta);...
                origin(1)-r*cos(alpha+theta),origin(2)+r*sin(alpha+theta);...
                origin(1)+r*cos(alpha-theta),origin(2)+r*sin(alpha-theta)];
            patch(vertexes(1:4,1),vertexes(1:4,2),fill,'edgecolor',[1 1 1])
            plot(vertexes(:,1),vertexes(:,2),'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Filled Diamond')
            r = diameter/2;
            vertexes = [origin(1)-r*cos(theta),origin(2)+r*sin(theta);...
                origin(1)+r*sin(theta),origin(2)+r*cos(theta);...
                origin(1)+r*cos(theta),origin(2)-r*sin(theta);...
                origin(1)-r*sin(theta),origin(2)-r*cos(theta);...
                origin(1)-r*cos(theta),origin(2)+r*sin(theta)];
            patch(vertexes(1:4,1),vertexes(1:4,2),fill,'edgecolor',[1 1 1])
            plot(vertexes(:,1),vertexes(:,2),'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Divided Diamond')
            r = diameter/2;
            vertexes = [origin(1)-r*cos(theta),origin(2)+r*sin(theta);...
                origin(1)+r*sin(theta),origin(2)+r*cos(theta);...
                origin(1)+r*cos(theta),origin(2)-r*sin(theta);...
                origin(1)-r*sin(theta),origin(2)-r*cos(theta);...
                origin(1)-r*cos(theta),origin(2)+r*sin(theta)];
            patch(vertexes(1:3,1),vertexes(1:3,2),fill,'edgecolor',[1 1 1])
            patch(vertexes(3:5,1),vertexes(3:5,2),'white','edgecolor',[1 1 1])
            plot([vertexes(:,1);vertexes(3,1)],[vertexes(:,2);vertexes(3,2)],'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Filled Triangle')
            r = sqrt(2)/2*diameter;
            vertexes = [origin(1)+r/sqrt(2)*sin(theta),origin(2)+r/sqrt(2)*cos(theta);...
                origin(1)+r*cos(pi/4+theta),origin(2)-r*sin(pi/4+theta);...
                origin(1)-r*cos(pi/4-theta),origin(2)-r*sin(pi/4-theta)];
            patch(vertexes(:,1),vertexes(:,2),fill,'edgecolor',[1 1 1])
            plot([vertexes(:,1);vertexes(1,1)],[vertexes(:,2);vertexes(1,2)],'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Divided Triangle')
            r = sqrt(2)/2*diameter;
            vertexes = [origin(1)+r/sqrt(2)*sin(theta),origin(2)+r/sqrt(2)*cos(theta);...
                origin(1)+r*cos(pi/4+theta),origin(2)-r*sin(pi/4+theta);...
                origin(1)-r/sqrt(2)*sin(theta),origin(2)-r/sqrt(2)*cos(theta);...
                origin(1)-r*cos(pi/4-theta),origin(2)-r*sin(pi/4-theta)];
            patch(vertexes(1:3,1),vertexes(1:3,2),fill,'edgecolor',[1 1 1])
            patch(vertexes([1,3,4],1),vertexes([1,3,4],2),'white','edgecolor',[1 1 1])
            plot([vertexes(:,1);vertexes(1,1)],[vertexes(:,2);vertexes(1,2)],'color',color,'linewidth',linewidth)
            plot([vertexes(1,1),vertexes(3,1)],[vertexes(1,2),vertexes(3,2)],'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Filled Star')
            vertexes = zeros(10,2);
            for i = 1:5
                vertexes(i*2-1,:) = origin + diameter/2*[sin(2*pi/5*(i-1)+theta),cos(2*pi/5*(i-1)+theta)];
                vertexes(i*2,:) = origin + diameter*sin(0.1*pi)/sin(0.7*pi)/2*[sin(0.4*pi*i-.2*pi+theta),cos(0.4*pi*i-.2*pi+theta)];
            end
            patch(vertexes(:,1),vertexes(:,2),fill,'edgecolor',[1 1 1])
            plot([vertexes(:,1);vertexes(1,1)],[vertexes(:,2);vertexes(1,2)],'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Flat Hexagon')
            r1 = sqrt(2)/4*diameter;
            r2 = diameter/2;
            vertexes = [r1*cos(pi/4-theta),r1*sin(pi/4-theta);...
                r2*cos(theta),-r2*sin(theta);...
                r1*cos(pi/4+theta),-r1*sin(pi/4+theta);...
                -r1*cos(pi/4-theta),-r1*sin(pi/4-theta);...
                -r2*cos(theta),r2*sin(theta);...
                -r1*cos(pi/4+theta),r1*sin(pi/4+theta);...
                r1*cos(pi/4-theta),r1*sin(pi/4-theta)];
            vertexes = vertexes + repmat(origin,size(vertexes,1),1);
            patch(vertexes(1:6,1),vertexes(1:6,2),fill,'edgecolor',[1 1 1])
            plot([vertexes(:,1);vertexes(1,1)],[vertexes(:,2);vertexes(1,2)],'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Pentagon')
            vertexes = zeros(5,2);
            for i = 1:5
                vertexes(i,:) = origin + diameter/2*[sin(2*pi/5*(i-1)+theta),cos(2*pi/5*(i-1)+theta)];
            end
            patch(vertexes(:,1),vertexes(:,2),fill,'edgecolor',[1 1 1])
            plot([vertexes(:,1);vertexes(1,1)],[vertexes(:,2);vertexes(1,2)],'color',color,'linewidth',linewidth)
        end
        
        if strcmpi(shape,'Flattened Diamond')
            r = diameter/2;
            vertexes = [0,r/sqrt(3);...
                r,0;...
                0,-r/sqrt(3);...
                -r,0;...
                0,r/sqrt(3)];
            vertexes = vertexes + repmat(origin,size(vertexes,1),1);
            patch(vertexes(:,1),vertexes(:,2),fill,'edgecolor',[1 1 1])
            plot([vertexes(:,1);vertexes(1,1)],[vertexes(:,2);vertexes(1,2)],'color',color,'linewidth',linewidth)
        end
    case 'CHAR'
        text(origin(1),origin(2),name,'FontSize',options.fontsize,'FontWeight','bold','HorizontalAlignment','center','BackgroundColor','w','Margin',1)
end


% shape,origin,diameter,fill,linewidth,rotation,name,identity)


end
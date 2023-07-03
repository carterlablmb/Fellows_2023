function AlexGUI
% This GUI allows a user to load a file output from AlexStepFindIntensity
% (all_fit_???.txt).  The user can identify points (relative frame,
% absolute frame, x,y (pixels or um?), intensity),  delete
% steps.  The user then presses Save,  new step positions are saved in an 8th column, DiscardChanges or
% presses the next/previous button to move onto the next image (which will
% automatically save any edits)


%% Clean up
clear all
%% Variables
%Set PretendNew to 1 to force program to read data as new:
PretendNew=0;
%% Define global variables
global DATA
global OLDNEW
global frm
global sfr
global int
global xpos
global ypos
global step
global FileStatus
%% More clean up:
close all
%% Load some data

%get file name via input box
[FileName1, Dir] = uigetfile({'all_fit*.txt'});
if FileName1==0  % loop closes script if cancel is pressed
    %error('File Open Cancelled')
    disp('File Open Cancelled');
    return  % can use this to exit gracefully without error messages
end
cd(Dir)
%cd (Dir) %move to the directory
DATA=dlmread(FileName1);
[~,FileName,~]=fileparts(FileName1);
FileName=erase(FileName,'all_fit_');
%may need to deal with NaN in the input here

if length(DATA(1,:))==7  ||  PretendNew==1 %New file or if you use PretendNew override
    frm    = DATA(:,1);  % read first column into vector t etc
    sfr = DATA(1,2)+1;
    xpos  = DATA(1,3);
    ypos  = DATA(1,4);
    int    = DATA(:,5);
    step   = DATA(:,6);
    OLDNEW = 1;
    
elseif length(DATA(1,:))==8  %Previously read file
    frm    = DATA(:,1);  % read first column into vector t etc
    sfr = DATA(1,2)+1;
    xpos  = DATA(1,3);
    ypos  = DATA(1,4);
    int    = DATA(:,5);
    step   = DATA(:,8);
    OLDNEW = 0;
else  %if DATA is not 7 or 8 columns long then stop
    disp('File not recognised:');
    return
end


%% Initialise UNDOSTEP and UNDOSCORE and set TextStatus
UNDOSTEP=zeros(length(frm),1);
WriteStatus=1 ; %Text is off
%mergestart=0;
%deletestart=0;
%stepstartid=0;
%LastID=0;
%plotter=1; %Good to plot...
FileList=dir(fullfile(Dir,'all_fit*.txt'));  %get a list of all files in current directory starting with all_fit
for w = 1:length(FileList)
    if strcmp(FileList(w).name,FileName1)==1
        FileStatus=w; % Value of 0 means the First button has not been pressed yet.
        break
    end
end

%% figure(1)


figure('Name','x y plot', ...
    'Numbertitle','on', ...
    'Menubar','none', ...
    'ToolBar','figure', ...
    'Position',[100 100 525 550], ...
    'KeyPressFcn',@keyStrokes,...
    'WindowButtonDownFcn',@clickButtonDown);
%'WindowButtonUpFcn',@clickButtonUp, ...
%'WindowButtonMotionFcn',@buttonMotion);
axes('Units','normalized')
movegui(gcf,'northwest')
%scrsz = get(0,'ScreenSize');
%set(gcf,'position',[10 scrsz(4)*(0.4) scrsz(3)*(500/1400) scrsz(4)*0.5])
%set(gca,'position',[.03  .05  .95  .85]);


%% Control panel
% radio buttons for choosing to draw or delete
ddpanel = uipanel('BackgroundColor',[0.9 0.9 0.9], ...
    'Position',[0.01 0.92 0.98 0.07]);

% text box to show x/y or vertex id
figure_title = uicontrol(ddpanel, ...
    'Style','text', ...
    'Units','normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.01 0.1 0.5 0.7], ...
    'FontSize',11, ...
    'ForegroundColor','k', ...
    'FontWeight','b', ...
    'BackgroundColor',[0.9 0.9 0.9]);


identify_radio = uicontrol(ddpanel, ...
    'Style','radiobutton', ...
    'String','I(d)ent', ...
    'FontSize',11, ...
    'Units','normalized', ...
    'BackgroundColor',[0.9 0.9 0.9], ...
    'Value',1, ...
    'Position', [0.45 0.1 0.2 0.8], ...
    'Callback',@clickIdentifyRadio);

deletestep_radio = uicontrol(ddpanel, ...
    'Style','radiobutton', ...
    'String','Delete(s)', ...
    'FontSize',11, ...
    'Units','normalized', ...
    'BackgroundColor',[0.9 0.9 0.9], ...
    'Value',0, ...
    'Position', [0.5 0.1 0.2 0.8], ...
    'Callback',@clickDeleteStepRadio);


% text box to show x/y or vertex id
xy_txt1 = uicontrol(ddpanel, ...
    'Style','text', ...
    'Units','normalized', ...
    'HorizontalAlignment','left',...
    'Position',[0.6 0.1 0.35 0.65], ...
    'FontSize',11, ...
    'FontWeight','b', ...
    'BackgroundColor',[0.9 0.9 0.9]);

% button to load all files with the format "all_fit*.txt in the current
% directory
uicontrol(ddpanel, ...
    'Style','pushbutton', ...
    'String','1st', ...
    'Units','normalized', ...
    'Position', [0.70 0.1 0.04 0.8], ...
    'Callback',@loadFirst);

% button to load next files with the format "all_fit*.txt in the list of
% those in the current directory
uicontrol(ddpanel, ...
    'Style','pushbutton', ...
    'String','Prv', ...
    'Units','normalized', ...
    'Position', [0.74 0.1 0.04 0.8], ...
    'Callback',@loadPrev);

% button to load next files with the format "all_fit*.txt in the list of
% those in the current directory
uicontrol(ddpanel, ...
    'Style','pushbutton', ...
    'String','Nxt', ...
    'Units','normalized', ...
    'Position', [0.78 0.1 0.04 0.8], ...
    'Callback',@loadNext);

uicontrol(ddpanel, ...
    'Style','checkbox', ...
    'String','Write', ...
    'Units','normalized', ...
    'Value',1, ...
    'Position', [0.82 0.1 0.03 0.8], ...
    'BackgroundColor',[0.9 0.9 0.9],...
    'Callback',@toggleWrite);

% button to save data and output a summary file
uicontrol(ddpanel, ...
    'Style','pushbutton', ...
    'String','Save', ...
    'Units','normalized', ...
    'Position', [0.86 0.1 0.06 0.8], ...
    'Callback',@saveGraph);

% button to undo last change to step
uicontrol(ddpanel, ...
    'Style','pushbutton', ...
    'String','Undo', ...
    'Units','normalized', ...
    'Position', [0.92 0.1 0.06 0.8], ...
    'Callback',@undoStepChange);


%% Initial plotting: clicking with the mouse will update plots

%Plot graphs in figure(1) 
plotInitialGraphs(frm,int,step)
%Put title of plot into title bar
set(figure_title,'String',['>' sprintf(FileName) ' Start: fr_' num2str(sfr) ' , x_' num2str(xpos,'%.1f') ' , y_' num2str(ypos,'%.1f')]);



%% Plot figures 1 and 2 (editable figures)
    function plotInitialGraphs(frm,int,step)
        
        
        %% set axis values (x=relative frame, y=intensity
        cla;
        xmin = round(min(frm))-10;
        xmax = round(max(frm))+10;
        ymin = round(min(int))-5;
        ymax = round(max(int))+5;
        ax1 = [xmin xmax ymin ymax];

        
        %% plot figure 1: xy plot
        figure(1)
        %Plot raw data
        hold off % first plot command will clear anything else.
        
        %1st and 5th coln of DATA are frame and intensity: plot appears much faster than
        %scatter so am using this command for points: see profile command
        %for speed info

            XYr =plot (frm,int,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);
  
        hold on % allows a second line to be put on the same plot
        %black line through points
        XYk=plot (frm,int,'k');
        %XYk=plot (1,1,'k');  %Leave this in to get the title to appear... Odd
        %Plot steps (column6)

        if OLDNEW==1
        XYSm=plot (frm,step,'s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8);
        else
        XYSm=plot (frm,step,'s','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
        end

        %Plot markers used in identify script
        h1=plot(1,1,'co','MarkerSize',10,'LineWidth',2);
        
        %% Add text labels to figure 1 with time of last point of step:

   
        %TitleName=regexprep(FileName, '_', '-');
        %title(TitleName)
        %label axes
        xlabel('relative frame')
        ylabel('intensity')
        
        %set grid lines
        %set(gca,'XTick',xmin:8:xmax); %because on this plot x is on the y axis.
        %set(gca,'YTick',ymin:6:ymax);
        grid on
        %sets axis to fit data
        %axis(ax1)
        
        
        %% Get handles for each figure for use in ButtonDown function
        %handles are attached to figures and called when the figure is clicked
        %after making myhandle, the handles for each marker are added to it.  They will be
        %called on mouse down
        
        %Handles for figure 1
        myhandle1=guihandles(1);
        myhandle1.marker=h1;
        myhandle1.XYr=XYr;
        myhandle1.XYk=XYk;
        myhandle1.XYSm=XYSm;
        guidata(1,myhandle1); %save the above data (called using fig 1 handle in ButtonDown and Redraw scripts
        
    end
%% Control clicking on radio buttons
    function clickIdentifyRadio(varargin)
        set(identify_radio,'Value',1);
        set(deletestep_radio,'Value',0);
    end

    function clickDeleteStepRadio(varargin)
        set(identify_radio,'Value',0);
        set(deletestep_radio,'Value',1);
    end
   

%% Function responding to key strokes when the figure is in focus.
    function keyStrokes(src,evnt)
        if evnt.Character == 'd'
            set(identify_radio,'Value',1);
            set(deletestep_radio,'Value',0);
        elseif  evnt.Character == 's'
            set(identify_radio,'Value',0);
            set(deletestep_radio,'Value',1);
        end
    end
%% Find nearest point

    function pt = getPoint()
        cpt = get(gca,'CurrentPoint');
        %pt = round(cpt(1));
        pt = cpt(1,1:2);
    end



%% Find the number of the nearest point to where you are
    function id = getNearestVertex(pt, frm, int)       
         xy=[frm, int];
         dsq = (sum((xy-pt(ones(size(xy,1),1),:)).^2,2));
         id = find(dsq==min(dsq));
         id = id(1);
    end

%% Button Down:

%If identify is hit, then just plots position and gives t,x,y
%If a

    function clickButtonDown(varargin)
        pt = getPoint;
        id = getNearestVertex(pt, frm, int);
        %  LastID = id;
        ax = axis;
        %get(gcbf, 'Name')  %will give name of figure doing the calling
        if pt(1)>=ax(1) && pt(1)<=ax(2) && pt(2)>=ax(3) && pt(2)<=ax(4)  %keeps point on axis area

            %Identify Point
            if get(identify_radio,'Value') % add point
                set(xy_txt1,'String',[num2str(frm(id)) ': ' num2str(int(id),'%.1f') ': ' num2str(step(id),'%.1f')])
                cyanMarker(id, DATA);  %plots cyan marker, deletes old one.


            %Delete Step.
            elseif get(deletestep_radio,'Value') %
                %Put value of step into UNDO matrices:
                %following works as long as UNDOSTEP is the correct length and is
                %already available (ie full of zeros...)
                set(xy_txt1,'String',[num2str(frm(id)) ' : ' num2str(int(id),'%.1f')])
                cyanMarker(id, DATA);  %plots cyan marker, deletes old one.
                UNDOSTEP=[step,UNDOSTEP];
                value=step(id);
                dellist = step==value;
                step(dellist)=NaN;
                
              
                %Call reploting function: which will delete plots and remake them
                plotReplot(frm,int,step);
                
    
                              
           end
        end
    end
%% Move marker buttons around

    function cyanMarker(id, DATA)
        myhandle1 = guidata(1); %get handles from fig1: see PlotGraph

        % remove marker from fig1
        delete(myhandle1.marker);

        figure(1)
        hold on
        myhandle1.marker=plot (DATA(id,1),DATA(id,5),'co','MarkerSize',10,'LineWidth',2);
        %hold off
        guidata(1,myhandle1)% save new handle to plot
    end

%% Save data on pressing the Save button

    function saveGraph(varargin)
        %  Replace input file with contents of DATA
        DATA=[DATA(:,1:7) step];
        dlmwrite(FileName1, DATA, 'delimiter', '\t')
        disp(['Saved data to 8th column of: ' FileName1])
    end

%% Undo the last change to steps
    function undoStepChange(varargin)
        %Change DATA back to the first column of UNDOSTEP
        if length(UNDOSTEP(1,:))>1  %only allow undo if there is more than one column
            step=UNDOSTEP(:,1);
            UNDOSTEP(:,1) = [];
            
            %%Call reploting function: which will delete plots and remake
            %%them, as well as calculating the new sum table.
            plotReplot(frm,int,step);
        end
    end


%% Plot figure 1
    function plotReplot(frm,int,step)
        %Get handles to objects in two figures so that can be deleted etc:
        myhandle1 = guidata(1); %get handles from fig1: see PlotGraph
        figure(1)
        hold on % allows a second line to be put on the same plo
       
            delete(myhandle1.XYr);
            myhandle1.XYr =plot (frm,int,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);

        % delete(myhandle1.XYk);
        %myhandle1.XYk=plot (1,1,'k');  %Leave this in to get the title to appear... Odd
         if OLDNEW==1
        delete(myhandle1.XYSm);
        myhandle1.XYSm=plot (frm,step,'s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8);
         else
        delete(myhandle1.XYSm);
        myhandle1.XYSm=plot (frm,step,'s','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',8);
         end
        %% Will not update text at the moment:
        %%%hold off
        guidata(1,myhandle1);% save new handle to plot
    end

  function toggleWrite(varargin)
        %myhandle1 = guidata(1); %get handles from fig1: see PlotGraph
        if WriteStatus==1  %if line is on
            WriteStatus=0;
        else
            WriteStatus=1;
        end
  end

%% Script to load the first file in the current director fitting all_fit*.txt
    function loadFirst(varargin)
        %If the Write checkbox is ticked, write the data on the current
        %file.
        if WriteStatus==1
            saveGraph
        end
            %FileName1 gets the first file on the list
            FileStatus=1;
            FileName1=FileList(FileStatus).name;
            %Reset all variables
            DATA=[];
            frm=[];
            int=[];
            step=[];
            UNDOSTEP=[];
        DATA=dlmread(FileName1);
        [~,FileName,~]=fileparts(FileName1);
                  FileName=erase(FileName,'all_fit_');            
            %get variables out of DATA, calculate them if necessary
            %may need to deal with NaN in the input here
           
            if length(DATA(1,:))==7  ||  PretendNew==1 %New file or if you use PretendNew override
                frm    = DATA(:,1);  % read first column into vector t etc
                sfr = DATA(1,2)+1;
                xpos  = DATA(1,3);
                ypos  = DATA(1,4);
                int    = DATA(:,5);
                step   = DATA(:,6);
                OLDNEW = 1;
                
            elseif length(DATA(1,:))==8  %Previously read file
                frm    = DATA(:,1);  % read first column into vector t etc
                sfr = DATA(1,2)+1;
                xpos  = DATA(1,3);
                ypos  = DATA(1,4);
                int    = DATA(:,5);
                step   = DATA(:,8);
                OLDNEW = 0;
            else  %if DATA is not 7 or 8 columns long then stop
                disp('File not recognised:');
                return
            end
            % Initialise UNDOSTEP and UNDOSCORE and set TextStatus
            
            UNDOSTEP=zeros(length(frm),1);
            % Initial plotting: clicking with the mouse will update plot
            %Plot graphs in figure(1)
            plotInitialGraphs(frm,int,step)
            %Put title of plot into title bar (and the start position)
            set(figure_title,'String',['>' sprintf(FileName) ' Start: fr_' num2str(sfr) ' , x_' num2str(xpos,'%.1f') ' , y_' num2str(ypos,'%.1f')]);
        end

%% LoadNext Button: loads the next file labelled all_fit*.txt in the current directory

    function loadNext(varargin)
        if  FileStatus<length(FileList)
            %If the Write checkbox is ticked, write the data on the current
            %file.
            if WriteStatus==1
                saveGraph
            end
            %FileName1 next file on the list
            FileStatus=FileStatus+1;
            FileName1=FileList(FileStatus).name;
            %Reset all variables
            DATA=[];
            frm=[];
            int=[];
            step=[];
            UNDOSTEP=[];
            %Read data
            DATA=dlmread(FileName1);
            [~,FileName,~]=fileparts(FileName1);
            %get variables out of DATA, calculate them if necessary
          FileName=erase(FileName,'all_fit_');            
            %get variables out of DATA, calculate them if necessary
            %may need to deal with NaN in the input here
            if length(DATA(1,:))==7  ||  PretendNew==1 %New file or if you use PretendNew override
                frm    = DATA(:,1);  % read first column into vector t etc
                sfr = DATA(1,2)+1;
                xpos  = DATA(1,3);
                ypos  = DATA(1,4);
                int    = DATA(:,5);
                step   = DATA(:,6);
                OLDNEW = 1;
                
            elseif length(DATA(1,:))==8  %Previously read file
                frm    = DATA(:,1);  % read first column into vector t etc
                sfr = DATA(1,2)+1;
                xpos  = DATA(1,3);
                ypos  = DATA(1,4);
                int    = DATA(:,5);
                step   = DATA(:,8);
                OLDNEW = 0;
            else  %if DATA is not 7 or 8 columns long then stop
                disp('File not recognised:');
                return
            end
            % Initialise UNDOSTEP 
            UNDOSTEP=zeros(length(frm),1);
            % Initial plotting: clicking with the mouse will update plot
            %Plot graphs in figure(1)
            plotInitialGraphs(frm,int,step)
            %Put title of plot into title bar
            set(figure_title,'String',['>' sprintf(FileName) ' Start: fr_' num2str(sfr) ' , x_' num2str(xpos,'%.1f') ' , y_' num2str(ypos,'%.1f')]);



        end
    end

%% LoadPrev Button: loads the previous file labelled all_fit*.txt in the current directory

    function loadPrev(varargin)
        if FileStatus>1 && FileStatus<=length(FileList)
            %If the Write checkbox is ticked, write the data on the current
            %file.
            if WriteStatus==1
                saveGraph
            end
            %FileName1 gets the previous file on the list
            FileStatus=FileStatus-1;
            FileName1=FileList(FileStatus).name;
            %Reset all variables
            DATA=[];
            frm=[];
            int=[];
            step=[];
            UNDOSTEP=[];
            %Read data
            DATA=dlmread(FileName1);
            [~,FileName,~]=fileparts(FileName1);
            FileName=erase(FileName,'all_fit_');
            %get variables out of DATA, calculate them if necessary
            %may need to deal with NaN in the input here
            if length(DATA(1,:))==7  ||  PretendNew==1 %New file or if you use PretendNew override
                frm    = DATA(:,1);  % read first column into vector t etc
                sfr = DATA(1,2)+1;
                xpos  = DATA(1,3);
                ypos  = DATA(1,4);
                int    = DATA(:,5);
                step   = DATA(:,6);
                OLDNEW = 1;
                
            elseif length(DATA(1,:))==8  %Previously read file
                frm    = DATA(:,1);  % read first column into vector t etc
                sfr = DATA(1,2)+1;
                xpos  = DATA(1,3);
                ypos  = DATA(1,4);
                int    = DATA(:,5);
                step   = DATA(:,8);
                OLDNEW = 0;
            else  %if DATA is not 7 or 8 columns long then stop
                disp('File not recognised:');
                return
            end
            % Initialise UNDOSTEP 
            UNDOSTEP=zeros(length(frm),1);
            % Initial plotting: clicking with the mouse will update plot
            %Plot graphs in figure(1)
            plotInitialGraphs(frm,int,step)
            %Put title of plot into title bar
           set(figure_title,'String',['>' sprintf(FileName) ' Start: fr_' num2str(sfr) ' , x_' num2str(xpos,'%.1f') ' , y_' num2str(ypos,'%.1f')]);

        end
    end

end





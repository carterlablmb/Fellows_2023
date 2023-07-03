function AlexHarvest
%% Open up an all_fit.txt file via a dialog box.
% The script will make a list of all all_fit.txt files in that directory.  It will open each one
% check that it contains 8 columns (i.e. has been checked via AlexGUI) and
% then pulls out the intensity values.  For each "step" in intensity it
% will write the value, the SD, the X,Y,real frame value and the name of
% the file.
DATA=[];
global OUTDATA
global OUTNAMES
OUTDATA=[];
OUTNAMES=[];

%%get file name via input box (all this really does is change to the
%%correct directory.
[Open, Dir] = uigetfile({'all_fit*.txt'});
if Open==0  % loop closes script if cancel is pressed
    %error('File Open Cancelled')
    disp('File Open Cancelled');
    return  % can use this to exit gracefully without error messages
end
cd(Dir)
%cd (Dir) %move to the directory
%get a list of all files in current directory starting with all_fit
FileList=dir(fullfile(Dir,'all_fit*.txt'));

for id = 1:length(FileList)
    FileName1=FileList(id).name;
    DATA=dlmread(FileName1);
    [~,FileName,~]=fileparts(FileName1);
    FileName=erase(FileName,'all_fit_');
    if length(DATA(1,:))==8 % check DATA contains 8 columns (steps are in col 8)
        %clear NaNs from column 8
        DATA = DATA(any(~isnan(DATA(:,8)),2),:);
        %Cut out redundant columns (move steps to col 5)
        DATA=DATA(:,[2:4 7:8]);
        for w=2:length(DATA(:,1))
            if DATA(w,5)==DATA(w-1,5)  % if value in col5 is the same as previous set line 1 to NaN
                DATA(w,1)=NaN;
            end

        end
        %get rid of lines containing NaNs in column 1 - cut down to just
        %steps
        DATA = DATA(any(~isnan(DATA(:,1)),2),:);
        disp([FileName ' - processed'])

        %Make a separate array containing the FileName for each
        NAMES=strings([length(DATA(:,1)),1]);
        NAMES(:)=FileName;

        %Append DATA into OUTDATA and NAMES into OUTNAMEs
        OUTDATA=cat(1,OUTDATA,DATA);
        OUTNAMES=cat(1,OUTNAMES,NAMES);
    else
        disp([FileName ' NOT PROCESSED - use AlexGUI to check it first'])
    end
end
%%Write out
%Default Name is the Name of current directory _harvest.csv
disp('Writing out steps from all_fit.txt files in the current directory')
[~,dirname]=fileparts(pwd);
[newfile,newpath,~] = uiputfile([dirname '_harvest.csv']);

OUTDATA=OUTDATA(:,[5 4 1 2 3]);
OUTPUT=[OUTDATA OUTNAMES];
OUTPUT=cat(1,["Step","Step SD","Frame","X postion (um)","Y position (um)","Filename"],OUTPUT);
foutname = [newpath,newfile];
writematrix(OUTPUT,foutname)
disp(['Wrote file: ' foutname])
end
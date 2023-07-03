function Bleach1_TifSplit()
%Reads all image.tif files in a directory as long as they have
%ident#image.csv file.
%The ident#image.csv file contains the output from TrackMate.
%ident is a marker of the type of track (diffusive, processive etc)

%The script will identify the area covered by each tracked spot and then
%cut it out from the movie and write it out as separate tifstack.  NB The
%frames will be relative - i.e the first image will be the frame that the
%spot starts beign tracked in.

%% Options for running programme
runall=0; % set run all to '1' to pick all .csv files in a directory
count=50; % update status ever 'count' lines of each .csv file.
padwidth=5; %padding round central pixel of a spot (too narrow and gaussian works less well)
offset=1; %to account for the fact matlab counts pixels from 1, imageJ from 0
defaultscale=1; %use this value if it is not possible to get it from the tifinfo
%% Load some data

%get file name via input box
[TrackMate, Dir] = uigetfile({'*.csv'});
if TrackMate==0  % loop closes script if cancel is pressed
    %error('File Open Cancelled')
    disp('File Open Cancelled');
    return  % can use this to exit gracefully without error messages
end
cd(Dir)

%if we want to run on all .csv files, collect those in the selected
%directory.  If not just run on the selected file.
if(runall==1)
    files = dir('*.csv');
else
    files=dir(TrackMate);
end
%Run through all the selected files
for k=1:length(files)
    TrackMate=files(k).name;
    %reset handles for each file.
    handles=[];
    handles.count=count;
    handles.width=padwidth;
    handles.offset=offset;
    %run analysis on each file.
    try
        launchTrack(TrackMate,handles);
    catch
        disp(['Following file could not run or ended early: ' TrackMate])
        return
    end

end
disp('Finished successfully')
%%Main function
    function launchTrack(TrackMate,handles)
            [~,f,~]=fileparts(TrackMate);
            %If there is an identifier (#) then get the text after it to be the image
            %name.  Otherwise take everything.
            a=split(f,'#');
            if length(a)==2
                %ident=[a{1}];
                imagename=[a{2}];
            elseif length(a)>2
                disp('.csv file must contain a single # character, marking an identifier')
                return
            else
                %ident=[];
                imagename=f;
            end

            handles.tiffile=[imagename,'.tif'];





            %% Variables
            %
            % pixelsize=0.10905;
            %tiffile='AlexMovie.tif';

            try  % Check if there is a tif file.
                info=[];
                info=imfinfo(handles.tiffile);
            catch
                disp(['No file:' handles.tiffile])
                return
            end
            handles.tiflength=size(info,1);
            scale=info.XResolution;
            warning('off'); %switches off tifread warning about pixel
            %If scale is [] then set it to default value instead
            if isempty(scale)
                scale=defaultscale;
            end

            % Create structure handles that can be accessed by all. Need these
            % variables to start other functions below.
            %handles=[];
            handles.filename=f;
            handles.scale=scale;
            handles.frame=[];
            handles.track=[];
            handles.x = [];
            handles.y = [];
            handles.cutsx= [];
            handles.cutsy= [];
            handles.cutex= [];
            handles.cutey= [];
            handles.intensity=[];
            handles.file=[];


            %% Open Trackmate file.

            fid =[];
            fid = fopen(TrackMate,'r');
            Header=textscan(fid,'%s',1,'delimiter','\n')  ;
            %The next use of textscan will now miss the top 4 lines
            Tmate=[]; %initialize Tmate
            Tmate=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
            %This will read in the 2 fields (strings) and put them into cell arrays.
            fclose(fid)
            %format Tmate as a double array
            Tmate=[Tmate{1:20}];
            Tmate=Tmate(:,[3,5:6]); %cut the format down to allow the rest to work
            Tmate=str2double(Tmate);

            %Tmate=sortrows(Tmate,6);  %sort Tmate based on the frame number
            Tmate(:,2:3)=(Tmate(:,2:3)*scale)+handles.offset; %convert trackmate distance measure to pixels.
            %Tmate(:,2:3)=(Tmate(:,2:3)*scale); %convert trackmate distance measure to pixels.

            %Tmate = Tmate(any(~isnan(Tmate),2),:);  % only keep rows that don't contain NaN (ie remove NaNs)
            %Tmate=Tmate(Tmate(:,6)<tiflength,:); % remove any rows with a frame number higher than the number of frames in the movie. NB using "<", instead of "<=" because Tframes start at 0
            %spotnumber=length(Tmate(:,1));

            %% Work out the parts of each Track to cut out
            handles.trackID=Tmate(:,1);

            for i=1:length(handles.trackID)
                %work out the start and end points on the tif for cutting out each frame.
                handles.cutsx(i)=round(Tmate(i,2))-handles.width;
                handles.cutsy(i)=round(Tmate(i,3))-handles.width;
                handles.cutex(i)=round(Tmate(i,2))+handles.width;
                handles.cutey(i)=round(Tmate(i,3))+handles.width;
                handles.startframe(i)=1;
                handles.endframe(i)=handles.tiflength;
            end

            %% For each track load each tif frame, cut out correct section and put into tifstore, then write out the new mini tif.
            % I think this is better than trying to load everything first (but might be
            % wrong).  Need to clean up handles.cutsx etc to remove any
            % values outside the image size.
            imwidth=info.Width;
            imheight=info.Height;

            idsx1 = handles.cutsx < 1;
            %idsx2 = handles.cutsx > imwidth; %unlikely sx > width, but just in case?
            idsy1 = handles.cutsy < 1;
            %idsy2 = handles.cutsy > imheight;

            %idex1 = handles.cutex < 1;
            idex2 = handles.cutex > imwidth;
            %idey1 = handles.cutey < 1;
            idey2 = handles.cutey > imheight;

            handles.cutsx(idsx1)=1;
            handles.cutsy(idsy1)=1;
            %handles.cutex(idex1)=2;
            %handles.cutey(idey1)=2;

            %handles.cutsx(idsx2)=info.Width-1;
            %handles.cutsy(idsy2)=info.Height-1;
            handles.cutex(idex2)=imwidth;
            handles.cutey(idey2)=imheight;



            try
                for j=1:length(handles.trackID)
                    tkID=handles.trackID(j);
                    sx=handles.cutsx(j);
                    sy=handles.cutsy(j);
                    ex=handles.cutex(j);
                    ey=handles.cutey(j);
                    sfm=handles.startframe(j);
                    efm=handles.endframe(j);
                    %preallocate a multdimensional array to take all the tif data

                    tifstore=zeros(ey-sy+1,ex-sx+1,efm-sfm+1);

                    %Run through handles, load correct tif-frame and place in a structured
                    %array for later writing out.
                    ct=1;
                    for k=sfm:efm
                        %Read in the first frame of the image
                        tf=imformats('tif');
                        im = double(feval(tf.read, handles.tiffile, k));
                        tifstore(:,:,ct)=im(sy:ey, sx:ex);
                        ct=ct+1;
                    end

                    %% Tifdata

                    [height, width, depth] = size(tifstore);
                    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                    tagstruct.ImageLength = height;
                    tagstruct.ImageWidth = width;
                    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                    % Grayscale image with real numbers
                    tagstruct.SamplesPerPixel = 1;
                    tagstruct.Compression = Tiff.Compression.None;
                    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
                    % tagstruct.SampleFormat = Tiff.SampleFormat.Int;
                    tagstruct.BitsPerSample = 16;

                    %write out tifstore
                    %open the tiff ('a' should append)
                    tfile=Tiff([f,'@',num2str(tkID),'.tif'], 'w');
                    %setTag(tfile,tagstruct)
                    %for d = 1:depth
                    %   write(tfile,uint16(tifstore(:, :, d)));
                    %end

                    for d = 1:depth
                        tfile.setTag(tagstruct);
                        tfile.write(uint16(tifstore(:, :, d)));
                        if d ~= depth
                            tfile.writeDirectory(); %this line allows file to be appended.
                        end
                    end
                    tfile.close();
                    writename=[f,'@',num2str(tkID),'.tif'];
                    disp(['Written: ' writename])

                end
            catch
                disp(['Failed to write:', f,'@',num2str(tkID),'.tif'])
            end
        end

end














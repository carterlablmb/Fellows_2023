function out=Alex2Bleach()
%Reads all image.tif files in a directory as long as they have
%ident#image.csv file.
%The ident#image.csv file contains the initial position of a spot that we
%want to measure the intensity of throughout the tifstack.
%ident is a marker of  
% 
%   the type of track (e.g. spots)

%The script will read in each frame and run Manu's gaussian fitter on all
%the tracks from TrackMate.  It will output files with "ident#image@n.txt"
%where n is the tracknumber.

%% Options for running programme
runall=0; % set run all to '1' to pick all .csv files in a directory
count=10; % update status ever 'count' lines of each .csv file.
padwidth=5; %padding round central pixel of a spot (too narrow and gaussian works less well)
offset=1; %to account for the fact matlab counts pixels from 1, imageJ from 0
defaultscale=1; %use this value if it is not possible to get it from the tifinfo
foldcutoff=10; %data more than X-fold above median is replaced by NaNs may want to change this!

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
handles.defaultscale=defaultscale;
handles.foldcutoff=foldcutoff;
%run analysis on each file.
try
launchTrack(TrackMate,handles);
catch
disp(['Could not run: ' TrackMate])  
end
    
end

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

tiffile=[imagename,'.tif'];


%% Variables
%
% pixelsize=0.10905;
%tiffile='AlexMovie.tif';

try  % Check if there is a tif file.
info=imfinfo(tiffile);
catch
    disp(['No file:' tiffile])
    return
end
tiflength=size(info,1);
scale=info.XResolution;
warning('off'); %switches off tifread warning about pixel
%If scale is [] then set it to default value instead
if isempty(scale)
    scale=handles.defaultscale;
end

% Create structure handles that can be accessed by all. Need these
% variables to start other functions below.
handles.filename=f;
handles.scale=scale;
handles.tiflength=tiflength;
handles.frame=[];
handles.track=[];
handles.x = [];
handles.y = [];
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
Tmate=Tmate(:,[1,5:6]); %cut the format down to allow the rest to work
Tmate=str2double(Tmate);

%Tmate=sortrows(Tmate,6);  %sort Tmate based on the frame number
Tmate(:,2:3)=(Tmate(:,2:3)*scale)+handles.offset; %convert trackmate distance measure to pixels.
%Tmate(:,2:3)=(Tmate(:,2:3)*scale); %convert trackmate distance measure to pixels.

%Tmate = Tmate(any(~isnan(Tmate),2),:);  % only keep rows that don't contain NaN (ie remove NaNs)
%Tmate=Tmate(Tmate(:,6)<tiflength,:); % remove any rows with a frame number higher than the number of frames in the movie. NB using "<", instead of "<=" because Tframes start at 0
spotnumber=length(Tmate(:,1));




%% Run through the Tmate array.
%% First frame
% For the first line get the initial guess position, read in the first
% image from the tif (to "im") the run the gaussian fit.  Put the frame, track and new x, y and
% amplitude into a new array "Gdata" that has been set up to be the same
% length as Tmate

%ouputlength=tiflength*spotnumber; %output length: in case we want to preassign the handles length?
%read in first tif image, find all spots, put gaussian values into handles.
% Then move onto the next tif image
ol=1;
i=1;
%while i<=100
while i<=tiflength
    j=1;
    while j<=spotnumber

        handles.track(ol)=(Tmate(j,1)); % First column gets TrackID
        handles.frame(ol)=i; % Second column gets Frame

        posGuess=[Tmate(j,2),Tmate(j,3)]; %posGuess gets the x,y postions from Tmate
        %Read in the first frame of the image
        tf=imformats('tif');
        im = double(feval(tf.read, tiffile, i));
        try          
            %Run Gaussian - if handles.width = 3 for the Window size (total
            %size of box) will be 7 (3+3+1)
            Wfit=1;
            if Wfit==0
                [pos,phot_count,~,~,~]=getPosGauss2d_ed(im, posGuess,  handles.width);
            else
                %Older fitting.
                guessx=posGuess(1);
                guessy=posGuess(2);
                temp=im(guessy-handles.width:guessy+handles.width , guessx-handles.width:guessx+handles.width);
                [fitpos,phot_count] = webbfit(temp,2);
                pos(1)=guessx-1-handles.width+fitpos(1);
                pos(2)=guessy-1-handles.width+fitpos(2);
            end
            %
            handles.x(ol)=pos(1);% GaussianFit X and Y position (in pixels for now)
            handles.y(ol)=pos(2);
            handles.intensity(ol)=phot_count; %Intensity
        catch % if it fails put NaN in.
            handles.x(ol)=NaN;
            handles.y(ol)=NaN;
            handles.intensity(ol)=NaN;
        end
        ol=ol+1;  %increment handles line number
        j=j+1;  % increment spot number
    end
    if mod(i,handles.count)== 0  %Display progress every 10 frames
        disp([handles.filename '.csv - frame: ' num2str(i)]);
    end
    i=i+1;

end



%% Clean up 
% Remove those spots for which handles.x(i) - Tmate(i,2) > 0 (rounded) and
% the same for y.

%handlelength=length(handles.x);
%medianintensity=median(handles.intensity);
%for c=1:handlelength
    %Remove outliers - by looking for spots which have shifted relative to
    %the input position or which are 10x bigger than the median value
%if(round(handles.x(c)-Tmate(c,2)))>0||(round(handles.y(c)-Tmate(c,3)))>0||handles.intensity(c)>(10*medianintensity)
    %disp(handles.intensity(c))
 %   handles.intensity(c)=NaN;
%end
%end

%% Write out data for each track

Tracks=unique(handles.track);
j=1;
while(j<=length(Tracks))
        trackid=Tracks(j);
        foutname = [handles.filename '@' num2str(trackid) '.txt'];
        fout =fopen(foutname,'w');
        rfr=1:length(handles.track(handles.track==trackid));  %relative frame 1 - length
        afr=handles.frame(handles.track==trackid);  %absolute frame
        xg=(handles.x(handles.track==trackid)-handles.offset)/handles.scale;
        yg=(handles.y(handles.track==trackid)-handles.offset)/handles.scale;
        %xg=handles.x(handles.track==trackid)/handles.scale;
        %yg=handles.y(handles.track==trackid)/handles.scale;
        inte=handles.intensity(handles.track==trackid);
        %clean up (make inte values NaN if they are > 10x bigger than the
        %median value) or less than zero
        %inte(inte>handles.foldcutoff*median(inte))=NaN;
        %inte(inte<0)=NaN;
        %figure(1); clf;  hold on;
        %plot(inte, '-*', 'Color', 'b');   plot(xg*1000,'-o','Color','r'); 
        %plot(yg*1000,'-o','Color', 'g');



        fprintf(fout, '%0.f\t%0.f\t%f\t%f\t%f\n', [rfr; afr; xg; yg; inte]);
        fclose(fout);
        j=j+1;
end



%% Gaussian fitting

function [pos, phot_count, a, normChi2,eccentricity]= getPosGauss2d_ed(im, posGuess,  windowSize)
%06092016 modified by Manu

% simple 2d gaussian fitting
% careful with (i,j) vs (x,y)!!

% setup some limits
posThresh = NaN; %limits turned off in freePosEllipseGaussFit_matlab
initguess = []; %no initial guess (ie freeGaussFitEllipse to calculate it
sigmaLim = [NaN NaN];% limits turned off in freePosEllipseGaussFit_matlab

%do the fit
% a: (A, sigma_x, sigma_y, b, Xo, Yo, theta )
% sigma_y is always the width along the major axis and theta is angle from y axis to the major axis
% phot_count = 2*pi*stdX*stdY*Amplitude (intensity over PSF volume)
[phot_count, a, normChi2, pos, eccentricity] = freeGaussFitEllipse( im, posGuess, windowSize,posThresh, sigmaLim, initguess);

%-------------------------------------------------------------------------------------------
function [phot_count,a,normChi2,pos,eccentricity ] = freeGaussFitEllipse( im, point_pos, windowSize,posLim, sigmaLim, initguess,varargin)
% function [phot_count pos normChi2 eccentricity a ] = freeGaussFitEllipse( im, point_pos, windowSize,posLim, sigmaLim, initguess)
% fit single 2d gaussian with fixed x, y position. wrapper to gauss fit tools function gaussfit_free_elliptical.cpp
%
% Inputs:
%   im - fit image should be of type double
%   point_pos - initial estimate of position
%   windowSize  - size of subimage to crop
%   posLim    - radius to allow shift from initial fit position
%   sigmaLim - [minwidth maxwidth]  - fit limits of psf width
%   initguess - [amplitude widthguess background X_POSim Y_POSim ];
% Outputs:
%   phot_count - volume of psf
%   a -  fit parameters, ie  [A, sigma_x, sigma_y, b, X, Y, theta ];
%    sigma_y is always the width along the major axis and theta is angle from y axis to the major axis
%
% NB ONLY ALLOWS SQUARE SUBIMAGES OTHERWISE RETURN 0
%
% Elliptical gaussian from cpp function:
%    xprime = (X-Xo)*cos(theta) - (Y-Yo)*sin(theta);
%    yprime = (X-Xo)*sin(theta) + (Y-Yo)*cos(theta);
%    e = exp((-(pow(xprime/xdenom,2)))-(pow(yprime/ydenom,2)));
%
%    x[element] = (A * e) + b;
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only�? license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
n = numel(varargin);
i = 1;
%defaults:
useCPPfit = false;
useAutoDetTol = false;
while i <= n
    if strcmp(varargin{i},'autoDetect')
        useAutoDetTol = true;
        i = i+1;
    elseif strcmp(varargin{i},'useMatlabFit')
        useCPPfit = false;
        i = i+1;
    else
        error('Unrecognised argument');
    end
end

if useAutoDetTol == false
    % default convergence tolerances
    verbose = false;
    epsilon1 =  10^-7; %gradient on lsq -
    epsilon2 =  10^-9; %gradient on fitParam - the most important
    epsilon3 = 0;  %absoluteValueLSQ - problem dependent - usually ~10^3 - NEVER USE IT!
    maxIter  = 100; %how fast it is in the absence of signal
else
    % for autodetection, convergence tolerances can be reduced (sub nanometre localization not required ofr accurate photon count
    % default convergence tolerances
    verbose = false;
    epsilon1 =  10^-2; %gradient on lsq -
    epsilon2 =  10^-2; %gradient on fitParam - the most important
    epsilon3 = 0;  %absoluteValueLSQ - problem dependent - usually ~10^3 - NEVER USE IT!
    maxIter  = 20; %how fast it is in the absence of signal
end

[sizey, sizex] = size(im);
X0=point_pos(1);
Y0=point_pos(2);

%round X0, Y0 to use as matrix locations
X0_int = round(X0);
Y0_int = round(Y0);
windowSize = round(windowSize); %radius should already be an integer anyway

% setup the limits of the cropped image
xstart =  max(1,X0_int-windowSize);
xfinish = min(sizex,X0_int+windowSize);
ystart =  max(1,Y0_int-windowSize);
yfinish = min(sizey,Y0_int+windowSize);

%crop to a small area around the point
fitIm = im( ystart:yfinish, xstart:xfinish);
fitIm = double(fitIm);%in case it's not already double
[sizeyFit,sizexFit] = size(fitIm);
% set up the point location in the cropped image coords
X_POSim = X0-xstart+1;
Y_POSim = Y0-ystart+1;

%set up the XY Lims
xLim = [(X_POSim - posLim) (X_POSim + posLim)];
yLim = [(Y_POSim - posLim) (Y_POSim + posLim)];
minwidth = sigmaLim(1);
maxwidth = sigmaLim(2);
%if an intial guess is not supplied, calculate it
if all(initguess==0)
    background = min(fitIm(:));
    amplitude = max(fitIm(:));
    widthguess = widthEstimate(fitIm)/2;
    if widthguess < minwidth
        widthguess = minwidth;
    elseif widthguess > maxwidth
        widthguess = maxwidth;
    end

    initguess = [amplitude widthguess background X_POSim Y_POSim ];
end
%do the fit
curvefitoptions = optimset( 'lsqcurvefit');
curvefitoptions = optimset( curvefitoptions,'Jacobian' ,'on','Display', 'off',  'TolX', epsilon2, 'TolFun', epsilon1,'MaxPCGIter',1,'MaxIter',maxIter);
[a, normChi2] =freePosEllipseGaussFit_matlab(fitIm,initguess ,xLim,yLim, sigmaLim,curvefitoptions);

% a: (A, sigma_x, sigma_y, b, Xo, Yo, theta )
I0 = a(1);
stdX = a(2);
stdY = a(3);
BG0 = a(4);
pos = [(a(5)+xstart-1) (a(6)+ystart-1)];
theta = a(7);
eccentricity = sqrt( 1 - (min(stdX,stdY)/max(stdX,stdY) )^2);
phot_count = 2*pi*stdX*stdY*I0;

%----------------------------------------------------------------------------------------------
function [fitParam, normChi2] = freePosEllipseGaussFit_matlab(inputIm,initguess ,xLim,yLim, sigmaLim,curvefitoptions);
% fits point spread function,
% F = (fitParam(1)*exp(    -(xprime).^2/(2*fitParam(2)^2)+(yprime).^2) /(2*fitParam(3)^2)   ) + fitParam(4))
%
%           xprime = (X-Xo)*cos(theta) - (Y-Yo)*sin(theta);
%           yprime = (X-Xo)*sin(theta) + (Y-Yo)*cos(theta);
%
% extra fit params fitParam(7) = theta, X0 and Y0 are fitParam(5) and (6)
%

A0start = initguess(1);
BGstart = initguess(3);
widthStart = initguess(2);
xStart = initguess(4);
yStart = initguess(5);
xMin = xLim(1);
xMax = xLim(2);
yMin = yLim(1);
yMax = yLim(2);
sigmaMin = sigmaLim(1);
sigmaMax = sigmaLim(2);

%set up the mesh, size of the input image for use in fitting
[sizey sizex] = size(inputIm);
num_pixels = sizey*sizex;
[X,Y]= meshgrid(1:sizex,1:sizey);
grid = [X Y];

AMPSCALEFACTOR =max(inputIm(:))/100;
if AMPSCALEFACTOR <= 0
    AMPSCALEFACTOR = 1;
end
%rescale the variables - to make magnitude of amplitude, background and width similar
inputIm = inputIm./AMPSCALEFACTOR;

% initguess input is  [amplitude widthguess background X_POSim Y_POSim ]
% initGuess7Vector output is [amplitude sx sy background X_POSim Y_POSim theta]
A0start = A0start./AMPSCALEFACTOR; %amplitude
BGstart = BGstart./AMPSCALEFACTOR; %backgound
if ( (initguess(2) < sigmaMin) || (initguess(2) > sigmaMax))%if given an out of bounds generate a
    initguess(2) = (sigmaMax+sigmaMin)/2; %sensible one- careful this isnt too close to true val tho
end
thetaStart = 0;
initGuess7Vector = [A0start widthStart widthStart BGstart xStart yStart thetaStart];
% A, sigma_x, sigma_y, b, Xo, Yo, theta

% Set fit limits on [amplitude widthx widthy background theta]
% dont set limits on theta but convert it to range 0->2pi afterwards
%MODIFIED to get rid of limits for SCCF calculation;110317SH
%lb = [0 sigmaMin sigmaMin 0 xMin yMin -inf];
%ub = [65535 sigmaMax sigmaMax 65535  xMax yMax inf];
lb = [];
ub = [];
%keyboard
%do the fit
try
    [fitParam, res] = ...
        lsqcurvefit(@(x, xdata) gauss2dw(x, xdata), ...
        initGuess7Vector ,grid ,inputIm ,...
        lb,ub,curvefitoptions);
catch ME
    if strcmp(ME.identifier,'optim:snls:InvalidUserFunction') % supplied absolutely empty image!
        fitParam = [0 0 0 0 0 0 0];
        res = 0;
    else
        rethrow(ME);
    end
end

%MODIFIED FOR SCCF CALCULATION! 110317SH
% modify fitParam to fit with "fitParam" output syntax from twotoneMain
% fitParam = (A, sigma_x, sigma_y, b, Xo, Yo, theta )
%if fitParam(1)< 0 %We know that negative amplitude values are patently unphysical so ignore them
%  fitParam(1) = 0;
%end

fitParam(1) = fitParam(1).*AMPSCALEFACTOR;%amplitude
fitParam(4) = fitParam(4).*AMPSCALEFACTOR;%background

fitParam(7) = mod(fitParam(7),2*pi);% transform theta onto the interval 0->2*pi
normChi2 = res/num_pixels;

%-------------------------------------------------------------------------
%---------------------Fitting Subfunctions--------------------------------
%-------------------------------------------------------------------------

function [F J] = gauss2dw(a, data)
% Used by the curve fitter to calculate values for a 2d gaussian
% with the x & y standard deviations equal
% and with fixed positions
% a(1) - A0
% a(2) - sX
% a(3) - sY
% a(4) - B
% a(5) - Xpos
% a(6) - Ypos
% a(7) - theta

%Initialise everything
[sizey sizex] = size(data);
sizex = sizex/2;

F = zeros(sizey, sizex);
X = F;
Y = F;

X = data(:,1:sizex);
Y = data(:,sizex+1:end);

xprime = (X-a(5))*cos(a(7)) - (Y-a(6))*sin(a(7));
yprime = (X-a(5))*sin(a(7)) + (Y-a(6))*cos(a(7));

% Only evaluate the exponential once:
expPart = exp( - ((xprime).^2 /(2*a(2)^2) + (yprime).^2 /(2*a(3)^2) ));

F = a(1)*expPart + a(4);

% compute the jacobian

% initialise everything
n = numel(F);
J = zeros(n,3); % initialise J
Ga1F = zeros(sizey, sizex);% dF/da(1)
Ga2F = Ga1F;% dF/da(2)
Ga3F = Ga1F;% dF/da(3)
Ga4F = Ga1F;% dF/da(4)
Ga5F = Ga1F;% dF/da(7)

% Calculate the grad_a1(F),  grad_a2(F), etc

Ga1F = expPart;

Ga2F = a(1).* expPart .*xprime.^2 .*a(2).^-3;% (A * e) * (pow(xprime,2) * pow(sigma_x,-3)); //dF/dsigma_x
Ga3F = a(1).* expPart .*yprime.^2 .*a(3).^-3;% (A * e) * (pow(yprime,2) * pow(sigma_y,-3)); //dF/dsigma_y

Ga4F = ones(size(X));
%dF/dX0 and dF/dY0 in cpp notation
%   jac[j++] = (A * e) * ( (xprime*cos(theta)*pow(sigma_x,-2)) + (yprime*sin(theta)*pow(sigma_y,-2)) ); //dF/dXo
%           jac[j++] = (A * e) * ( (yprime*cos(theta)*pow(sigma_y,-2)) - (xprime*sin(theta)*pow(sigma_x,-2)) ); //dF/dYo
Ga5F = a(1).* expPart .* ...
    (  xprime.*a(2).^(-2).*cos(a(7)) + yprime.*a(3).^(-2)*sin(a(7)) );
Ga6F = a(1).* expPart .* ...
    ( -xprime.*a(2).^(-2).*sin(a(7)) + yprime.*a(3).^(-2)*cos(a(7)) );

%dF/da(7) in c++ notation:
% (-A * e) *( (-xprime * pow(sigma_x,-2)) * ((X-Xo)*sin(theta) + (Y-Yo)*cos(theta)) + (yprime*pow(sigma_y,-2))*((X-Xo)*cos(theta) - (Y-Yo)*sin(theta)) );
Ga7F = -a(1).* expPart.* ...
    (  (-xprime).*a(2).^(-2).*((X-a(5))*sin(a(7))+ (Y-a(6))*cos(a(7))) ...
    + (yprime).*a(3).^(-2).*((X-a(5))*cos(a(7))- (Y-a(6))*sin(a(7))) );

% Form the jacobian, see the printed notes on getGaussFit for derivation
J = [Ga1F(:) Ga2F(:) Ga3F(:) Ga4F(:) Ga5F(:) Ga6F(:) Ga7F(:)];

%----------------------------------------------------------------------------------------------
function widthEst = widthEstimate(m)
%function to better estimate width for initial guess
[sizey sizex] = size(m);
vx = sum(m);
vy = sum(m');

vx = vx.*(vx>0);
vy = vy.*(vy>0);

x = [1:sizex];
y = [1:sizey];

cx = sum(vx.*x)/sum(vx);
cy = sum(vy.*y)/sum(vy);

sx = sqrt(sum(vx.*(abs(x-cx).^2))/sum(vx));
sy = sqrt(sum(vy.*(abs(y-cy).^2))/sum(vy));

widthEst = 0.5*(sx + sy) ;

%%Older webbfit gaussian  PSFwidth=round(wavelength/2xNA)pixelsize(nm)=2(!)
function [fitpos N] = webbfit(im,PSFwidth)
% 
% A = size(im);
% HalfWidth = (A(1) - 1) ./ 2;
% x = -HalfWidth:HalfWidth;
% [X Y] = ndgrid(x, x);
[X,Y] = ndgrid(1:size(im,1), 1:size(im,2));

% M = size(X);

fitpos = [0 0];
for i = 1:200
    oldpos = fitpos;
    % compute the mask
    gaussmask = exp(-0.5*((X-fitpos(1)).^2)/(PSFwidth) -0.5*((Y-fitpos(2)).^2)/(PSFwidth));
    denom = sum(sum(im.*gaussmask));
    fitpos(1) = sum(sum(X.*im.*gaussmask))/denom;
    fitpos(2) = sum(sum(Y.*im.*gaussmask))/denom; 
end
N = sum(sum(denom))/sum(sum(gaussmask.*gaussmask));
% ChiSq = sum(sum((MolImg - NumSigPhot .*
% squeeze(MaskArray(MaskElements(1), MaskElements(2), :, :))).^2));

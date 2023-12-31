function Bleach3_StepFindIntensity
%This is a reworking of Elizabeth Villa's stepfinder.  It reads in the
%output from Alex2Track and searches for steps in intensity measurement.  It outputs a new file starting all_fit_ containing an additional two
%columns with the average and STDev of intensity for the different
%plateaus.

%It should be run in the directory in which the output files from
%Alex2Track are located.  It will automatically run on all of them.

%The following are the default values needed for the stepFinder script.

minStep = 15;        % minimum accepted jump (nm)
W = 8;              % window size
pass = 10;          % number of passes of the filter
usemedian=1;     % set to 0 to calculate the mean values for plateaus


files = dir('*.txt');
%Strip out those starting with all_fit
idx = cellfun('isempty',strfind({files.name},'all_fit'));
files=files(idx);
for i=1:length(files)
    try
    stepFinder(files(i).name,minStep,W,pass);
    catch
        disp('PROCESSING STOPPED - FILE FORMAT INCORRECT')
    end
end

    function stepFinder(filename,minStep,W,pass)
        %  Fits steps to a set of data points from single molecule processive motors.
        %
        %     usage:  stepFinder(filename,minStep,W,pass,proj,pixelsize)
        %
        %  where
        %   filename: is the name of a file containing x and y values of the
        %  data (that should previously been rotated along the axis of the
        %  movement and cleaned from noise);
        %    minStep: is the minimum step the motor is expected to take
        %    W:       window size for the averaging. Should be chosen so at most
        %  one step is taken during this window
        %    pass:    how many times (passes) the filter should run the data
        %  through (needs to be automated for convergence).
        %    proj:    which coordinate to look into (values: 1 or 2). In the case
        %  of data passed by alignTrace, the axis along the axoneme is 1.
        %    pixelsize:  calibration measurement for the microscope.
        %    timestep:   the time between frames in seconds (e.g., 0.1 for 100ms)
        %
        %  outputs three files:
        %     fit_filename    the fitted trace
        %     step_filename   histogram of steps (useful for later dealing with
        % multiple data sets by just adding these files directly)
        %     dwell_filename  histogram of dwell times
        %
        %  Notes: A fit to a double exponential is included, but I am not sure how
        %         well it is working.
        %         Need to include the timestamps directly from the metadata, or at
        %         least, read it as a parameter (right now is hard coded :-P)
        %
        % stepFinder v.0.1
        % Elizabeth Villa, Physiology Course 2007
        % villa@ks.uiuc.edu
        
        
        fprintf(1,'****** Filename %s ******\n',filename);
        
        % Clean the data a little. At this stage it only gets rids of values that
        % are very far from the trace, so that they are not counted in the PCA
        % analysis to align the trace, or as steps.
        % Needs to be upgraded, perhaps combined with filterFluct.
        
        
            %get intensity data from file
            %filename='Diff#HaloACT_220216_1002@0.txt';
            A= textread(filename);
            %remove Nans
            B = A(any(~isnan(A(:,5)),2),:); 
            fm=B(:,2);
            I=B(:,5);

        
        % The filter, based on Smith, D.A.,  Phil. Trans. R. Soc. Lond. B (1998) 353,
        % 1969-1981.
        [I0, I, Px] = filterX(I,W,pass);  %I0 is the original data, I is the filtered trace, Px are spikes corresponding to steps in filtered trace
        
        %  Complicated filter
        
        NT = length(Px);
        halfE = round(W/2);
        
        % Threshold for taking a step
        % Note the value of sigma should be close to one coming out of the filter
        % (Needs to be automatically looked into and throw a warning if not)
        
        %Pthresh = (2*minStep/3*sigma);
        Pthresh = (2*minStep/3);
        
        %%%%%%%%%%%%% X %%%%%%%%%%%%%
        
        minpx = findminima(Px);
        maxpx = findmaxima(Px);
        
        % Arrays of maxima and minima
        Pmax=zeros(1,length(Px));
        Pmax(maxpx) = Px(maxpx);
        Pmin=zeros(1,length(Px));
        Pmin(minpx) = Px(minpx);
        
        % Clean so that there are no two peaks inside a window of width halfE*2
        for i=halfE*2:length(Pmax)-halfE*2
            if (Pmax(i) < max(Pmax(i:i+halfE)) || (Pmax(i) < max(Pmax(i-halfE:i))))
                Pmax(i) = 0;
            end
            if (Pmin(i) > min(Pmin(i:i+halfE)) || (Pmin(i) > min(Pmin(i-halfE:i))))
                Pmin(i) = 0;
            end
            if 1
                if Pmax(i) == Pmax(i+1)
                    Pmax(i+1) = 0;
                end
                if Pmin(i) == Pmin(i+1)
                    Pmin(i+1) = 0;
                end
            end
        end
        
        % Get the filtered stepping probability
        PedgeX = zeros(1,length(Px));
        PedgeX = Pmin + Pmax;
        
        % Find the intervals between jumps
        intervalsX =  [1 find(abs(PedgeX)>Pthresh) length(PedgeX)];
        %size(intervalsX);
        
        % The actual fit
        %
        if usemedian==0
            for i=1:length(intervalsX)-1
                if i == 1
                    xfit(intervalsX(i):intervalsX(i+1)) = mean(I0(intervalsX(i)+1:intervalsX(i+1)));
                    sfit(intervalsX(i):intervalsX(i+1)) = std(I0(intervalsX(i)+1:intervalsX(i+1)));
                else
                    xfit(intervalsX(i)+1:intervalsX(i+1)) = mean(I0(intervalsX(i)+1:intervalsX(i+1)));
                    sfit(intervalsX(i)+1:intervalsX(i+1)) = std(I0(intervalsX(i)+1:intervalsX(i+1)));
                end
            end
        elseif usemedian==1
            for i=1:length(intervalsX)-1
                if i == 1
                    xfit(intervalsX(i):intervalsX(i+1)) = median(I0(intervalsX(i)+1:intervalsX(i+1)));
                    sfit(intervalsX(i):intervalsX(i+1)) = std(I0(intervalsX(i)+1:intervalsX(i+1)));
                else
                    xfit(intervalsX(i)+1:intervalsX(i+1)) = median(I0(intervalsX(i)+1:intervalsX(i+1)));
                    sfit(intervalsX(i)+1:intervalsX(i+1)) = std(I0(intervalsX(i)+1:intervalsX(i+1)));
                end
            end
        else
            disp('Set usemedian to 0 or 1')
        end
        
      
        
        % Show the fit of x
        %figure(1); clf;  hold on;
        %plot(I0, '-*', 'Color', 'k');   plot(xfit,'-o','Color','r'); 
        %plot(sfit,'-o','Color', 'b');

        % Make a column containing NaN and fill relevant parts with xfit
        A(:,6:7)=NaN;
        %Fill those rows that have an xfit,sfit value with those values
        %(leave rest as NaN).  Recreates the original 
        for j=1:length(fm)
            idx=find(A(:,2)==fm(j));
            A(idx,6:7)=[xfit(j),sfit(j)];
        end
        
        tk=A(:,1)';
        frm=A(:,2)';
        x=A(:,3)';
        y=A(:,4)';
        I=A(:,5)';
        Ifit=A(:,6)';
        Sfit=A(:,7)';
        
        
        % All the fitted data to a file
        
        foutname = ['all_fit_' filename];
        fout =fopen(foutname,'w');
        fprintf(fout, '%0.f\t%0.f\t%f\t%f\t%f\t%f\t%f\n', [tk; frm; x; y; I; Ifit; Sfit]);
        fclose(fout);
        
        function[x0, x, pp, s1] = filterX(x0,W,pass)
            %  Filters the data in filename to sharpen the edges
            %  Implementation of Smith, D.E., XXX
            %     usage:  filterX(filename,minStep,W,pass,proj);
            %  where    filename    the name of the txt file,
            %           minStep     the minimal step the motor can make
            %           W           the window size
            %           pass        number of passes of the filter (needs automation!)
            %           proj        column of the axis along the axoneme (1,2)
            
            
     
            r = 10;             % Coefficient to which elevate the variances for the filter
            halfE = round(W/2); % Window in which no two maxima/minima can coexist
            NT = length(x0);
            % figure(1);
            %plot(x0, '-o', 'Color', 'b'); hold on;
            %set(gca,'YGrid','on','YTick', [0:minStep:2*median(x0)], 'XGrid', 'on','XTick',[0:20:length(x0)])
            
            for p=1:pass
                
                fprintf(1, '..%i..',p);
                if p == 1
                    x = x0;
                end
                
                % Extend time series by W points at each end
                % xpad = padarray(x,W);
                
                xpad = zeros(1,NT+2*W);
                xpad(1:W) = x(1);
                xpad(W+1:NT+W)= x(1:NT);
                xpad((NT+W+1):(NT+2*W+1)) = x(NT);
                
                
                % Calculate forward and backward running averages and standard deviations
                for i=1+W:NT+1+W
                    xavefor(i) = mean(xpad(i:i+W));
                    xvarfor(i) =  var(xpad(i:i+W));
                    xavebak(i) = mean(xpad(i-W:i));
                    xvarbak(i) =  var(xpad(i-W:i));
                end
                
                
                
                % The
                rsp = xvarfor.^r;
                rsm = xvarbak.^r;
                
                gm = rsp./(rsp+rsm);
                gp = rsm./(rsp+rsm);
                
                % Calculate the standard deviation of the trace
                if p == 1
                    s = zeros(1,length(gp));
                    s1 = sqrt(gp.*xvarfor+gm./xvarbak);  % signal from noise
                    for i=1:length(s)-3
                        save(i) = mean(s(i:i+3));
                    end
                    crap = find(save>2);
                end
                
                
                % The filtered signal!
                x = gp.*xavefor + gm.*xavebak;
                x = x(W+1:NT+W);
                x(crap) = 0;
                
                % figure(2); plot(x,'Color', [p/pass 0 0]); hold on;
                
                % In the last pass get the vvalue of sigma (the average noise signal)
                if p == pass
                    fprintf(1,'done\n');
                    s = sqrt(gp.*xvarfor+gm./xvarbak);  % signal from noise
                    sigma = mean(isfinite(s));          % average of signal from noise
                    if sigma < 1.1
                        pppad = (xavefor - xavebak);
                        pp = pppad(1+W:NT+W);               % unpadded Y
                        fprintf(1, 'Sigma is %f\n', sigma);
                    else
                        fprintf(2,'Your noise is still high (std= %f), consider filtering for longer times\n', sigma);
                    end
                    %plotting stuff
                    %close all
                    %figure(1); plot(x0); hold on, plot(x), 
                    %figure(2); plot(pp)
                    %figure(3); plot(s1)
                    
                end
            end
            fprintf('Filtering done\n');
        end 
        function minima = findminima(x)
            %FINDMAXIMA  Find location of local maxima
            %  From David Sampson
            %  See also FINDMINIMA
            
            % Unwrap to vector
            x = x(:);
            % Identify whether signal is rising or falling
            upordown = sign(diff(x));
            % Find points where signal is rising before, falling after
            minflags = [upordown(1)>0; diff(upordown)>0; upordown(end)<0];
            minima   = find(minflags);
        end
        function maxima = findmaxima(x)
            %FINDMAXIMA  Find location of local maxima
            %  From David Sampson
            %  See also FINDMINIMA
            
            % Unwrap to vector
            x = x(:);
            % Identify whether signal is rising or falling
            upordown = sign(diff(x));
            % Find points where signal is rising before, falling after
            maxflags = [upordown(1)<0; diff(upordown)<0; upordown(end)>0];
            maxima   = find(maxflags);
        end
    end







end



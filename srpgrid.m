function [finalpos,finalsrp]=srpgrid(flens, mic_loc, fs, lsb, usb)
%% This function uses SRP-PHAT with  grid search to locate a single source using a frame of data and 4 microphones.
%% Inputs: 
%%% 1) s is "A FRAME" of data (L x 4), L should be a power of 2
%%% 2) mic_loc is the microphone 3D-locations (4 x 3) ( in meters)
%%% 3) fs: sampling rate (Hz)
%%% 4) lsb: a row-vector of the lower rectangular search boundary (meters)
%%% 5) usb: row-vector of the upper rectangular search boundary (m)
%%% It also calls other subroutine: fe.
%% Outputs:
%%% 1) finalpos: estimated location of the source 
%%% 2) finalsrp: srp-phat value of the point source
%%% 3) finalfe: number of fe's evaluated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Usage example:
%%% [finalpos,finalsrp,finalfe]=srplems(s, mic_loc, fs, lsb, usb)
%%% 4 microphones ---> mic_loc is a (4x3) matrix.
%%% A 2048-sample framelength ---> s is a (2048 x 24) matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off all

if nargin < 5, usb=[10 0 3]; end
if nargin < 4, lsb=[-1 -1 0]; end
if nargin < 3, fs=48000; end


L = size(flens,1); %%% determine frame-length
M = size(mic_loc,1); %%% number of microphones (also number of channels of the input X).
np=M*(M-1)/2; %%%number of independent pairs

%steplength = L/4;  %%% non-overlapped length of a frame (75% overlap)
dftsize = L;       %%% dft size equals to frame size
temperatureC=24.0;
speedofsound=331.4*sqrt(1.0+(temperatureC/273));
magiconst=10*fs/speedofsound;  
interval = 1;
%% Determine the maximum end-fire length (in samples) of a microphone pairs:
mdist=pdist(mic_loc);
efmax=max(mdist);
efs=2*(fix(efmax*fs/speedofsound)); %%% End-fire is symmetric about the 0th sample so multiplying by 2.
%efs=801;
hefs=round(efs/2); %%%half of the end-fire region
%% Get the linear-indices of 'np' independent mic-pairs:
% w=1:M;
% wn=[0:M-1]'*M;
% fm1=repmat(wn,1,M);
% fm2=repmat(w,M,1);
% mm=fm1+fm2;
% tr=triu(mm,1);%%Get the upper half above the main diagonal.
% gidM=nonzeros(tr'); %%%keep only non-zero values -> linear indices of the pairs
% clear fm1 fm2 w wn

%% Doing the GCC-PHAT:

sf=fft(flens,dftsize);                    %%%FFT of the original signals  
csf=conj(sf);
yv=zeros(np,efs);                     %%%%Initialize yv to store the SRP-PHAT samples

p=1;

for k=1:M-1
      su1mic=sf(:,k)*ones(1,M);      %%%Create M copies of each signal's FFT
      prodall=su1mic.*csf;           %%%%Calculate the cross-power spectrum: = fft(x1).*conj(fft(x2))
      ss=prodall(:,k+1:M);           %%%% ss will be the cross-power spectra of microphone pairs (i,i+1), (i,i+2)...(i,M)
      ss=ss./(abs(ss)+eps);          %%%% PHAT weighting
      
      ssifft=real(ifft(ss,dftsize)); %%%% Get the GCC-PHAT ssifft, which is the IFFT of the cross-power spectra
      %newssifft=[ssifft(end-hefs:end,:); ssifft(1:efs-hefs-1,:)]; %%% Only select 'efs' samples (the beginning+end portions)
      newssifft=[ssifft(end-hefs+1:end,:); ssifft(1:efs-hefs,:)];
      newssifftr = newssifft';          
      yv(p:p+M-k-1,:)=newssifftr;    %%%%Store in yv
      p=p+M-k;                       %%%% Update the current index of yv
      clear su1mic prodall ss ssifft newssifft newssifftr
end

%% cubic interpolation (factor of 10) on GCC

xx=1:.1:efs;
x=1:efs;
yintp=spline(x,yv,xx);
yintpt=yintp';


efsintp=length(xx)/2;

row1=([0:np-1]*2*efsintp)'; 

[bestpos,bestsrp]=grid(lsb,usb,interval,magiconst,mic_loc,yintpt,row1,efsintp,3);

finalpos = bestpos;  %%%Final source location estimate
finalsrp = bestsrp/np;  %%%Final source's SRP-PHAT value (normalized by number of pairs)


    function [bestpos,bestsrp] = grid(bstart,bend,interval,magiconst,mic_loc,yintpt,row1,efsintp,dimen)
        %%%This function does grid search algorithm.
        bstart
        bend
        gridgen = round(abs(bstart-bend)./interval);
        num = prod(gridgen);
        pointvec = zeros(num,dimen);
        yval = zeros(1,num);
        for i= 1:gridgen(1)
            for t = 1:gridgen(2)
                for n = 1:gridgen(3)
                    point = [bstart(1)+i*interval bstart(2)+t*interval bstart(3)+n*interval]
                    pointvec(i,:) = point;
                    yval(i) = srpscore(point,magiconst,mic_loc,yintpt,row1,efsintp);
                end
            end
        end

        [v,ind]=max(yval);
        
        bestsrp=v;
        bestpos=pointvec(ind,:);
       
        [az,el,r] = cart2sph(pointvec(:,1),pointvec(:,2),pointvec(:,3));
        
        surf([pointvec(:,1),pointvec(:,2),yval']);
        %surf([pointvec(:,1),pointvec(:,2),yval']);
        %clear searchVol srp_vec
    end

%% calculate SRP-PHAT score for a point:
    
    function inc = increase(gridgen,dimen)
        temp = gridgen;
        inc = zeros(1,3);
        for j=1:dimen
            inc(j)= round(temp/dimen);
            temp = temp - inc(j);
        end
    end
    
        
    function yval1 = srpscore(point,magiconst,mic_loc,yintpt,row1,efsintp)
        %%%This function evaluates its SRP-PHAT value.
        M=size(mic_loc,1); %%%number of mics

        %%%generate a random point in the search space:
        x=point;   %This vector x defines the coordinates (x,y,z) of the rand point x
        %%%Find the distances from M-microphones to the point:
        a1=ones(M,1);
        xx1=a1*x;
        xdiff1=xx1-mic_loc;
        dists=sqrt(sum(xdiff1.*xdiff1,2));


        %%%%Differences in distances:
        ddiffs_ones=ones(M,1)*dists';
        ddm=ddiffs_ones-ddiffs_ones';

        %%% Calculate the TDOA index:
      %  v=ddm(gidM);
        v=nonzeros(tril(ddm,0));

        %ddiffsi32=int32(round(magiconst*v+efsintp));
        %ddiffsi32=floor(magiconst*v+efsintp)+1; 
        ddiffsi32=round(magiconst*v+efsintp+1.5); 
        row=row1+ddiffsi32;   
        
        %%%Pull out the GCC-PHAT value corresponding to that TDOA:
        v1=yintpt(row);
        %%% SRP-PHAT value of the point:
        yval1=sum(v1);
    end
end

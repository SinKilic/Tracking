%% Loading of file and generation of frame array
% Load the file with the ScanR analysis data including the time-frame and 
%the x and y coordinates
%ScanOut_name = input('Enter the filename of the csv file with ScanR output: ', 's');
%ScanOut_t = readtable(strcat(ScanOut_name,'.txt'));
[filename,path] = uigetfile('*.txt*',  'All Files (*.txt*)');
ScanOut_t = readtable([path,filename]);
ScanOut = table2array(ScanOut_t);
% Automatically determine the number of rows and columns
r_Scn = size(ScanOut,1);
c_Scn = size(ScanOut,2);
% Input the columns in which the time, x and y coordinates are located and
% the number of extra columns required (Default 6)
c_time = input('Enter the column with the time: ');
c_x = input('Enter the column with the x-coordinate: ');
c_y = input('Enter the column with the y-coordinate: ');
xtr_c = 6;
% Automatically assign the slide number in the last frame
time_fr = ScanOut(r_Scn,c_time);
% Initiate arrays/extra arrays to mark the usage of rows in the different
% cycles
ScanOut(:, c_Scn+1:c_Scn+xtr_c) = zeros(r_Scn,xtr_c);
track = zeros(r_Scn,c_Scn+xtr_c);
% Determine the number of cells in the last frame by summing the number of
% frame with that number
cell_max = sum(ScanOut(:,c_time) == time_fr);
% Create a threedimensional array with 1st dim being the max number of rows
% per frame, the second dimension being the same number as in the track and
% ScanR output and the third dimension being the number of frames
frame_array = zeros(cell_max,c_Scn+xtr_c,time_fr);
% Cycle through the time frames, cell rows and ScanR output rows
for j = 1:time_fr
    for k = 1:cell_max
        for i = 1:r_Scn
            % If the table information is in the same frame and has not
            % been used previously assign, 
            if ScanOut(i,c_time) == j && ScanOut(i,c_Scn+1) ~=1 && frame_array(k,c_Scn+1,j) ~=1
                frame_array(k,1:c_Scn,j) = ScanOut(i,1:c_Scn);
                frame_array(k,c_Scn+1,j) = 1;
                ScanOut(i,c_Scn+1) = 1;
            end
        end
    end
end
%%
frame_array_backup = frame_array;
%% Formation of initial tracks based on distance constraint
tic
mxdist = 100;
TRLINE = 1;
TRKID = 1;
for m = 1:time_fr
    for s = 1:cell_max
        if frame_array(s,2,m) == 0 % If the row is empty, stop searching in it
            break
        elseif frame_array (s,c_Scn+2,m) ~=1 % Only search in rows not previously used
            curr_cell = frame_array(s,:,m); % Assign current cell
            frame_array(s,c_Scn+2,m) = 1; % Tag object that it was used
            track(TRLINE,:) = curr_cell; % Copy cell data into tracking array
            track(TRLINE,c_Scn+3) = TRKID; % Copy TRKID to column 
            TRLINE = TRLINE+1; % Move to next line in track array
        end
            for j = 2:time_fr
                for i = 1:cell_max
                    if frame_array(i,2,j) == 0 % Search only in non-empty rows
                        break
                    end
                    dist = sqrt((curr_cell(1,c_x)-frame_array(i,c_x,j))^2+(curr_cell(1,c_y)-frame_array(i,c_y,j))^2);
                    frame_array(i,c_Scn+4,j) = dist;
                end
                A = frame_array(:,c_Scn+4,j);
                B = find(A==min(A(A>0)));
                if size(B,1)>1
                    B(2:size(B,1)) = [];
                end
                if frame_array(B,c_Scn+4,j) < mxdist && ((curr_cell(1,c_time) == frame_array(B,c_time,j)-1)||(curr_cell(1,c_time) == frame_array(B,c_time,j)-2)) &&  frame_array(B,c_Scn+2,j) ~= 1
                    curr_cell = frame_array(B,:,j);
                    track(TRLINE,:) = curr_cell;
                    frame_array(B,c_Scn+2,j) = 1;
                    track(TRLINE,c_Scn+3) = TRKID;
                    TRLINE = TRLINE + 1;
                end
            end
       TRKID = TRKID + 1;
    end
end
track(:,[1 c_Scn+3]) = track(:,[c_Scn+3 1]);
for i=size(track,1):-1:1
    if track(i,c_time)==0
        track(i,:) = [];
    end
end
% Rename tracks in a consecutive fashion (Get index sorted without numbers)
TRKID = 1;
for i = 1:size(track,1)
     if i == 1
         track(i,size(track,2)) = TRKID;
     elseif i == size(track,1) && track(i-1,1)
         TRKID = TRKID + 1;
         track(i,size(track,2)) = TRKID;
     elseif track(i,1) ~= track(i-1,1)
         TRKID = TRKID + 1;
         track(i,size(track,2)) = TRKID;
     else
         track(i,size(track,2)) = TRKID;
     end
end
track(:,[1 size(track,2)]) = track(:,[size(track,2) 1]);
track(:,size(track,2)) = 0;
toc
%% Load pre-existing csv file
[filename,path] = uigetfile('*.csv*',  'All Files (*.csv*)');
track = dlmread([path,filename]);
c_time = input('Enter the column with the time: ');
c_x = input('Enter the column with the x-coordinate: ');
c_y = input('Enter the column with the y-coordinate: ');
%% Output coordinates for loading to imagej
Coord = array2table(track(:,[1 c_time c_x c_y]));
%% Output csv file
dlmwrite("GFP-MDC1_2_slimmed.csv",track);
%% Track renaming
TRKCURR = 134;
TRKNEW = 2030;
for i = 1:size(track,1)
    if track(i,1) == TRKCURR
        track(i,1) = TRKNEW;
    end
end
track = sortrows(track,[1 c_time]);
disp('Tracks renamed');
%% Exchange tracks from one point onward
TPEXCH = 4;
TRK1 = 29;
TRK2 = 540;
for i = 1:size(track,1)
    if track(i,1) == TRK1 && track(i,c_time) >= TPEXCH
        track(i,1) = TRK2;
    elseif track(i,1) == TRK2 && track(i,c_time) >= TPEXCH
        track(i,1) = TRK1;
    end
end
track = sortrows(track,[1 c_time]);
disp('Tracks exchanged');
%% Merge tracks from one time-point onward
TPMERGE = 36;
TRKOLD = 210;
TRKMERGE = 194;
for i = 1:size(track,1)
    if track(i,1) == TRKOLD && track(i,c_time) >= TPMERGE
        track(i,1) = TRKMERGE;
    end
end
track = sortrows(track,[1 c_time]);
disp('Tracks merged');
%% Remove track in time interval
TRKINTREM = 603;
TPREM1 = 91;
TPREM2 = 91;
for i=size(track,1):-1:1
    if track(i,1) == TRKINTREM && track(i,c_time) >= TPREM1 && track(i,c_time) <= TPREM2
        track(i,:) = [];
    end
end
track = sortrows(track,[1 c_time]);
disp('Track in interval removed');
%% Remove tracks completely
TRKRMV = 541;
for i =size(track,1):-1:1
    if track(i,1) == TRKRMV
        track(i,:) = [];
    end
end
track = sortrows(track,[1 c_time]);
disp('Track removed');
%% Remove unassigned tracks below 1000 and subtract 1000 from remaining
for i = size(track,1):-1:1
    if track(i,1) < 1000
        track(i,:) = [];
    end
end
for i = 1:size(track,1)
    track(i,1) = track(i,1)-1000;
end
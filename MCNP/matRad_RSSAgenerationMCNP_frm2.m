function [runfileName, WSSAfileName] = matRad_RSSAgenerationMCNP_frm2(MLCleafPositions, MLCidentifier, runfileName, SSD4kernelCalc, WSSAoption, npsKernelGen, particleType, transportMedium)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate dose for predefined MLC shape either from dose calculation in
% matRad or predefined.
%
% call
%   [control_MLCfieldWSSA, control_MLCfieldPhantom] =
%                                   kernelGenerationMCNP(resultGUI, beamOfInterest, MLCshapeOfInterest, WSSAoption)
%
% input
%   MLCleafPositions:   leaf positions have to be defined from left to
%                       right in a 2-by-#leafs sized matrix where upper row
%                       defines leaf bank A and lower row defines leaf bank B
%   MLCidentifier:      string identifier to select beam shaping device
%   runfileName:        name of runfile to be generated
%   WSSAoption:         boolean to control generation of WSSA file
%
% output
%   control_MLCfieldWSSA:   control variable to check if WSSA was created
%                           sucessfully
%
% References
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 04/2020
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check function input
if nargin < 2
    error('Input not defined properly. Please make sure to characterize leaf positions and the name of the predefined beam shaping device you would like to simualte.')
end


%% Get paths
pathRunfiles = pwd;
% pathWSSAfiles = strcat(matRad_getMATRADdirectory_hardCoded, filesep, 'MCNP', filesep, 'WSSA_depot', filesep);

%% Generate runfiles for given leaf positions
listMLCid = ["frm2MLC1"; "frm2MLC3"];   % Please list all names for known beam shaping devices here

if isempty(listMLCid==MLCidentifier); error('Beam shaping device unknown! Cannot calculate radiation field...'), end
switch MLCidentifier
    case 'frm2MLC1'
        disp('*****')
        disp('Calculation of radiation field performed for MLC1...')
        disp('*****')
        
        % Check if MLC positions extend maximum field width
        maxFieldSize_MLC1 = 20; % vertical field size in cm
        if max(max(abs(MLCleafPositions)))>maxFieldSize_MLC1/2
            error('Leaf positions out of range, maximum vertical field size of MLC1 exceeded!')
        end
        
        % Write runfiles for FRM2 MLC1 for given leafpositions
        runfileName = makeRunfileFRM2MLC1(MLCleafPositions, pathRunfiles, runfileName);
        
        % Write runfiles for water phantom PTW
        sourceOpeningAngle = 3.39;
        [waterPhantomRunfileName, WSSAfileName] = makeWaterPhantomFRM2PTW(runfileName, SSD4kernelCalc, npsKernelGen, particleType, transportMedium, sourceOpeningAngle, WSSAoption);
        
        
    case'frm2MLC3'
        disp('*****')
        disp('Calculation of radiation field performed for MLC3...')
        disp('*****')
        
        % Check if MLC positions extend maximum field width
        maxFieldSize_MLC3 = 18; % vertical field size in cm
        if max(max(abs(MLCleafPositions)))>maxFieldSize_MLC3/2
            error('Leaf positions out of range, maximum vertical field size of MLC3 exceeded!')
        end
        
        % Write runfiles for FRM2 MLC1 for given leafpositions
        runfileName = makeRunfileFRM2MLC3(MLCleafPositions, pathRunfiles, runfileName);
        
        % Write runfiles for water phantom PTW
        sourceOpeningAngle = 1.63;
        [waterPhantomRunfileName, WSSAfileName] = makeWaterPhantomFRM2PTW(runfileName, SSD4kernelCalc, npsKernelGen, particleType, transportMedium, sourceOpeningAngle, WSSAoption);
end

% Concatenate runfiles
runfileName = concatenateRunfile(runfileName, waterPhantomRunfileName);



% Clean up
fragmentList = dir('*block*');
for counter = 1:length(fragmentList)
    delete(fragmentList(counter).name)
end


%% Function definition for FRM 2 MLC 1 MCNP runfile goes here...
    function runfileName = makeRunfileFRM2MLC1(MLCleafPositions, pathRunfiles, runfileName)
        fileID_runfile_blockA = fopen(strcat(pathRunfiles,filesep,runfileName, '_blockA.txt'), 'w');
        fileID_runfile_blockB = fopen(strcat(pathRunfiles,filesep,runfileName, '_blockB.txt'), 'w');
        fileID_runfile_blockC = fopen(strcat(pathRunfiles,filesep,runfileName, '_blockC.txt'), 'w');
        
        % Write block A: Cells
        fprintf(fileID_runfile_blockA, 'MCNP Runfile for WSSA/matRad Kernel Generation for MLC1 at FRM2\n');
        fprintf(fileID_runfile_blockA, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockA, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockA, 'C Block A: Cells\n');
        fprintf(fileID_runfile_blockA, 'C ***************************************************************\n');
        
        fprintf(fileID_runfile_blockA, '1000 0 7000 9 (-1:-2:-3:-4:-5:-7:-8) $ Propagation tunnel\n');
        fprintf(fileID_runfile_blockA, '9 0 -9\n'); % RSSA helper cell
        fprintf(fileID_runfile_blockA, '1001 0  (-6 $ Start MLC $ Inside of MLC for fluence tally\n');
        fprintf(fileID_runfile_blockA, '        501 $ Bank A leaf 1\n');  % leaf bank A
        fprintf(fileID_runfile_blockA, '        502 $ Bank A leaf 2\n');
        fprintf(fileID_runfile_blockA, '        503 $ Bank A leaf 3\n');
        fprintf(fileID_runfile_blockA, '        504 $ Bank A leaf 4\n');
        fprintf(fileID_runfile_blockA, '        505 $ Bank A leaf 5\n');
        fprintf(fileID_runfile_blockA, '        506 $ Bank A leaf 6\n');
        fprintf(fileID_runfile_blockA, '        507 $ Bank A leaf 7\n');
        fprintf(fileID_runfile_blockA, '        508 $ Bank A leaf 8\n');
        fprintf(fileID_runfile_blockA, '        509 $ Bank A leaf 9\n');
        fprintf(fileID_runfile_blockA, '        510 $ Bank A leaf 10\n');
        fprintf(fileID_runfile_blockA, '        511 $ Bank A leaf 11\n');
        fprintf(fileID_runfile_blockA, '        512 $ Bank A leaf 12\n');
        fprintf(fileID_runfile_blockA, '        513 $ Bank A leaf 13\n');
        fprintf(fileID_runfile_blockA, '        514 $ Bank A leaf 14\n');
        fprintf(fileID_runfile_blockA, '        515 $ Bank A leaf 15\n');
        fprintf(fileID_runfile_blockA, '        516 $ Bank A leaf 16\n');
        fprintf(fileID_runfile_blockA, '        517 $ Bank A leaf 17\n');
        fprintf(fileID_runfile_blockA, '        518 $ Bank A leaf 18\n');
        fprintf(fileID_runfile_blockA, '        519 $ Bank A leaf 19\n');
        fprintf(fileID_runfile_blockA, '        520 $ Bank A leaf 20\n');
        
        fprintf(fileID_runfile_blockA, '        601 $ Bank B leaf 1\n');  % leaf bank B
        fprintf(fileID_runfile_blockA, '        602 $ Bank B leaf 2\n');
        fprintf(fileID_runfile_blockA, '        603 $ Bank B leaf 3\n');
        fprintf(fileID_runfile_blockA, '        604 $ Bank B leaf 4\n');
        fprintf(fileID_runfile_blockA, '        605 $ Bank B leaf 5\n');
        fprintf(fileID_runfile_blockA, '        606 $ Bank B leaf 6\n');
        fprintf(fileID_runfile_blockA, '        607 $ Bank B leaf 7\n');
        fprintf(fileID_runfile_blockA, '        608 $ Bank B leaf 8\n');
        fprintf(fileID_runfile_blockA, '        609 $ Bank B leaf 9\n');
        fprintf(fileID_runfile_blockA, '        610 $ Bank B leaf 10\n');
        fprintf(fileID_runfile_blockA, '        611 $ Bank B leaf 11\n');
        fprintf(fileID_runfile_blockA, '        612 $ Bank B leaf 12\n');
        fprintf(fileID_runfile_blockA, '        613 $ Bank B leaf 13\n');
        fprintf(fileID_runfile_blockA, '        614 $ Bank B leaf 14\n');
        fprintf(fileID_runfile_blockA, '        615 $ Bank B leaf 15\n');
        fprintf(fileID_runfile_blockA, '        616 $ Bank B leaf 16\n');
        fprintf(fileID_runfile_blockA, '        617 $ Bank B leaf 17\n');
        fprintf(fileID_runfile_blockA, '        618 $ Bank B leaf 18\n');
        fprintf(fileID_runfile_blockA, '        619 $ Bank B leaf 19\n');
        fprintf(fileID_runfile_blockA, '        620) $ Bank B leaf 20\n');
        
        fprintf(fileID_runfile_blockA, '9999 0 1 2 3 4 5 6 7 8 $ Graveyard\n');
        
        fprintf(fileID_runfile_blockA, '1002 0 -6 $ Exclude MLC\n');
        fprintf(fileID_runfile_blockA, '        (-501: $ Bank A leaf 1\n');  % leaf bank A
        fprintf(fileID_runfile_blockA, '        -502: $ Bank A leaf 2\n');
        fprintf(fileID_runfile_blockA, '        -503: $ Bank A leaf 3\n');
        fprintf(fileID_runfile_blockA, '        -504: $ Bank A leaf 4\n');
        fprintf(fileID_runfile_blockA, '        -505: $ Bank A leaf 5\n');
        fprintf(fileID_runfile_blockA, '        -506: $ Bank A leaf 6\n');
        fprintf(fileID_runfile_blockA, '        -507: $ Bank A leaf 7\n');
        fprintf(fileID_runfile_blockA, '        -508: $ Bank A leaf 8\n');
        fprintf(fileID_runfile_blockA, '        -509: $ Bank A leaf 9\n');
        fprintf(fileID_runfile_blockA, '        -510: $ Bank A leaf 10\n');
        fprintf(fileID_runfile_blockA, '        -511: $ Bank A leaf 11\n');
        fprintf(fileID_runfile_blockA, '        -512: $ Bank A leaf 12\n');
        fprintf(fileID_runfile_blockA, '        -513: $ Bank A leaf 13\n');
        fprintf(fileID_runfile_blockA, '        -514: $ Bank A leaf 14\n');
        fprintf(fileID_runfile_blockA, '        -515: $ Bank A leaf 15\n');
        fprintf(fileID_runfile_blockA, '        -516: $ Bank A leaf 16\n');
        fprintf(fileID_runfile_blockA, '        -517: $ Bank A leaf 17\n');
        fprintf(fileID_runfile_blockA, '        -518: $ Bank A leaf 18\n');
        fprintf(fileID_runfile_blockA, '        -519: $ Bank A leaf 19\n');
        fprintf(fileID_runfile_blockA, '        -520: $ Bank A leaf 20\n');
        
        fprintf(fileID_runfile_blockA, '        -601: $ Bank B leaf 1\n');  % leaf bank B
        fprintf(fileID_runfile_blockA, '        -602: $ Bank B leaf 2\n');
        fprintf(fileID_runfile_blockA, '        -603: $ Bank B leaf 3\n');
        fprintf(fileID_runfile_blockA, '        -604: $ Bank B leaf 4\n');
        fprintf(fileID_runfile_blockA, '        -605: $ Bank B leaf 5\n');
        fprintf(fileID_runfile_blockA, '        -606: $ Bank B leaf 6\n');
        fprintf(fileID_runfile_blockA, '        -607: $ Bank B leaf 7\n');
        fprintf(fileID_runfile_blockA, '        -608: $ Bank B leaf 8\n');
        fprintf(fileID_runfile_blockA, '        -609: $ Bank B leaf 9\n');
        fprintf(fileID_runfile_blockA, '        -610: $ Bank B leaf 10\n');
        fprintf(fileID_runfile_blockA, '        -611: $ Bank B leaf 11\n');
        fprintf(fileID_runfile_blockA, '        -612: $ Bank B leaf 12\n');
        fprintf(fileID_runfile_blockA, '        -613: $ Bank B leaf 13\n');
        fprintf(fileID_runfile_blockA, '        -614: $ Bank B leaf 14\n');
        fprintf(fileID_runfile_blockA, '        -615: $ Bank B leaf 15\n');
        fprintf(fileID_runfile_blockA, '        -616: $ Bank B leaf 16\n');
        fprintf(fileID_runfile_blockA, '        -617: $ Bank B leaf 17\n');
        fprintf(fileID_runfile_blockA, '        -618: $ Bank B leaf 18\n');
        fprintf(fileID_runfile_blockA, '        -619: $ Bank B leaf 19\n');
        fprintf(fileID_runfile_blockA, '        -620) $ Bank B leaf 20\n');
        
        fprintf(fileID_runfile_blockB, '\n');
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockB, 'C Block B: Surfaces\n');
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        
        fprintf(fileID_runfile_blockB, '1 RPP 0 62 -50 50 -50 50  $ Outer boundaries of simulation volume\n');
        fprintf(fileID_runfile_blockB, '4 RPP 288 418 -50 50 -50 50\n');
        fprintf(fileID_runfile_blockB, '8 RPP 493 700 -50 50 -50 50\n');
        fprintf(fileID_runfile_blockB, '9 RPP 492.5 492.9 -10.9 10.9 -14.9 14.9\n');  % RSSA generation helper surface
        
        fprintf(fileID_runfile_blockB, '2 ARB 62 -8 -9 $ Inner boundaries of transport tunnnel\n');
        fprintf(fileID_runfile_blockB, '        62 -8 9\n');
        fprintf(fileID_runfile_blockB, '        62 8 9\n');
        fprintf(fileID_runfile_blockB, '        62 8 -9\n');
        fprintf(fileID_runfile_blockB, '        274 -9 -12\n');
        fprintf(fileID_runfile_blockB, '        274 -9 12\n');
        fprintf(fileID_runfile_blockB, '        274 9 12\n');
        fprintf(fileID_runfile_blockB, '        274 9 -12\n');
        fprintf(fileID_runfile_blockB, '        1234 2367 1458 3478 1256 5678\n');
        fprintf(fileID_runfile_blockB, '3 RPP 274 288 -15 15 -17.5 17.5\n');
        fprintf(fileID_runfile_blockB, '5 RPP 418 428 -10 10 -15 15\n');
        fprintf(fileID_runfile_blockB, '6 RPP 428 478 -10 10 -15 15\n');
        fprintf(fileID_runfile_blockB, '7 RPP 478 493 -11 11 -15 15\n');
        
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockB, 'C MLC Leafs\n');
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockB, '501 501 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '502 502 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '503 503 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '504 504 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '505 505 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '506 506 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '507 507 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '508 508 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '509 509 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '510 510 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '511 511 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '512 512 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '513 513 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '514 514 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '515 515 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '516 516 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '517 517 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '518 518 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '519 519 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '520 520 RPP 0 50 -6 6 -.75 .75\n');
        
        fprintf(fileID_runfile_blockB, '601 601 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '602 602 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '603 603 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '604 604 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '605 605 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '606 606 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '607 607 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '608 608 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '609 609 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '610 610 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '611 611 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '612 612 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '613 613 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '614 614 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '615 615 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '616 616 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '617 617 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '618 618 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '619 619 RPP 0 50 -6 6 -.75 .75\n');
        fprintf(fileID_runfile_blockB, '620 620 RPP 0 50 -6 6 -.75 .75\n');
        
        fprintf(fileID_runfile_blockB, '420 px 533 $ Detector plane\n');
        
        fprintf(fileID_runfile_blockC, '\n');
        fprintf(fileID_runfile_blockC, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockC, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockC, 'C Block C\n');
        fprintf(fileID_runfile_blockC, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockC, 'C Block C: Tranformations\n');
        fprintf(fileID_runfile_blockC, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockC, ['TR501 428 ', num2str(6 + MLCleafPositions(1,1)), ' 14.25 1 0 0 0 1 0 0 0 1 $ Block A\n']);
        fprintf(fileID_runfile_blockC, ['TR502 428 ', num2str(6 + MLCleafPositions(1,2)), ' 12.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR503 428 ', num2str(6 + MLCleafPositions(1,3)), ' 11.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR504 428 ', num2str(6 + MLCleafPositions(1,4)), ' 9.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR505 428 ', num2str(6 + MLCleafPositions(1,5)), ' 8.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR506 428 ', num2str(6 + MLCleafPositions(1,6)), ' 6.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR507 428 ', num2str(6 + MLCleafPositions(1,7)), ' 5.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR508 428 ', num2str(6 + MLCleafPositions(1,8)), ' 3.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR509 428 ', num2str(6 + MLCleafPositions(1,9)), ' 2.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR510 428 ', num2str(6 + MLCleafPositions(1,10)), ' .75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR511 428 ', num2str(6 + MLCleafPositions(1,11)), ' -.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR512 428 ', num2str(6 + MLCleafPositions(1,12)), ' -2.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR513 428 ', num2str(6 + MLCleafPositions(1,13)), ' -3.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR514 428 ', num2str(6 + MLCleafPositions(1,14)), ' -5.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR515 428 ', num2str(6 + MLCleafPositions(1,15)), ' -6.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR516 428 ', num2str(6 + MLCleafPositions(1,16)), ' -8.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR517 428 ', num2str(6 + MLCleafPositions(1,17)), ' -9.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR518 428 ', num2str(6 + MLCleafPositions(1,18)), ' -11.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR519 428 ', num2str(6 + MLCleafPositions(1,19)), ' -12.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR520 428 ', num2str(6 + MLCleafPositions(1,20)), ' -14.25 1 0 0 0 1 0 0 0 1\n']);
        
        fprintf(fileID_runfile_blockC, ['TR601 428 ', num2str(MLCleafPositions(2,1) - 6), ' 14.25 1 0 0 0 1 0 0 0 1 $ Block B\n']);
        fprintf(fileID_runfile_blockC, ['TR602 428 ', num2str(MLCleafPositions(2,2) - 6), ' 12.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR603 428 ', num2str(MLCleafPositions(2,3) - 6), ' 11.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR604 428 ', num2str(MLCleafPositions(2,4) - 6), ' 9.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR605 428 ', num2str(MLCleafPositions(2,5) - 6), ' 8.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR606 428 ', num2str(MLCleafPositions(2,6) - 6), ' 6.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR607 428 ', num2str(MLCleafPositions(2,7) - 6), ' 5.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR608 428 ', num2str(MLCleafPositions(2,8) - 6), ' 3.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR609 428 ', num2str(MLCleafPositions(2,9) - 6), ' 2.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR610 428 ', num2str(MLCleafPositions(2,10) - 6), ' .75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR611 428 ', num2str(MLCleafPositions(2,11) - 6), ' -.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR612 428 ', num2str(MLCleafPositions(2,12) - 6), ' -2.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR613 428 ', num2str(MLCleafPositions(2,13) - 6), ' -3.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR614 428 ', num2str(MLCleafPositions(2,14) - 6), ' -5.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR615 428 ', num2str(MLCleafPositions(2,15) - 6), ' -6.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR616 428 ', num2str(MLCleafPositions(2,16) - 6), ' -8.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR617 428 ', num2str(MLCleafPositions(2,17) - 6), ' -9.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR618 428 ', num2str(MLCleafPositions(2,18) - 6), ' -11.25 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR619 428 ', num2str(MLCleafPositions(2,19) - 6), ' -12.75 1 0 0 0 1 0 0 0 1\n']);
        fprintf(fileID_runfile_blockC, ['TR620 428 ', num2str(MLCleafPositions(2,20) - 6), ' -14.25 1 0 0 0 1 0 0 0 1\n']);
        
        fclose(fileID_runfile_blockA);
        fclose(fileID_runfile_blockB);
        fclose(fileID_runfile_blockC);
        
    end

%% Function definition for FRM 2 MLC 3 MCNP runfile goes here...
    function runfileName = makeRunfileFRM2MLC3(MLCleafPositions, pathRunfiles, runfileName)
        
        % Prepare leaf dimenions s. th. leafs fit into available area
        horizontalDimMLCExitWindow = 265;
        leafDimensions_MLC3 = ...
            [ 18.5 14.5 14.5 14.5 14.5 7.5 11.5 11.5 11.5 18.2 11.5 11.5 11.5 7.5 14.5 14.5 14.5 14.5 18.5;...  % Back side
            14.9 14.5 14.5 14.5 14.5 17.5 11.5 11.5 11.5 8.2 11.5 11.5 11.5 17.5 14.5 14.5 14.5 14.5 14.8];     % Front side
        
        % add additional margin around leafs s.th. entrance window has width = 26.5cm
        leafDimensions_MLC3(1,:) = leafDimensions_MLC3(1,:) + .5;
        additionalMarginOuterLeafs = (horizontalDimMLCExitWindow - sum(leafDimensions_MLC3(1,:)))/2;
        leafDimensions_MLC3(1,1) = leafDimensions_MLC3(1,1) + additionalMarginOuterLeafs;
        leafDimensions_MLC3(1,end) = leafDimensions_MLC3(1,end) + additionalMarginOuterLeafs;
        
        % add additional margin around leafs s.th. entrance window has width = 26.5cm
        leafDimensions_MLC3(2,:) = leafDimensions_MLC3(2,:) + .35;
        additionalMarginOuterLeafs = (horizontalDimMLCExitWindow - sum(leafDimensions_MLC3(2,:)))/2;
        leafDimensions_MLC3(2,1) = leafDimensions_MLC3(2,1) + additionalMarginOuterLeafs;
        leafDimensions_MLC3(2,end) = leafDimensions_MLC3(2,end) + additionalMarginOuterLeafs;
        
        leafDimensions_MLC3 = leafDimensions_MLC3/10;   % rescale to cm for MCNP
        horizontalDimMLCExitWindow = horizontalDimMLCExitWindow/10;
        
        fileID_runfile_blockA = fopen(strcat(pathRunfiles,filesep,runfileName, '_blockA.txt'), 'w');
        fileID_runfile_blockB = fopen(strcat(pathRunfiles,filesep,runfileName, '_blockB.txt'), 'w');
        fileID_runfile_blockC = fopen(strcat(pathRunfiles,filesep,runfileName, '_blockC.txt'), 'w');
        
        % Write block A: Cells
        fprintf(fileID_runfile_blockA, 'MCNP Runfile for WSSA Generation for MLC3 at FRM2\n');
        fprintf(fileID_runfile_blockA, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockA, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockA, 'C Block A: Cells\n');
        fprintf(fileID_runfile_blockA, 'C ***************************************************************\n');
        
        fprintf(fileID_runfile_blockA, '1000 0 7000 (-1:-2:-3:-4:-5:-8) $ Propagation tunnel\n');
        fprintf(fileID_runfile_blockA, '1001 0  (-6 $ Start MLC $ Inside of MLC for fluence tally\n');
        fprintf(fileID_runfile_blockA, '        501 $ Bank A leaf 1\n');  % leaf bank A
        fprintf(fileID_runfile_blockA, '        502 $ Bank A leaf 2\n');
        fprintf(fileID_runfile_blockA, '        503 $ Bank A leaf 3\n');
        fprintf(fileID_runfile_blockA, '        504 $ Bank A leaf 4\n');
        fprintf(fileID_runfile_blockA, '        505 $ Bank A leaf 5\n');
        fprintf(fileID_runfile_blockA, '        506 $ Bank A leaf 6\n');
        fprintf(fileID_runfile_blockA, '        507 $ Bank A leaf 7\n');
        fprintf(fileID_runfile_blockA, '        508 $ Bank A leaf 8\n');
        fprintf(fileID_runfile_blockA, '        509 $ Bank A leaf 9\n');
        fprintf(fileID_runfile_blockA, '        510 $ Bank A leaf 10\n');
        fprintf(fileID_runfile_blockA, '        511 $ Bank A leaf 11\n');
        fprintf(fileID_runfile_blockA, '        512 $ Bank A leaf 12\n');
        fprintf(fileID_runfile_blockA, '        513 $ Bank A leaf 13\n');
        fprintf(fileID_runfile_blockA, '        514 $ Bank A leaf 14\n');
        fprintf(fileID_runfile_blockA, '        515 $ Bank A leaf 15\n');
        fprintf(fileID_runfile_blockA, '        516 $ Bank A leaf 16\n');
        fprintf(fileID_runfile_blockA, '        517 $ Bank A leaf 17\n');
        fprintf(fileID_runfile_blockA, '        518 $ Bank A leaf 18\n');
        fprintf(fileID_runfile_blockA, '        519 $ Bank A leaf 19\n');
        
        fprintf(fileID_runfile_blockA, '        601 $ Bank B leaf 1\n');  % leaf bank B
        fprintf(fileID_runfile_blockA, '        602 $ Bank B leaf 2\n');
        fprintf(fileID_runfile_blockA, '        603 $ Bank B leaf 3\n');
        fprintf(fileID_runfile_blockA, '        604 $ Bank B leaf 4\n');
        fprintf(fileID_runfile_blockA, '        605 $ Bank B leaf 5\n');
        fprintf(fileID_runfile_blockA, '        606 $ Bank B leaf 6\n');
        fprintf(fileID_runfile_blockA, '        607 $ Bank B leaf 7\n');
        fprintf(fileID_runfile_blockA, '        608 $ Bank B leaf 8\n');
        fprintf(fileID_runfile_blockA, '        609 $ Bank B leaf 9\n');
        fprintf(fileID_runfile_blockA, '        610 $ Bank B leaf 10\n');
        fprintf(fileID_runfile_blockA, '        611 $ Bank B leaf 11\n');
        fprintf(fileID_runfile_blockA, '        612 $ Bank B leaf 12\n');
        fprintf(fileID_runfile_blockA, '        613 $ Bank B leaf 13\n');
        fprintf(fileID_runfile_blockA, '        614 $ Bank B leaf 14\n');
        fprintf(fileID_runfile_blockA, '        615 $ Bank B leaf 15\n');
        fprintf(fileID_runfile_blockA, '        616 $ Bank B leaf 16\n');
        fprintf(fileID_runfile_blockA, '        617 $ Bank B leaf 17\n');
        fprintf(fileID_runfile_blockA, '        618 $ Bank B leaf 18\n');
        fprintf(fileID_runfile_blockA, '        619) $ Bank B leaf 19\n');
        
        
        fprintf(fileID_runfile_blockA, '9999 0 1 2 3 4 5 6 8 $ Graveyard\n');
        
        fprintf(fileID_runfile_blockA, '1002 0 -6 $ Exclude MLC\n');
        fprintf(fileID_runfile_blockA, '        (-501: $ Bank A leaf 1\n');  % leaf bank A
        fprintf(fileID_runfile_blockA, '        -502: $ Bank A leaf 2\n');
        fprintf(fileID_runfile_blockA, '        -503: $ Bank A leaf 3\n');
        fprintf(fileID_runfile_blockA, '        -504: $ Bank A leaf 4\n');
        fprintf(fileID_runfile_blockA, '        -505: $ Bank A leaf 5\n');
        fprintf(fileID_runfile_blockA, '        -506: $ Bank A leaf 6\n');
        fprintf(fileID_runfile_blockA, '        -507: $ Bank A leaf 7\n');
        fprintf(fileID_runfile_blockA, '        -508: $ Bank A leaf 8\n');
        fprintf(fileID_runfile_blockA, '        -509: $ Bank A leaf 9\n');
        fprintf(fileID_runfile_blockA, '        -510: $ Bank A leaf 10\n');
        fprintf(fileID_runfile_blockA, '        -511: $ Bank A leaf 11\n');
        fprintf(fileID_runfile_blockA, '        -512: $ Bank A leaf 12\n');
        fprintf(fileID_runfile_blockA, '        -513: $ Bank A leaf 13\n');
        fprintf(fileID_runfile_blockA, '        -514: $ Bank A leaf 14\n');
        fprintf(fileID_runfile_blockA, '        -515: $ Bank A leaf 15\n');
        fprintf(fileID_runfile_blockA, '        -516: $ Bank A leaf 16\n');
        fprintf(fileID_runfile_blockA, '        -517: $ Bank A leaf 17\n');
        fprintf(fileID_runfile_blockA, '        -518: $ Bank A leaf 18\n');
        fprintf(fileID_runfile_blockA, '        -519: $ Bank A leaf 19\n');
        
        fprintf(fileID_runfile_blockA, '        -601: $ Bank B leaf 1\n');  % leaf bank B
        fprintf(fileID_runfile_blockA, '        -602: $ Bank B leaf 2\n');
        fprintf(fileID_runfile_blockA, '        -603: $ Bank B leaf 3\n');
        fprintf(fileID_runfile_blockA, '        -604: $ Bank B leaf 4\n');
        fprintf(fileID_runfile_blockA, '        -605: $ Bank B leaf 5\n');
        fprintf(fileID_runfile_blockA, '        -606: $ Bank B leaf 6\n');
        fprintf(fileID_runfile_blockA, '        -607: $ Bank B leaf 7\n');
        fprintf(fileID_runfile_blockA, '        -608: $ Bank B leaf 8\n');
        fprintf(fileID_runfile_blockA, '        -609: $ Bank B leaf 9\n');
        fprintf(fileID_runfile_blockA, '        -610: $ Bank B leaf 10\n');
        fprintf(fileID_runfile_blockA, '        -611: $ Bank B leaf 11\n');
        fprintf(fileID_runfile_blockA, '        -612: $ Bank B leaf 12\n');
        fprintf(fileID_runfile_blockA, '        -613: $ Bank B leaf 13\n');
        fprintf(fileID_runfile_blockA, '        -614: $ Bank B leaf 14\n');
        fprintf(fileID_runfile_blockA, '        -615: $ Bank B leaf 15\n');
        fprintf(fileID_runfile_blockA, '        -616: $ Bank B leaf 16\n');
        fprintf(fileID_runfile_blockA, '        -617: $ Bank B leaf 17\n');
        fprintf(fileID_runfile_blockA, '        -618: $ Bank B leaf 18\n');
        fprintf(fileID_runfile_blockA, '        -619) $ Bank B leaf 19\n');
        
        fprintf(fileID_runfile_blockB, '\n');
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockB, 'C Block B: Surfaces\n');
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        
        fprintf(fileID_runfile_blockB, '1 RPP 0 62 -50 50 -50 50  $ Outer boundaries of simulation volume\n');
        fprintf(fileID_runfile_blockB, '4 RPP 288 423 -50 50 -50 50\n');
        fprintf(fileID_runfile_blockB, '8 RPP 523 700 -50 50 -50 50\n');
        
        fprintf(fileID_runfile_blockB, '2 ARB 62 -8 -9 $ Inner boundaries of transport tunnnel\n');
        fprintf(fileID_runfile_blockB, '        62 -8 9\n');
        fprintf(fileID_runfile_blockB, '        62 8 9\n');
        fprintf(fileID_runfile_blockB, '        62 8 -9\n');
        fprintf(fileID_runfile_blockB, '        274 -9 -12\n');
        fprintf(fileID_runfile_blockB, '        274 -9 12\n');
        fprintf(fileID_runfile_blockB, '        274 9 12\n');
        fprintf(fileID_runfile_blockB, '        274 9 -12\n');
        fprintf(fileID_runfile_blockB, '        1234 2367 1458 3478 1256 5678\n');
        fprintf(fileID_runfile_blockB, '3 RPP 274 288 -15 15 -17.5 17.5\n');
        fprintf(fileID_runfile_blockB, '5 RPP 423 463 -12 12 -15 15\n');
        fprintf(fileID_runfile_blockB, '6 RPP 463 523 -10 10 -13.25 13.25\n');
        
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockB, 'C MLC Leafs\n');
        fprintf(fileID_runfile_blockB, 'C ***************************************************************\n');
        % leaf 1
        leafCounter = 1;
        fprintf(fileID_runfile_blockB, [num2str(500 + leafCounter), ' ', num2str(500 + leafCounter), ' ARB 0 -5 ', ...
            num2str(-horizontalDimMLCExitWindow/2), '\n']);
        fprintf(fileID_runfile_blockB, ['        0 -5 ', num2str(-horizontalDimMLCExitWindow/2 +leafDimensions_MLC3(1,leafCounter)), '\n']);
        fprintf(fileID_runfile_blockB, ['        0 5 ', num2str(-horizontalDimMLCExitWindow/2 +leafDimensions_MLC3(1,leafCounter)), '\n']);
        fprintf(fileID_runfile_blockB, ['        0 5 ', num2str(-horizontalDimMLCExitWindow/2), '\n']);
        fprintf(fileID_runfile_blockB, ['        60 -5 ', num2str(-horizontalDimMLCExitWindow/2), '\n']);
        fprintf(fileID_runfile_blockB, ['        60 -5 ', num2str(-horizontalDimMLCExitWindow/2 +leafDimensions_MLC3(2,leafCounter)), '\n']);
        fprintf(fileID_runfile_blockB, ['        60 5 ', num2str(-horizontalDimMLCExitWindow/2 +leafDimensions_MLC3(2,leafCounter)), '\n']);
        fprintf(fileID_runfile_blockB, ['        60 5 ', num2str(-horizontalDimMLCExitWindow/2), '\n']);
        fprintf(fileID_runfile_blockB, '        1234 2367 1458 3478 1256 5678\n');
        
        % leafs 2:end
        for leafCounter = 2:size(leafDimensions_MLC3,2)
            fprintf(fileID_runfile_blockB, [num2str(500 + leafCounter), ' ', num2str(500 + leafCounter), ' ARB 0 -5 ', ...
                num2str(-horizontalDimMLCExitWindow/2 + sum(leafDimensions_MLC3(1,1:leafCounter-1))), '\n']);
            fprintf(fileID_runfile_blockB, ['        0 -5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(1,1:leafCounter))), '\n']);
            fprintf(fileID_runfile_blockB, ['        0 5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(1,1:leafCounter))), '\n']);
            fprintf(fileID_runfile_blockB, ['        0 5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(1,1:leafCounter-1))), '\n']);
            fprintf(fileID_runfile_blockB, ['        60 -5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(2,1:leafCounter-1))), '\n']);
            fprintf(fileID_runfile_blockB, ['        60 -5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(2,1:leafCounter))), '\n']);
            fprintf(fileID_runfile_blockB, ['        60 5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(2,1:leafCounter))), '\n']);
            fprintf(fileID_runfile_blockB, ['        60 5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(2,1:leafCounter-1))), '\n']);
            fprintf(fileID_runfile_blockB, '        1234 2367 1458 3478 1256 5678\n');
        end
        
        % leaf 1
        leafCounter = 1;
        fprintf(fileID_runfile_blockB, [num2str(600 + leafCounter), ' ', num2str(600 + leafCounter), ' ARB 0 -5 ', ...
            num2str(-horizontalDimMLCExitWindow/2), '\n']);
        fprintf(fileID_runfile_blockB, ['        0 -5 ', num2str(-horizontalDimMLCExitWindow/2 +leafDimensions_MLC3(1,leafCounter)), '\n']);
        fprintf(fileID_runfile_blockB, ['        0 5 ', num2str(-horizontalDimMLCExitWindow/2 +leafDimensions_MLC3(1,leafCounter)), '\n']);
        fprintf(fileID_runfile_blockB, ['        0 5 ', num2str(-horizontalDimMLCExitWindow/2), '\n']);
        fprintf(fileID_runfile_blockB, ['        60 -5 ', num2str(-horizontalDimMLCExitWindow/2), '\n']);
        fprintf(fileID_runfile_blockB, ['        60 -5 ', num2str(-horizontalDimMLCExitWindow/2 +leafDimensions_MLC3(2,leafCounter)), '\n']);
        fprintf(fileID_runfile_blockB, ['        60 5 ', num2str(-horizontalDimMLCExitWindow/2 +leafDimensions_MLC3(2,leafCounter)), '\n']);
        fprintf(fileID_runfile_blockB, ['        60 5 ', num2str(-horizontalDimMLCExitWindow/2), '\n']);
        fprintf(fileID_runfile_blockB, '        1234 2367 1458 3478 1256 5678\n');
        
        % leafs 2:end
        for leafCounter = 2:size(leafDimensions_MLC3,2)
            fprintf(fileID_runfile_blockB, [num2str(600 + leafCounter), ' ', num2str(600 + leafCounter), ' ARB 0 -5 ', ...
                num2str(-horizontalDimMLCExitWindow/2 + sum(leafDimensions_MLC3(1,1:leafCounter-1))), '\n']);
            fprintf(fileID_runfile_blockB, ['        0 -5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(1,1:leafCounter))), '\n']);
            fprintf(fileID_runfile_blockB, ['        0 5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(1,1:leafCounter))), '\n']);
            fprintf(fileID_runfile_blockB, ['        0 5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(1,1:leafCounter-1))), '\n']);
            fprintf(fileID_runfile_blockB, ['        60 -5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(2,1:leafCounter-1))), '\n']);
            fprintf(fileID_runfile_blockB, ['        60 -5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(2,1:leafCounter))), '\n']);
            fprintf(fileID_runfile_blockB, ['        60 5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(2,1:leafCounter))), '\n']);
            fprintf(fileID_runfile_blockB, ['        60 5 ', num2str(-horizontalDimMLCExitWindow/2 + ...
                sum(leafDimensions_MLC3(2,1:leafCounter-1))), '\n']);
            fprintf(fileID_runfile_blockB, '        1234 2367 1458 3478 1256 5678\n');
        end
        
        
        fprintf(fileID_runfile_blockB, '420 px 533 $ Detector plane\n');
        
        fprintf(fileID_runfile_blockC, '\n');
        fprintf(fileID_runfile_blockC, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockC, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockC, 'C Block C\n');
        fprintf(fileID_runfile_blockC, 'C ***************************************************************\n');
        fprintf(fileID_runfile_blockC, 'C Block C: Tranformations\n');
        fprintf(fileID_runfile_blockC, 'C ***************************************************************\n');
        for leafCounter = 1:size(leafDimensions_MLC3,2)
            if leafCounter == 1
                fprintf(fileID_runfile_blockC, ['TR', num2str(500 + leafCounter), ' 463 ', num2str(5 + MLCleafPositions(1,leafCounter)), ' ', ...
                    '0' ,...
                    ' 1 0 0 0 1 0 0 0 1 $ Block A\n']);
            else
                fprintf(fileID_runfile_blockC, ['TR', num2str(500 + leafCounter), ' 463 ', num2str(5 + MLCleafPositions(1,leafCounter)), ' ', ...
                    '0' ,...
                    ' 1 0 0 0 1 0 0 0 1 $ Block A\n']);
            end
        end
        
        for leafCounter = 1:size(leafDimensions_MLC3,2)
            if leafCounter == 1
                fprintf(fileID_runfile_blockC, ['TR', num2str(600 + leafCounter), ' 463 ', num2str(MLCleafPositions(2,1) - 5), ' ', ...
                    '0' ,...
                    ' 1 0 0 0 1 0 0 0 1 $ Block B\n']);
            else
                fprintf(fileID_runfile_blockC, ['TR', num2str(600 + leafCounter), ' 463 ', num2str(MLCleafPositions(2,leafCounter) - 5), ' ', ...
                    '0' ,...
                    ' 1 0 0 0 1 0 0 0 1 $ Block B\n']);
            end
        end
        
        
        fclose(fileID_runfile_blockA);
        fclose(fileID_runfile_blockB);
        fclose(fileID_runfile_blockC);
        
    end

%% Function definition for FRM 2 PTW water phantom (52x63.5x63.5cm^3)
    function [waterPhantomRunfileName, WSSAfileName] = makeWaterPhantomFRM2PTW(runfileName, SSD4kernelCalc, npsKernelGen, particleType, transportMedium, sourceOpeningAngle, WSSAoption) % SSD defined from wall to water phantom surface
        phantomDimensions = [63.5 52 63.5];
        wallPosition = 493;
        diameterCentralDetector = 1.5;
        waterPhantomRunfileName = [runfileName, '_waterPhantom'];
        fileID_waterPhantom_blockA = fopen(strcat(pwd,filesep, [waterPhantomRunfileName, '_blockA.txt']), 'w');
        fileID_waterPhantom_blockB = fopen(strcat(pwd,filesep, [waterPhantomRunfileName, '_blockB.txt']), 'w');
        fileID_waterPhantom_blockC = fopen(strcat(pwd,filesep, [waterPhantomRunfileName, '_blockC.txt']), 'w');
        
        fprintf(fileID_waterPhantom_blockA, 'C ***************************************************************\n');
        fprintf(fileID_waterPhantom_blockA, 'C Block A: Cells Waterphantom\n');
        fprintf(fileID_waterPhantom_blockA, 'C ***************************************************************\n');
        
        fprintf(fileID_waterPhantom_blockB, 'C ***************************************************************\n');
        fprintf(fileID_waterPhantom_blockB, 'C Block B: Surfaces Waterphantom\n');
        fprintf(fileID_waterPhantom_blockB, 'C ***************************************************************\n');
        
        switch particleType
            case 'neutronField'
                % Get list of available tabulated neutron spectra
                spectralInformation.pathLocation = fullfile(matRad_getMATRADdirectory_hardCoded,'MCNP', 'SpectralInformation', filesep);
                spectralInformation.neutronSpectrum = dir([spectralInformation.pathLocation, 'spectrum_neutrons_*']);
                
                % Check if there is more than one neutron spectrum available in
                % ../MATRAD/MCNP/SpectralInformation
                if size(spectralInformation.neutronSpectrum,1)~=1
                    [spectralInformation.neutronIndex,~] = listdlg('PromptString','Please select neutron spectrum:',...
                        'SelectionMode','single',...
                        'ListString',{spectralInformation.neutronSpectrum.name});
                else
                    spectralInformation.neutronIndex = 1;
                end
                
                % Read spectral information from selected file with first column as
                % energy and second column as spectral information
                fid_neutronSpectrum = fopen([spectralInformation.neutronSpectrum(spectralInformation.neutronIndex).folder, filesep, spectralInformation.neutronSpectrum(spectralInformation.neutronIndex).name], 'r');
                spectralInformation.spectrumValuesNeutrons = fscanf(fid_neutronSpectrum, '%f', [2,inf]);
                spectralInformation.spectrumValuesNeutrons = spectralInformation.spectrumValuesNeutrons';
                fclose(fid_neutronSpectrum);
                
                % Define source card, note: VEC=reference vector for the direction
                % sampling, DIR=cosine of angle between VEC and partice direction,
                % in case DIR=-1 a monodirectional source in counter direction of
                % VEC.
                % ERG=d3 used to define spectrum according to information read from
                % tabulated data in ..\MATRAD\MCNP\SpectralInformation
                
                source.sourceCard_0 = 'SDEF\n        POS=0 0 0\n        X=0 Y=d1 Z=d2\n        VEC=1 0 0\n        DIR=d4 PAR=1 ERG=d3 $wgt=3.5088e+03  $ changed from HB (wgt=30.65e6) since primer states use 1/fsa2 (fsa2=0.000285)\n';
                source.sourceCard_1_i = 'SI1 -7.5 7.5\n';  % Initial position and source extension
                source.sourceCard_1_p = 'SP1 0 1\n';
                source.sourceCard_2_i = 'SI2 -7.5 7.5\n';  % ...
                source.sourceCard_2_p = 'SP2 0 1\n';
                source.energyCard_3_i_0 = 'SI3 H\n';        % Energy bins
                source.energyCard_3_i = '        %8d\n';
                source.energyCard_3_p_0 = 'SP3 D\n';        % Spectral information
                source.energyCard_3_p = '        %8d\n';
                source.energyCard_4_i = ['SI4   -1 ', num2str(cosd(sourceOpeningAngle)), ' 1\n'];
                source.energyCard_4_p = ['SP4    0 ', num2str(1-4*pi*(sin((sourceOpeningAngle/360)*2*pi/2)^2)/(4*pi)),' ', num2str(4*pi*(sin((sourceOpeningAngle/360)*2*pi/2)^2)/(4*pi)), '\n'];
                source.energyCard_4_b = 'SB4 0 0 1\n';
                
                % Write source block
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C C: Source\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                
                % Write initial source position and extension
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_0);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_1_i);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_1_p);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_2_i);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_2_p);
                
                % Write spectral distribution
                fprintf(fileID_waterPhantom_blockC, source.energyCard_3_i_0);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_3_i, ...
                    spectralInformation.spectrumValuesNeutrons(:,1));
                fprintf(fileID_waterPhantom_blockC, source.energyCard_3_p_0);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_3_p, ...
                    spectralInformation.spectrumValuesNeutrons(:,2));
                fprintf(fileID_waterPhantom_blockC, source.energyCard_4_i);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_4_p);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_4_b);
                
            case 'photonField'
                % Get list of available tabulated photon spectra
                spectralInformation.pathLocation = fullfile(matRad_getMATRADdirectory_hardCoded,'MCNP', 'SpectralInformation', filesep);
                spectralInformation.photonSpectrum = dir([spectralInformation.pathLocation, 'spectrum_photons_*']);
                
                % Check if there is more than one photon spectra available in
                % ../MATRAD/MCNP/SpectralInformation
                if size(spectralInformation.photonSpectrum,1)~=1
                    [spectralInformation.photonIndex,~] = listdlg('PromptString','Please select photon spectrum:',...
                        'SelectionMode','single',...
                        'ListString',{spectralInformation.photonSpectrum.name});
                else
                    spectralInformation.photonIndex = 1;
                end
                
                % Read spectral information from selected file with first column as
                % energy and second column as spectral information
                fid_photonSpectrum = fopen([spectralInformation.photonSpectrum(spectralInformation.photonIndex).folder, filesep, spectralInformation.photonSpectrum(spectralInformation.photonIndex).name], 'r');
                spectralInformation.spectrumValuesPhotons = fscanf(fid_photonSpectrum, '%f', [2,inf]);
                spectralInformation.spectrumValuesPhotons = spectralInformation.spectrumValuesPhotons';
                fclose(fid_photonSpectrum);
                
                % Define source card, note: VEC=reference vector for the direction
                % sampling, DIR=cosine of angle between VEC and partice direction,
                % in case DIR=-1 a monodirectional source in counter direction of
                % VEC.
                % ERG=d3 used to define spectrum according to information read from
                % tabulated data in ..\MATRAD\MCNP\SpectralInformation
                
                source.sourceCard_0 = 'SDEF\n        POS=0 0 0\n        X=0 Y=d1 Z=d2\n        VEC=1 0 0\n        DIR=d4 PAR=2 ERG=d3        $ wgt=3.5088e+03  $ changed from HB (wgt=30.65e6) since primer states use 1/fsa2 (fsa2=0.000285)\n';
                source.sourceCard_1_i = 'SI1 -7.5 7.5\n';  % Initial position and source extension
                source.sourceCard_1_p = 'SP1 0 1\n';
                source.sourceCard_2_i = 'SI2 -7.5 7.5\n';  % ...
                source.sourceCard_2_p = 'SP2 0 1\n';
                source.energyCard_3_i_0 = 'SI3 H\n';        % Energy bins
                source.energyCard_3_i = '        %8d\n';
                source.energyCard_3_p_0 = 'SP3 D\n';        % Spectral information
                source.energyCard_3_p = '        %8d\n';
                source.energyCard_4_i = ['SI4   -1 ', num2str(cosd(sourceOpeningAngle)), ' 1\n'];
                source.energyCard_4_p = ['SP4    0 ', num2str(1-4*pi*(sin((sourceOpeningAngle/360)*2*pi/2)^2)/(4*pi)),' ', num2str(4*pi*(sin((sourceOpeningAngle/360)*2*pi/2)^2)/(4*pi)), '\n'];
                source.energyCard_4_b = 'SB4 0 0 1\n';
                
                % Write source block
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C C: Source\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                
                % Write initial source position and extension
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_0);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_1_i);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_1_p);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_2_i);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_2_p);
                
                % Write spectral distribution
                fprintf(fileID_waterPhantom_blockC, source.energyCard_3_i_0);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_3_i, ...
                    spectralInformation.spectrumValuesPhotons(:,1));
                fprintf(fileID_waterPhantom_blockC, source.energyCard_3_p_0);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_3_p, ...
                    spectralInformation.spectrumValuesPhotons(:,2));
                fprintf(fileID_waterPhantom_blockC, source.energyCard_4_i);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_4_p);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_4_b);
                
            case 'mixedField'
                % Get list of available tabulated neutron spectra
                spectralInformation.pathLocation = fullfile(matRad_getMATRADdirectory_hardCoded,'MCNP', 'SpectralInformation', filesep);
                spectralInformation.neutronSpectrum = dir([spectralInformation.pathLocation, 'spectrum_neutrons_*']);
                spectralInformation.photonSpectrum = dir([spectralInformation.pathLocation, 'spectrum_photons_*']);
                
                % Check if there is more than one neutron spectrum available in
                % ../MATRAD/MCNP/SpectralInformation
                if size(spectralInformation.neutronSpectrum,1)~=1
                    [spectralInformation.neutronIndex,~] = listdlg('PromptString','Please select neutron spectrum:',...
                        'SelectionMode','single',...
                        'ListString',{spectralInformation.neutronSpectrum.name});
                else
                    spectralInformation.neutronIndex = 1;
                end
                % Same for photons
                if size(spectralInformation.photonSpectrum,1)~=1
                    [spectralInformation.photonIndex,~] = listdlg('PromptString','Please select photon spectrum:',...
                        'SelectionMode','single',...
                        'ListString',{spectralInformation.photonSpectrum.name});
                else
                    spectralInformation.photonIndex = 1;
                end
                
                % Read spectral information from selected file with first column as
                % energy and second column as spectral information
                % Neutrons
                fid_neutronSpectrum = fopen([spectralInformation.neutronSpectrum(spectralInformation.neutronIndex).folder, filesep, spectralInformation.neutronSpectrum(spectralInformation.neutronIndex).name], 'r');
                spectralInformation.spectrumValuesNeutrons = fscanf(fid_neutronSpectrum, '%f', [2,inf]);
                spectralInformation.spectrumValuesNeutrons = spectralInformation.spectrumValuesNeutrons';
                fclose(fid_neutronSpectrum);
                % Photons
                fid_photonSpectrum = fopen([spectralInformation.photonSpectrum(spectralInformation.photonIndex).folder, filesep, spectralInformation.photonSpectrum(spectralInformation.photonIndex).name], 'r');
                spectralInformation.spectrumValuesPhotons = fscanf(fid_photonSpectrum, '%f', [2,inf]);
                spectralInformation.spectrumValuesPhotons = spectralInformation.spectrumValuesPhotons';
                fclose(fid_photonSpectrum);
                
                % Define source card, note: VEC=reference vector for the direction
                % sampling, DIR=cosine of angle between VEC and partice direction,
                % in case DIR=-1 a monodirectional source in counter direction of
                % VEC.
                % ERG=d3 used to define spectrum according to information read from
                % tabulated data in ..\MATRAD\MCNP\SpectralInformation
                
                source.sourceCard_0 = 'SDEF\n        POS=0 0 0\n        X=0 Y=d1 Z=d2\n        VEC=1 0 0\n        PAR=d3 ERG=fpar=d4 DIR=d7\n';
                source.sourceCard_1_i = 'SI1 -7.5 7.5\n';  % Initial position and source extension
                source.sourceCard_1_p = 'SP1 0 1\n';
                source.sourceCard_2_i = 'SI2 -7.5 7.5\n';  % ...
                source.sourceCard_2_p = 'SP2 0 1\n';
                source.particleDistribution_i = 'SI3 L 1 2\n';
                source.particleDistribution_p = 'SP3 3.2 2.9\n';    % 3.2e8 n/cm^2/s and 2.9 gammas/cm^2/s
                source.particleEnergyDistributions = 'DS4 S 5 6\n';
                source.neutronEnergyCard_5_i_0 = 'SI5 H\n';        % Energy bins neutrons
                source.neutronEnergyCard_5_i = '        %8d\n';
                source.neutronEnergyCard_5_p_0 = 'SP5 D\n';        % Spectral information neutrons
                source.neutronEnergyCard_5_p = '        %8d\n';
                source.photonEnergyCard_6_i_0 = 'SI6 H\n';        % Energy bins photons
                source.photonEnergyCard_6_i = '        %8d\n';
                source.photonEnergyCard_6_p_0 = 'SP6 D\n';        % Spectral information photons
                source.photonEnergyCard_6_p = '        %8d\n';
                source.energyCard_7_i = ['SI7   -1 ', num2str(cosd(sourceOpeningAngle)), ' 1\n'];
                source.energyCard_7_p = ['SP7    0 ', num2str(1-4*pi*(sin((sourceOpeningAngle/360)*2*pi/2)^2)/(4*pi)),' ', num2str(4*pi*(sin((sourceOpeningAngle/360)*2*pi/2)^2)/(4*pi)), '\n'];
                source.energyCard_7_b = 'SB7 0 0 1\n';
                
                % Write source block
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C C: Source\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                
                % Write initial source position and extension
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_0);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_1_i);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_1_p);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_2_i);
                fprintf(fileID_waterPhantom_blockC, source.sourceCard_2_p);
                
                % Write spectral distribution
                fprintf(fileID_waterPhantom_blockC, source.particleDistribution_i);
                fprintf(fileID_waterPhantom_blockC, source.particleDistribution_p);
                fprintf(fileID_waterPhantom_blockC, source.particleEnergyDistributions);
                % Neutrons
                fprintf(fileID_waterPhantom_blockC, source.neutronEnergyCard_5_i_0);
                fprintf(fileID_waterPhantom_blockC, source.neutronEnergyCard_5_i, ...
                    spectralInformation.spectrumValuesNeutrons(:,1));
                fprintf(fileID_waterPhantom_blockC, source.neutronEnergyCard_5_p_0);
                fprintf(fileID_waterPhantom_blockC, source.neutronEnergyCard_5_p, ...
                    spectralInformation.spectrumValuesNeutrons(:,2));
                % Photons
                fprintf(fileID_waterPhantom_blockC, source.photonEnergyCard_6_i_0);
                fprintf(fileID_waterPhantom_blockC, source.photonEnergyCard_6_i, ...
                    spectralInformation.spectrumValuesPhotons(:,1));
                fprintf(fileID_waterPhantom_blockC, source.photonEnergyCard_6_p_0);
                fprintf(fileID_waterPhantom_blockC, source.photonEnergyCard_6_p, ...
                    spectralInformation.spectrumValuesPhotons(:,2));
                fprintf(fileID_waterPhantom_blockC, source.energyCard_7_i);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_7_p);
                fprintf(fileID_waterPhantom_blockC, source.energyCard_7_b);
                
        end       
        
        
        fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
        switch transportMedium
            case 'water'
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C Block C: Waterphantom\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                
                fprintf(fileID_waterPhantom_blockA, '7000 1 -1 -7000 7999 $ Water phantom w/o central PDD detector\n');
                fprintf(fileID_waterPhantom_blockB, ['7000 RPP ', num2str(wallPosition + SSD4kernelCalc), ' ', ...
                    num2str(wallPosition + SSD4kernelCalc + phantomDimensions(1)), ' ', ...
                    num2str(-phantomDimensions(2)/2), ' ', num2str(phantomDimensions(2)/2), ' ', ...
                    num2str(-phantomDimensions(3)/2), ' ', num2str(phantomDimensions(3)/2), ' $ Water phantom outer contour\n']);
                fprintf(fileID_waterPhantom_blockB, ['7999 CX ', num2str(diameterCentralDetector/2), '\n']);
                
                fprintf(fileID_waterPhantom_blockC, 'MODE N P E H D T S A #\n');
                fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'PHYS:N 101 0 0 3J 1 -1 3J 0 0\n');
                fprintf(fileID_waterPhantom_blockC, 'PHYS:H 101 101 -1 J 0 J 1 3J 0 0 0 0.917\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:n J 0\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:p J 1.0E-03\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:e J 1.0E-03\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:h J 1.0E-06\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:d J 1.0E-06\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:t J 1.0E-06\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:s J 1.0E-06\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:a J 1.0E-06\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:# J 1.0E-06\n');
                fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
            case 'none'
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C Block C: Waterphantom\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                
                fprintf(fileID_waterPhantom_blockA, '7000 1 -1 -7000 $ Water phantom\n');
                fprintf(fileID_waterPhantom_blockB, ['7000 RPP ', num2str(wallPosition + SSD4kernelCalc), ' ', ...
                    num2str(wallPosition + SSD4kernelCalc + phantomDimensions(1)), ' ', ...
                    num2str(-phantomDimensions(2)/2), ' ', num2str(phantomDimensions(2)/2), ' ', ...
                    num2str(-phantomDimensions(3)/2), ' ', num2str(phantomDimensions(3)/2), ' $ Water phantom outer contour\n']);
                
                fprintf(fileID_waterPhantom_blockC, 'MODE N P\n');
                fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'PHYS:N 101 0 0 3J 1 -1 3J 0 0\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:n J 0\n');
                fprintf(fileID_waterPhantom_blockC, 'CUT:p J 1.0E-03\n');
                fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
        end
        fprintf(fileID_waterPhantom_blockC, ['NPS ', num2str(npsKernelGen), '\n']);
        fprintf(fileID_waterPhantom_blockC, ['PRDMP J ', num2str(npsKernelGen/2),' 1 1\n']);
        
        fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
        fprintf(fileID_waterPhantom_blockC, 'C Transport medium\n');
        fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
        switch transportMedium
            case 'water'
                fprintf(fileID_waterPhantom_blockC, 'm1   1001.66c   -0.112102272\n');
                fprintf(fileID_waterPhantom_blockC, '        1002.80c   -2.57667e-05\n');
                fprintf(fileID_waterPhantom_blockC, '        8016.80c   -0.887512659\n');
                fprintf(fileID_waterPhantom_blockC, '        8017.80c   -0.000359302\n');
                fprintf(fileID_waterPhantom_blockC, 'mt1  lwtr.20t\n');
                fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
            case 'none'
                fprintf(fileID_waterPhantom_blockC, 'm1   7014.70c -0.755000000\n');
                fprintf(fileID_waterPhantom_blockC, '        8016.70c -0.232000000\n');
                fprintf(fileID_waterPhantom_blockC, '        18040.70c -0.013000000\n');
                fprintf(fileID_waterPhantom_blockC, 'C ************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'void\n');
        end
        
        if WSSAoption
            fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
            fprintf(fileID_waterPhantom_blockC, 'C Block C: Write WSSA\n');
            fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                    switch transportMedium
                        case 'water'
                                        fprintf(fileID_waterPhantom_blockC, 'SSW 9.1 (-9) PTY=N P\n');
                        case 'none'
                                        fprintf(fileID_waterPhantom_blockC, 'SSW 9.1 (-9) PTY=N P\n');
                    end
            WSSAfileName = 'WSSAGenerated';
        else
            WSSAfileName = 'WSSAnotGenerated';
        end
        
        
        fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
        fprintf(fileID_waterPhantom_blockC, 'C Block C: Tallies Waterphantom\n');
        fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
        switch transportMedium
            case 'water'
                numberTalliesPDD = [35 150 60];
                tallyResolution = [.1 .2 .5];
                tallyCounter = 1;
                fprintf(fileID_waterPhantom_blockB, [num2str(7000+tallyCounter),' PX ',num2str(wallPosition + SSD4kernelCalc), '\n']);
                
                while tallyCounter <= numberTalliesPDD(1)
                    % Write cell cards
                    fprintf(fileID_waterPhantom_blockA, [num2str(7000+tallyCounter),' 1 -1 -7999 ',num2str(7000+tallyCounter), ' -',num2str(7000+tallyCounter+1), '\n']);
                    % Write surface cards
                    fprintf(fileID_waterPhantom_blockB, [num2str(7000+tallyCounter+1),' PX ',num2str(wallPosition + SSD4kernelCalc + (tallyCounter)*tallyResolution(1)), '\n']);
                    % Write tallies
                    fprintf(fileID_waterPhantom_blockC, ['+F',num2str(tallyCounter), '16 ', num2str(7000+tallyCounter), '\n']);
                    fprintf(fileID_waterPhantom_blockC, ['SF',num2str(tallyCounter), '16 7999\n']);
                    tallyCounter = tallyCounter +1;
                end
                
                while tallyCounter <= numberTalliesPDD(1)+numberTalliesPDD(2)
                    % Write cell cards
                    fprintf(fileID_waterPhantom_blockA, [num2str(7000+tallyCounter),' 1 -1 -7999 ',num2str(7000+tallyCounter), ' -',num2str(7000+tallyCounter+1), '\n']);
                    % Write surface cards
                    fprintf(fileID_waterPhantom_blockB, [num2str(7000+tallyCounter+1),' PX ',num2str(wallPosition + SSD4kernelCalc + numberTalliesPDD(1)*tallyResolution(1) + (tallyCounter-numberTalliesPDD(1))*tallyResolution(2)), '\n']);
                    % Write tallies
                    fprintf(fileID_waterPhantom_blockC, ['+F',num2str(tallyCounter), '16 ', num2str(7000+tallyCounter), '\n']);
                    fprintf(fileID_waterPhantom_blockC, ['SF',num2str(tallyCounter), '16 7999\n']);
                    tallyCounter = tallyCounter +1;
                end
                
                while tallyCounter <= sum(numberTalliesPDD)
                    % Write cell cards
                    fprintf(fileID_waterPhantom_blockA, [num2str(7000+tallyCounter),' 1 -1 -7999 ',num2str(7000+tallyCounter), ' -',num2str(7000+tallyCounter+1), '\n']);
                    % Write surface cards
                    fprintf(fileID_waterPhantom_blockB, [num2str(7000+tallyCounter+1),' PX ',num2str(wallPosition + SSD4kernelCalc + numberTalliesPDD(1)*tallyResolution(1) + numberTalliesPDD(2)*tallyResolution(2) + (tallyCounter-numberTalliesPDD(1)-numberTalliesPDD(2))*tallyResolution(3)), '\n']);
                    % Write tallies
                    fprintf(fileID_waterPhantom_blockC, ['+F',num2str(tallyCounter), '16 ', num2str(7000+tallyCounter), '\n']);
                    fprintf(fileID_waterPhantom_blockC, ['SF',num2str(tallyCounter), '16 7999\n']);
                    tallyCounter = tallyCounter +1;
                end
                
                
                % Fluence tallies inside MLC
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C Fluence Tallies:\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                
                fprintf(fileID_waterPhantom_blockC, 'F10014:N 1001\n'); %n neutrons
                fprintf(fileID_waterPhantom_blockC, 'E10014\n');
                energyResolution = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1:1e-1:.9 1:.25:11];
                for counterEnergyRes = 1:length(energyResolution)
                    fprintf(fileID_waterPhantom_blockC, ['        ', num2str(energyResolution(counterEnergyRes)), '\n']);
                end
                fprintf(fileID_waterPhantom_blockC, 'F10024:P 1001\n'); % photons
                fprintf(fileID_waterPhantom_blockC, 'E10024\n');
                energyResolution = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1:1e-1:.9 1:.25:11];
                for counterEnergyRes = 1:length(energyResolution)
                    fprintf(fileID_waterPhantom_blockC, ['        ', num2str(energyResolution(counterEnergyRes)), '\n']);
                end
                
                % TMESH tally inside water phantom
                meshTally.typeCard = 'TMESH\n';
                meshTally.geometry = 'RMESH3 %s\n';
                meshTally.corA = 'CORA3 %.4f %dI %.4f\n';
                meshTally.corB = 'CORB3 %.4f %dI %.4f\n';
                meshTally.corC = 'CORC3 %.4f %dI %.4f\n';
                tallyKeyword= 'TOTAL';
                meshTally.resolution = .25;
                
                
                %Write to text file
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C TMESH Tally:\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, meshTally.typeCard);
                fprintf(fileID_waterPhantom_blockC, meshTally.geometry, tallyKeyword);
                fprintf(fileID_waterPhantom_blockC, meshTally.corA, wallPosition + SSD4kernelCalc, (phantomDimensions(1)/meshTally.resolution-1), wallPosition + SSD4kernelCalc + phantomDimensions(1));
                fprintf(fileID_waterPhantom_blockC, meshTally.corB, -phantomDimensions(2)/2, (phantomDimensions(2)/meshTally.resolution-1), phantomDimensions(2)/2);
                fprintf(fileID_waterPhantom_blockC, meshTally.corC, -phantomDimensions(3)/2, (phantomDimensions(3)/meshTally.resolution-1), phantomDimensions(3)/2);
                fprintf(fileID_waterPhantom_blockC, 'ENDMD\n');
            case 'none'
                disp('Water phantom not set up in simulation.')
                numberTalliesPDD = [35 150 60];
                % Fluence tallies inside MLC
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C Fluence Tallies:\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                
                fprintf(fileID_waterPhantom_blockC, 'F10014:N 1001\n'); %n neutrons
                fprintf(fileID_waterPhantom_blockC, 'E10014\n');
                energyResolution = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1:1e-1:.9 1:.25:11];
                for counterEnergyRes = 1:length(energyResolution)
                    fprintf(fileID_waterPhantom_blockC, ['        ', num2str(energyResolution(counterEnergyRes)), '\n']);
                end
                fprintf(fileID_waterPhantom_blockC, 'F10024:P 1001\n'); % photons
                fprintf(fileID_waterPhantom_blockC, 'E10024\n');
                energyResolution = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1:1e-1:.9 1:.25:11];
                for counterEnergyRes = 1:length(energyResolution)
                    fprintf(fileID_waterPhantom_blockC, ['        ', num2str(energyResolution(counterEnergyRes)), '\n']);
                end
                
                % TMESH tally inside water phantom
                meshTally.typeCard = 'TMESH\n';
                meshTally.geometry = 'RMESH3 %s\n';
                meshTally.corA = 'CORA3 %.4f %dI %.4f\n';
                meshTally.corB = 'CORB3 %.4f %dI %.4f\n';
                meshTally.corC = 'CORC3 %.4f %dI %.4f\n';
                tallyKeyword= 'TOTAL';
                meshTally.resolution = .25;
                
                
                %Write to text file
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, 'C TMESH Tally:\n');
                fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
                fprintf(fileID_waterPhantom_blockC, meshTally.typeCard);
                fprintf(fileID_waterPhantom_blockC, meshTally.geometry, tallyKeyword);
                fprintf(fileID_waterPhantom_blockC, meshTally.corA, wallPosition + SSD4kernelCalc, (phantomDimensions(1)/meshTally.resolution-1), wallPosition + SSD4kernelCalc + phantomDimensions(1));
                fprintf(fileID_waterPhantom_blockC, meshTally.corB, -phantomDimensions(2)/2, (phantomDimensions(2)/meshTally.resolution-1), phantomDimensions(2)/2);
                fprintf(fileID_waterPhantom_blockC, meshTally.corC, -phantomDimensions(3)/2, (phantomDimensions(3)/meshTally.resolution-1), phantomDimensions(3)/2);
                fprintf(fileID_waterPhantom_blockC, 'ENDMD\n');
                
        end
        
        
        fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
        fprintf(fileID_waterPhantom_blockC, 'C Block C: Physics\n');
        fprintf(fileID_waterPhantom_blockC, 'C ***************************************************************\n');
        switch transportMedium
            case 'water'
                fprintf(fileID_waterPhantom_blockC, 'IMP:N,P,E,H,D,T,S,A,# 1 1 1 0 0\n');
                for importanceCounter = 1:sum(numberTalliesPDD)+1
                    fprintf(fileID_waterPhantom_blockC, '        1\n');
                end

            case 'none'
                fprintf(fileID_waterPhantom_blockC, 'IMP:N,P 1 1 1 0 0 1\n');               
        end
        
        fprintf(fileID_waterPhantom_blockC, 'RAND GEN=2 SEED=43 STRIDE=10000000001\n');
        
        
        fclose(fileID_waterPhantom_blockA);
        fclose(fileID_waterPhantom_blockB);
        fclose(fileID_waterPhantom_blockC);
        
    end

%% Function to concatenate runfiles
    function runfileName = concatenateRunfile(runfileName, waterPhantomRunfileName)
        blockA_0 = [runfileName, '_blockA.txt'];
        blockB_0 = [runfileName, '_blockB.txt'];
        blockC_0 = [runfileName, '_blockC.txt'];
        
        blockA_waterPhantom =  [waterPhantomRunfileName, '_blockA.txt'];
        blockB_waterPhantom = [waterPhantomRunfileName, '_blockB.txt'];
        blockC_waterPhantom = [waterPhantomRunfileName, '_blockC.txt'];
        
        if ismac || isunix
            system(['cat ', blockA_0,' >> ', runfileName]);
            system(['cat ', blockA_waterPhantom,' >> ', runfileName]);
            system(['cat ', blockB_0,' >> ', runfileName]);
            system(['cat ', blockB_waterPhantom,' >> ', runfileName]);
            system(['cat ', blockC_0,' >> ', runfileName]);
            system(['cat ', blockC_waterPhantom,' >> ', runfileName]);
        elseif ispc
            system(['type ', blockA_0,' >> ', runfileName]);
            system(['type ', blockA_waterPhantom,' >> ', runfileName]);
            system(['type ', blockB_0,' >> ', runfileName]);
            system(['type ', blockB_waterPhantom,' >> ', runfileName]);
            system(['type ', blockC_0,' >> ', runfileName]);
            system(['type ', blockC_waterPhantom,' >> ', runfileName]);
        else
            disp('Platform not supported but you can concatenate the blocks to one runfile by hand.')
        end
    end
end
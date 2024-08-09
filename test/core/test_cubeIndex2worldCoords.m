function test_suite = test_cubeIndex2worldCoords

test_functions = localfunctions();

initTestSuite;

function test_empty_cubeIndex2worldCoords
    ct = helper_getTestCt();
    assertExceptionThrown(@()matRad_cubeIndex2worldCoords()); 
    assertExceptionThrown(@()matRad_cubeIndex2worldCoords([],ct));
    assertExceptionThrown(@()matRad_cubeIndex2worldCoords([100 100 100],ct));  %Out of bounds

%with ct.dicom
function test_Dicom_cubeIndex2worldCoords
    ct = helper_getTestCt();
    ct.dicomInfo.ImagePositionPatient = [-10 -10 -10];
    index = [31 32 33] ;
    expected = [83 80 86];
    assertEqual(matRad_cubeIndex2worldCoords(index , ct), expected);
    assertEqual(matRad_cubeIndex2worldCoords(sub2ind(ct.cubeDim,index(1),index(2),index(3)) , ct), expected);

    %Multiple indices
    index = [31 32 33; 33 31 32];
    expected = [83 80 86; 80 86 83];
    assertEqual(matRad_cubeIndex2worldCoords(index , ct), expected);
    assertEqual(matRad_cubeIndex2worldCoords(sub2ind(ct.cubeDim,index(:,1),index(:,2),index(:,3)) , ct), expected);
    
    %Wrong dimension error
    assertExceptionThrown(@()matRad_cubeIndex2worldCoords([1 1 1 1],ct));
    assertExceptionThrown(@()matRad_cubeIndex2worldCoords([1 1],ct));


% without ct.dicomInfo
function test_noDicom_cubeIndex2worldCoords
    ct = helper_getTestCt();

    index = [31 32 33];
    expected = [3 0 6];
    assertEqual(matRad_cubeIndex2worldCoords(index , ct), expected);
    assertEqual(matRad_cubeIndex2worldCoords(sub2ind(ct.cubeDim,index(1),index(2),index(3)) , ct), expected);

    index = [1 2 3; 3 2 1];
    expected = [-87 -90 -84; -87 -84 -90];

    assertEqual(matRad_cubeIndex2worldCoords(index , ct), expected);
    assertEqual(matRad_cubeIndex2worldCoords(sub2ind(ct.cubeDim,index(:,1),index(:,2),index(:,3)) , ct), expected);

% without ct.dicomInfo
function test_grid_cubeIndex2worldCoords
    ct = helper_getTestCt();

    ct.dimensions = ct.cubeDim;
    ct = rmfield(ct,'cubeDim');

    index = [31 32 33];
    expected = [3 0 6];
    assertEqual(matRad_cubeIndex2worldCoords(index , ct), expected);
    assertEqual(matRad_cubeIndex2worldCoords(sub2ind(ct.dimensions,index(1),index(2),index(3)) , ct), expected);

    index = [1 2 3; 3 2 1];
    expected = [-87 -90 -84; -87 -84 -90];

    assertEqual(matRad_cubeIndex2worldCoords(index , ct), expected);
    assertEqual(matRad_cubeIndex2worldCoords(sub2ind(ct.dimensions,index(:,1),index(:,2),index(:,3)) , ct), expected);
      


function ct = helper_getTestCt()
    ct.resolution.x = 3;
    ct.resolution.y = 3;
    ct.resolution.z = 3;
    ct.cubeDim = [60 60 60 ];

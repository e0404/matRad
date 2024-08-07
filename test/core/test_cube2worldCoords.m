function test_suite = test_cube2worldCoords

test_functions = localfunctions();

initTestSuite;

function test_empty_cube2worldCoords
    ct = get_testCtHelper();
    assertExceptionThrown(@()matRad_cube2worldCoords());
    assertExceptionThrown(@()matRad_cube2worldCoords([],ct));
    assertExceptionThrown(@()matRad_cube2worldCoords([100 100 100],ct));
    


%with ct.dicom
function test_Dicom_cube2worldCoords
    ct = get_testCtHelper();
    ct.dicomInfo.ImagePositionPatient = [-10 -10 -10];
    vCoord = [31 31 31] ;
    expected = [ 80 80 80];
    assertEqual(matRad_cube2worldCoords(vCoord , ct), expected);

% without ct.dicomInfo
function test_noDicom_cube2worldCoords
    ct = get_testCtHelper();
    vCoord = [31 31 31] ;
    expected = [ 0 0 0];
    assertEqual(matRad_cube2worldCoords(vCoord , ct), expected);   
    vCoord = [1 2 3];
    expected = [-90 -87 -84] ;
    assertEqual(matRad_cube2worldCoords(vCoord , ct), expected);   
      


function ct = get_testCtHelper()
    ct.resolution.x = 3;
    ct.resolution.y = 3;
    ct.resolution.z = 3;
    ct.cubeDim = [60 60 60 ];

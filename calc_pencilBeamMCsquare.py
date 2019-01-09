#!/usr/bin/env python3

import argparse
import sys, getopt
import scipy.io as sio
import numpy as np
import os.path

def main(argv):

    parser = argparse.ArgumentParser(description='interface between matRad and MCsquare for Dij computation')
    parser.add_argument('--exportCT', action='store_true', help='export CT cube')
    parser.add_argument('--bixelNb', default=-1, type=int, help='bixel number')
    args = parser.parse_args()

    if args.exportCT:
        print('matRad-MCsquare interface: Export CT')
        # EXPORT PATIENT DATA
        exportPatientData()

    if args.bixelNb > 0:
        print('matRad-MCsquare interface: Create plan')

        # EXPORT PLAN
        (nbHistories,nbThreads,seed) = exportPlan(args.bixelNb)

        # EXPORT CONFIG FILE
        exportConfig(args.bixelNb, Num_Primaries = nbHistories, Num_Threads = nbThreads, RNG_Seed = seed)

        # CALCULATE
        submitSimulation(args.bixelNb)

        # TRANSFORM OUTPUT
        transformOutput(args.bixelNb, totalNbParticles = nbHistories)


def transformOutput(bixelNb,totalNbParticles):

    import SimpleITK as sitk

    def load_itk(filename):
        # Reads the image using SimpleITK
        itkimage = sitk.ReadImage(filename)

        # Convert the image to a  numpy array first and then shuffle the dimensions to get axis in the order z,y,x
        image = sitk.GetArrayFromImage(itkimage)

        # Convert index from z,y,-x to -x,y,z
        image = np.swapaxes(image,0,2)
        # Convert index from -x,y,z to x,y,z
        image = np.flip(image,axis=0)

        # Read the origin of the image, will be used to convert the coordinates from world to voxel and vice versa.
        origin = np.array(list(itkimage.GetOrigin()))

        # Read the spacing along each dimension
        spacing = np.array(list(itkimage.GetSpacing()))

        return image, origin, spacing

    print("Reading MCsquare results")

    # MCsquare
    #   Dose: output in eV/g/proton
    #           eV/g/proton -> Gy/proton: convFactor = 1.602176e-19 * 1000.
    #   Energy: output in eV/proton
    #           eV/proton -> MeV/proton: convFactor = 1e-6

    bixelDose = {}

    tallies = [('Dose','physicalDose')]
    convFactor = [1.602176e-19 * 1000. * totalNbParticles]

    filename = 'Patient.mhd'
    [ct, origin, spacing] = load_itk(filename)
    # Consider matRad axis ordering (Y,X,Z) and rescale data
    ct = np.swapaxes(ct,0,1)

    dtypes = [ ("bixelDose",'f8',ct.shape) for (quantity,label) in tallies ]
    bixelDose['bixelDose'] = np.zeros(1, dtype=dtypes)

    for ix, (quantity,label) in enumerate(tallies):

        filename = 'Outputs_bixelNb_{}/{}.mhd'.format(bixelNb,quantity)
        [data, origin, spacing] = load_itk(filename)
        data = convFactor[ix]*data

        # Consider matRad axis ordering (Y,X,Z) and rescale data
        data = np.swapaxes(data,0,1)
        bixelDose['bixelDose'] = data

    # EXPORT 
    print("Creating MATLAB workspace with plan for visualization")
    sio.savemat('bixelDose_{}.mat'.format(bixelNb),bixelDose)

def exportPatientData():

    # IMPORT MATRAD PLAN
    print("Reading matRad workspace")
    inputfile = 'matRad_data.mat'
    data = sio.loadmat(inputfile)

    print("Exporting patient data")

    import SimpleITK as sitk

    # RSP cube
    cubeRSP = data['ct']['cube'][0,0].item()
    # matRad index order (Y,X,Z) -> MC-square (X,Y,Z)
    cubeRSP = np.swapaxes(cubeRSP,0,1)
    print("CUBE RSP SHAPE {}".format(cubeRSP.shape))
    # convert RSP to HU
    cubeHU = convertRSPtoHU(cubeRSP)
    print("CUBE HU SHAPE {}".format(cubeRSP.shape))

    # Create image cube
    cubeHU = np.flip(cubeHU,axis=0)
    cubeHU = np.swapaxes(cubeHU,0,2) # the method sitk.GetImageFromArray swap the axes z <-> x
    image    = sitk.GetImageFromArray(cubeHU)
    image.SetSpacing([data['ct']['resolution'][0,0]['x'][0,0].item(),
                      data['ct']['resolution'][0,0]['y'][0,0].item(),
                      data['ct']['resolution'][0,0]['z'][0,0].item()])

    sitk.WriteImage(image, 'Patient.mhd')
    print("IMAGE SIZE {}".format(image.GetSize()))

def convertRSPtoHU(data):

    print("Converting RSP to HU")

    # Convert density to HU using inverse conversion used in MC-square
    # Using the scanner
    ##Scanners/matRad_water

    ## ===================
    ## HU	density g/cm3
    ## ===================
    #-1050	0.0001
    #-999		0.001
    #   0		1.000
    #1000		2.000
    #1500   4.000
    array_densities = [0.,0.0001,0.001,1.000,2.000,4.000]
    array_HU = [-1050,-1050,-999,0,1000,1500]

    from scipy.interpolate import interp1d
    f_RSP2HU = interp1d(array_densities,array_HU,kind="linear")

#    return (f_RSP2HU(data)).astype(int)
    return f_RSP2HU(data)


# Create plan file
def exportPlan(bixelNb):

    # IMPORT MATRAD PLAN
    print("Reading matRad workspace")
    inputfile = 'matRad_data.mat'
    data = sio.loadmat(inputfile)

    try:
      MCparameters = data['MCparameters']
    except:
      print('MCparameters structure not available')
      sys.exit(2)

    try:
      pln = data['pln']
    except:
      print('pln structure not available')
      sys.exit(2)

    bixel = MCparameters['bixels'].item()[0,bixelNb-1]
    beamNb = int(bixel['pln']['beamNb'].item().squeeze())
    nbHistories = bixel['ncase'].item()
    seed = int(bixel['rngSeed'].item())
    nbThreads = int(MCparameters['nbThreads'].item())

    outputfile = 'MCpencilbeam_bixelNb_{}.txt'.format(bixelNb)

    # Dumb initializations
    PlanName = 'matRad_bixel'
    NumberOfFractions = 1
    FractionID = 1
    NumberOfFields = 1
    FieldsID = '###FieldsID\n{}\n'.format(1)
    TotalMetersetWeightOfAllFields = None

    nbParticles = []
    isoCenter = []

    if 'propStf' in pln.dtype.names:
      print('Trying to read iso-center from pln.propStf')
      nbOfIsoCenters = pln['propStf'][0,0]['isoCenter'].shape[0]
      for isoCenterIx in range(0,nbOfIsoCenters):
          isoCenterList = pln['propStf'][0,0]['isoCenter'].item()[isoCenterIx].tolist()
          isoCenter.append({ 'x': isoCenterList[0]-0.5*data['ct']['resolution'][0,0]['x'][0,0].item(),
                        'y': isoCenterList[1]-0.5*data['ct']['resolution'][0,0]['y'][0,0].item(),
                        'z': isoCenterList[2]-0.5*data['ct']['resolution'][0,0]['z'][0,0].item()})

    elif 'isoCenter'in pln.dtype.names:
      print('Trying to read iso-center from pln.isoCenter')
      nbOfIsoCenters = pln['isoCenter'].item().shape[0]
      for isoCenterIx in range(0,nbOfIsoCenters):
          isoCenterList = pln['isoCenter'].item()[isoCenterIx].tolist()
          isoCenter.append({ 'x': isoCenterList[0]-0.5*data['ct']['resolution'][0,0]['x'][0,0].item(),
                        'y': isoCenterList[1]-0.5*data['ct']['resolution'][0,0]['y'][0,0].item(),
                        'z': isoCenterList[2]-0.5*data['ct']['resolution'][0,0]['z'][0,0].item()})

    NumberOfFields = 1

    TotalMetersetWeightOfAllFields = 0
    for fieldIx in range(NumberOfFields):
        nbParticles.append(nbHistories)
    TotalMetersetWeightOfAllFields = np.sum(nbParticles)

    plan = """#TREATMENT-PLAN-DESCRIPTION
#PlanName
{}
#NumberOfFractions
{}
##FractionID
{}
##NumberOfFields
{}
{}
#TotalMetersetWeightOfAllFields
{}
 
""".format( PlanName,
            NumberOfFractions,
            FractionID,
            NumberOfFields,
            FieldsID,
            TotalMetersetWeightOfAllFields)

    # LOOP OVER ALL FIELDS
    for fieldIx in range(NumberOfFields):

        gantryAngles = pln['propStf'][0,0]['gantryAngles'].item().squeeze()
        couchAngles = pln['propStf'][0,0]['couchAngles'].item().squeeze()
        if type(gantryAngles) == list:
            matRad_gantry = gantryAngles[beamNb-1]
            matRad_couch = couchAngles[beamNb-1]
        else:
            matRad_gantry = gantryAngles
            matRad_couch = couchAngles

        FieldID = fieldIx+1
        FinalCumulativeMeterSetWeight = nbParticles[fieldIx]
        GantryAngle = -matRad_gantry + 180
        while GantryAngle > 360:
            GantryAngle -= 360
        PatientSupportAngle = matRad_couch
        IsocenterPositionX = isoCenter[fieldIx]['x']
        IsocenterPositionY = isoCenter[fieldIx]['y']
        IsocenterPositionZ = isoCenter[fieldIx]['z']
        NumberOfControlPoints = 1

        plan += """#FIELD-DESCRIPTION
###FieldID
{}
###FinalCumulativeMeterSetWeight
{}
###GantryAngle
{}
###PatientSupportAngle
{}
###IsocenterPosition
{}           {}           {}
###NumberOfControlPoints
{}

#SPOTS-DESCRIPTION
""".format( FieldID,
            FinalCumulativeMeterSetWeight,
            GantryAngle,
            PatientSupportAngle,
            IsocenterPositionX,
            IsocenterPositionY,
            IsocenterPositionZ,
            NumberOfControlPoints)

        # LOOP OVER ALL CONTROL POINTS
        for controlPointIx in range(NumberOfControlPoints):

            ControlPointIndex = controlPointIx + 1
            SpotTunnedID = 1
            CumulativeMetersetWeight = nbHistories
            Energy = bixel['stf']['energy'].item().squeeze()
            NbOfScannedSpots = 1

            plan += """####ControlPointIndex
{}
####SpotTunnedID
{}
####CumulativeMetersetWeight
{}
####Energy (MeV)
{}
####NbOfScannedSpots
{}
####X Y Weight
""".format( ControlPointIndex,
            SpotTunnedID,
            CumulativeMetersetWeight,
            Energy,
            NbOfScannedSpots)


            # LOOP OVER ALL CONTROL POINTS
            for SpotIx in range(NbOfScannedSpots):

                X = -1. * bixel['stf']['posY'].item().squeeze()
                Y = -1. * bixel['stf']['posX'].item().squeeze()
                Weight = nbHistories

                plan += "{} {} {}\n".format( X, Y, Weight)

    with open(outputfile,'w') as f:
        f.write(plan)
    f.close()

    return (TotalMetersetWeightOfAllFields,nbThreads,seed)

# Create config file
def exportConfig(bixelNb, Num_Primaries, Num_Threads, RNG_Seed):

    print("Exporting config file")

    CT_File = 'Patient.mhd'

    HU_Density_Conversion_File = 'Scanners/matRad_water/HU_Density_Conversion.txt'
    HU_Material_Conversion_File = 'Scanners/matRad_water/HU_Material_Conversion.txt'
    BDL_Machine_Parameter_File = 'BDL/BDL_matrad.txt'
    BDL_Plan_File = 'MCpencilbeam_bixelNb_{}.txt'.format(bixelNb)

    config = """
######################
# Configuration file #
######################

### Simulation parameters:
Num_Threads 	{}
RNG_Seed	{}
Num_Primaries 	{}
E_Cut_Pro	0.5
D_Max		0.2
Epsilon_Max	0.25
Te_Min		0.05

### Input files
CT_File 			{}
HU_Density_Conversion_File	{}
HU_Material_Conversion_File	{}
BDL_Machine_Parameter_File 	{}
BDL_Plan_File 			{}

### Physical parameters
Simulate_Nuclear_Interactions	True
Simulate_Secondary_Protons	True
Simulate_Secondary_Deuterons	True
Simulate_Secondary_Alphas	True

### 4D simulation 
4D_Mode				False
4D_Dose_Accumulation		False

### Robustness simulation
Robustness_Mode			False
Simulate_nominal_plan		True
#Systematic_Setup_Error		0.25 0.25 0.25
#Random_Setup_Error		0.1  0.1  0.1
#Systematic_Range_Error		3.0

### Beamlet simulation
Beamlet_Mode			False
Beamlet_Parallelization		False

### Output parameters
Output_Directory		Outputs_bixelNb_{}

Energy_ASCII_Output		False
Energy_MHD_Output		False
Energy_Sparse_Output		False
Dose_ASCII_Output		False
Dose_MHD_Output			True
Dose_Sparse_Output		False
LET_ASCII_Output		False
LET_MHD_Output			False
LET_Sparse_Output		False

Densities_Output		False
Materials_Output		False

Compute_DVH			False

Dose_Sparse_Threshold		0
Energy_Sparse_Threshold		0
LET_Sparse_Threshold		0

Score_PromptGammas		False
PG_LowEnergyCut 		0.0
PG_HighEnergyCut		50.0
PG_Spectrum_NumBin 		150
PG_Spectrum_Binning 		0.1

LET_Calculation_Method		StopPow

Dose_to_Water_conversion	Disabled

Dose_Segmentation			False
Density_Threshold_for_Segmentation	0.01
""".format( Num_Threads,
            RNG_Seed,
            Num_Primaries,
            CT_File,
            HU_Density_Conversion_File,
            HU_Material_Conversion_File,
            BDL_Machine_Parameter_File,
            BDL_Plan_File,
            bixelNb)

    with open('config_bixelNb_{}.txt'.format(bixelNb),'w') as f:
        f.write(config)
    f.close()

def submitSimulation(bixelNb):
    from subprocess import call
    print("Submitting jobs")
    ### call("./MCsquare_linux config_bixelNb_{}.txt".format(bixelNb),shell=True)
    call("MCSquare_windows.exe config_bixelNb_{}.txt".format(bixelNb),shell=True)

if __name__ == "__main__":
   main(sys.argv[1:])

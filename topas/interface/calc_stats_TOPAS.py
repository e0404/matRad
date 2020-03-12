#!/usr/bin/env python

import argparse
import sys, getopt
import scipy.io as sio
import topas2numpy as topas
import numpy as np
import os.path

def main(argv):

    try:
      opts, args = getopt.getopt(argv,"h:i:o:topasOutputType",["input=","output="])
    except getopt.GetoptError:
      print('calc_stats_TOPAS.py -matrad <matRad workspace> -i <inputfile> -o <outputfile> -topasOutputType <csv,bin>')
      sys.exit(2)
    parser = argparse.ArgumentParser(description='tool to calc stats for TOPAS simulations')
    parser.add_argument('-i', required=True, type=str, help='input file with MC parameters')
    parser.add_argument('-o', default='MCdata.mat', required=False, type=str, help='output file, default MCdata.mat')
    parser.add_argument('-topasOutputType', default='bin', required=False, type=str, help='output file, default bin')
    args = parser.parse_args()

    for opt, arg in opts:
      if opt == '-h':
         print('calc_stats_TOPAS.py -i <inputfile> -o <outputfile>')
         sys.exit()

    inputfile = args.i
    outputfile = args.o
    topasOutputType = args.topasOutputType

    # Dose-averaged quantities
    doseAverageQuantities = ['alpha','beta','M1','F2','m2','f3','LET']

    print('Input file: ', inputfile)
    print('Output file: ', outputfile)

    MCparam = sio.loadmat(inputfile)

    nbRuns = int(MCparam['MCheader']['nbRuns'].item().squeeze().tolist())
    nbFields = int(MCparam['MCheader']['nbFields'].item().squeeze().tolist())
    simLabel = MCparam['MCheader']['simLabel'].item().squeeze().tolist()
    tallies = [quantity.item() for quantity in MCparam['MCheader']['tallies'].item().squeeze().tolist()]
    # ensure physicalDose is read first as it is required to compute dose-average quantities
    tallies.remove('physicalDose')
    tallies = [quantity for alist in [['physicalDose'],tallies] for quantity in alist]
    nbHistories = np.array(MCparam['MCheader']['nbHistories'].item().squeeze().tolist())
    nbHistories = nbHistories.reshape(nbHistories.size)
    nbParticles = np.array(MCparam['MCheader']['nbParticles'].item().squeeze().tolist())
    nbParticles = nbParticles.reshape(nbParticles.size)
    cubeDim = tuple([int(dim) for dim in MCparam['MCheader']['cubeDim'].item().squeeze().tolist()])
    RSP = np.array(MCparam['MCheader']['RSP'].item())
    voxelDimensions = np.array(MCparam['MCheader']['voxelDimensions'].item().squeeze().tolist())
    voxel_volume_cm3 = 1.0e-3*np.product(voxelDimensions)

    resultMC = {}

    dtypeSum = [ ("MC_%s" % quantity,'f8',cubeDim) for quantity in tallies ]
    dtypeFields = [ ("MC_%s_f%d" % (quantity,ifield+1),'f8',cubeDim) for quantity in tallies for ifield in range(nbFields) ]
    dtypeFields_STD = [ ("MC_%s_f%d_STD" % (quantity,ifield+1),'f8',cubeDim) for quantity in tallies for ifield in range(nbFields) ]
    dtypes = dtypeSum + dtypeFields + dtypeFields_STD

    print(dtypes)
    resultMC['resultMC'] = np.zeros(1, dtype=dtypes)

    for quantity in tallies:
        print("Adding quantity %s" % quantity)
        quantityKey = "MC_%s" % quantity

        for iField in range(nbFields):
            print("-- Reading data for beam# %d" % (iField+1))
            quantityFieldKey = "MC_%s_f%d" % (quantity,iField+1)
            quantityFieldKey_STD = "MC_%s_f%d_STD" % (quantity,iField+1)

            data = []
            for iRun in range(nbRuns):

                filename = "simData_%s_field%d_run%d_%s.%s" % (simLabel,iField+1,iRun+1,quantity,'bin')
                if not os.path.isfile(filename):
                    filename = "simData_%s_field%d_run%d_%s.%s" % (simLabel,iField+1,iRun+1,quantity,'csv')
                if not os.path.isfile(filename):
                    raise Exception('File %s does not exists!' % filename)

                print("  -- Reading data for run# %d: %s" % (iRun+1,filename))
                binnedResult = topas.BinnedResult(filename)
                if 'Sum' in binnedResult.statistics:
                    dataRun = binnedResult.data['Sum']
                    if quantity[0:12] == 'sqrt_betaHCP':
                        dataRun = np.square(binnedResult.data['Sum'])
                    if quantity in ['beta']:
                        dataRun = np.square(binnedResult.data['Sum'])
                    if quantity in ['physicalDose',
                                    'AbsNbSmallClust',
                                    'AbsNbLargeClust',
                                    'AbsNbSmallClust_V2',
                                    'AbsNbLargeClust_V2']:
                        dataRun *= (float(nbParticles[iField])/float(nbHistories[iField]))
                        # Convert units: 10^6 clusters/voxel -> 10^6 clusters/mug
                        if quantity in ['AbsNbSmallClust',
                                        'AbsNbLargeClust',
                                        'AbsNbSmallClust_V2',
                                        'AbsNbLargeClust_V2']:
                            #voxel_volume_cm3
                            M_mug = np.multiply(RSP,1.0e+6*voxel_volume_cm3)
                            # RSP has axis ordering (Y,X,Z), but TOPAS data from run has ordering (X,Y,Z)
                            M_mug = np.swapaxes(M_mug,0,1)
                            idx = M_mug>0
                            dataRun[idx] = np.divide(dataRun[idx],M_mug[idx])
                            idx = M_mug<=0
                            dataRun[idx] = 0.

                    data.append(dataRun)
                else:
                    print("Sum not available")

            if len(data) == nbRuns:
                print("  -- Calculating stats for %d runs" % nbRuns)
                data_MEAN = np.mean(np.asarray(data),0)
                # standard deviation of the sample
                data_STD = np.std(np.asarray(data),0,ddof=1)
                # standard deviation of the mean
#                data_STD = np.std(np.asarray(data),0,ddof=1)
#                data_STD = 1./np.sqrt(float(nbRuns-1)) * data_STD

                # Consider matRad axis ordering (Y,X,Z) and rescale data
                data_MEAN = np.swapaxes(data_MEAN,0,1)
                data_STD = np.swapaxes(data_STD,0,1)

                resultMC['resultMC'][quantityFieldKey] = data_MEAN
#                resultMC['resultMC'][quantityFieldKey_STD] = np.divide(data_STD,data_MEAN)
                resultMC['resultMC'][quantityFieldKey_STD] = data_STD
                # dose-average quantities
                if quantity in doseAverageQuantities:
                    resultMC['resultMC'][quantityKey] += np.multiply(resultMC['resultMC']["MC_physicalDose_f%d" % (iField+1)],data_MEAN)
                else:
                    resultMC['resultMC'][quantityKey] += data_MEAN

        if quantity in doseAverageQuantities:
            idx = resultMC['resultMC']['MC_physicalDose'] > 0
            resultMC['resultMC'][quantityKey][idx] = np.divide(resultMC['resultMC'][quantityKey][idx],resultMC['resultMC']['MC_physicalDose'][idx])


    sio.savemat(outputfile,resultMC)


if __name__ == '__main__':
   main(sys.argv[1:])

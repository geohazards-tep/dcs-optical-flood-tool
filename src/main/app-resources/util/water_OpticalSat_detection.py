class ImageWork():
    def __init__(self,satType = ''):
        import os
        self.sep = os.sep
        if satType == '':
            self.satType = 'S2R'
        elif satType.upper() == 'S2R':
            self.satType = satType.upper()
        elif satType.upper() == 'L8R':
            self.satType = satType.upper()
        else:
            print 'Error in satTyp. Usable SatTypes are "S2R" or "L8R"'
            return
        self.noData = 999999999

    def noDataToNan(self, imageDict, npType='float32'):
        import numpy as np
        import copy
        if isinstance(imageDict, dict):
            tmp = self.setImageDict()
            tmp['array'] = copy.deepcopy((imageDict['array']).astype(npType))
            tmp['array'][tmp['array'] == self.noData] = np.nan
            self.copyImageProperty(tmp, imageDict)
            return copy.deepcopy(tmp)

    def nanToNoData(self, imageDict, npType='float32'):
        import numpy as np
        import copy
        if isinstance(imageDict, dict):
            tmp = self.setImageDict()
            tmp['array'] = copy.deepcopy((imageDict['array']).astype(npType))
            tmp['array'][np.isnan(tmp['array'])] = self.noData
            self.copyImageProperty(tmp, imageDict)
            return copy.deepcopy(tmp)

    # IMAGE READ/WRITE #

    def getGeoImage(self, path):
        '''
        Import Image data file
        :param path: Image path
        :return: GDALDataset or None
        '''
        import os
        from osgeo import gdal
        if os.path.isfile(path):
            try:
                dataSet = gdal.Open(path, gdal.GA_ReadOnly)
                return dataSet
            except:
                print 'Problem in import image:\n'
                return None
        else:
            print 'It is not a file please check path'
            return None

    def getNpTOAImage(self, path,noData=0):
        from osgeo import gdal
        import numpy as np
        import copy
        import os
        gdalImage = self.getGeoImage(path)
        gdalBand = gdalImage.GetRasterBand(1)
        transform = gdalImage.GetGeoTransform()
        projection = gdalImage.GetProjection()
        rasterArray = np.array(gdalBand.ReadAsArray())
        if 'L8R' in self.satType.upper():
            path1, file = os.path.split(path)
            bandKey = {'B1':'B01','B2':'B02', 'B3':'B03','B4':'B04','B5':'B05','B6':'B06','B7':'B07','B8':'B08','B9':'B09'}
            band = file[22:24]
            multBandKey = 'MULT_' + bandKey[band]
            addBandKey = 'ADD_' + bandKey[band]

            meta = self.metaL8(path1)
            if meta.has_key(multBandKey):
                toaRefl = np.add(
                    np.multiply(
                        meta[multBandKey],
                        rasterArray
                    ),
                    meta[addBandKey]
                )
                rasterArray = np.divide(
                    toaRefl,
                    np.sin(meta['ELEVATION'] * (np.pi / 180))
                )
        elif 'S2R' in self.satType.upper():
            rasterArray = np.multiply(rasterArray, 0.0001)
        imageFile = os.path.basename(path)
        data = self.extractData(imageFile)
        # Updata no data value in array with new value
        rasterArray[rasterArray == noData] = self.noData
        return copy.deepcopy({'array': rasterArray, 'xSize': gdalBand.XSize, 'ySize': gdalBand.YSize,
                              'gdalTransform': transform, 'gdalProjection': projection, 'yyyymmdd': data})

    def writeNpBandAsImage(self, dictionaryImage, fileOutName, path='', gdalType='', gdalDrive='GTiff'):
        from osgeo import gdal, gdal_array
        import os
        if self.is_imageDict(dictionaryImage):
            if path == '':
                path = os.getcwd()
            if gdalType == '':
                gdalType = gdal_array.NumericTypeCodeToGDALTypeCode(dictionaryImage['array'].dtype)
            else:
                gdalType = gdal_array.NumericTypeCodeToGDALTypeCode(gdalType)

            if os.path.isdir(path):
                driver = gdal.GetDriverByName('GTiff')
                output = driver.Create(path + self.sep + fileOutName, dictionaryImage['xSize'],
                                       dictionaryImage['ySize'], 1, gdalType)
                outputband = output.GetRasterBand(1)
                outputband.SetNoDataValue(self.noData)
                outputband.WriteArray(dictionaryImage['array'])
                output.SetGeoTransform(dictionaryImage['gdalTransform'])
                output.SetProjection(dictionaryImage['gdalProjection'])
                return True

    # IMAGE READ/WRITE SUPPORT #

    def setImageDict(self):
        import copy
        return copy.deepcopy({'array': None, 'xSize': None, 'ySize': None,
                              'gdalTransform': None, 'gdalProjection': None, 'yyyymmdd': None})

    def copyImageProperty(self, dictionaryImage1, dictionaryImage2):
        if self.is_imageDict(dictionaryImage2) and self.is_imageDict(dictionaryImage1):
            dictionaryImage1['xSize'] = dictionaryImage2['xSize']
            dictionaryImage1['ySize'] = dictionaryImage2['ySize']
            dictionaryImage1['gdalTransform'] = dictionaryImage2['gdalTransform']
            dictionaryImage1['gdalProjection'] = dictionaryImage2['gdalProjection']
            dictionaryImage1['yyyymmdd'] = dictionaryImage2['yyyymmdd']

    # INPUT FUNCTION #

    def buildeVirtualStackImage(self, path):
        import os
        import copy
        from osgeo import gdal
        store = 0
        storeImage =[]
        storeOutName = []
        noJamp = True
        nomeOut = True
        print "buildVirtualStackImage path: ", path
	if self.satType.upper() == 'S2R':
            bands = ['B02','B03','B04','B08','B11','B12']
            attribute0 = 'S2'
            attribute1 = '.JP2'
        elif self.satType.upper() == 'L8R':
            bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
            attribute0 = 'LC8'
            attribute1 = '.TIF'
        bandIn = []
        for root, dirs, files in os.walk(path):
            allImageFiles = [x for x in files if attribute0 in x.upper() and x.upper().endswith(attribute1)]
            for band in bands:
                selectImageFile = [x for x in allImageFiles if band in x.upper()]
                if not selectImageFile == []:
                    bandIn.append(self.getNpTOAImage(root + self.sep +selectImageFile[0]))
                    if nomeOut:
                        outNameImage = selectImageFile[0].split('_')
                        outNameImage[-1] = 'WaterMask.tif'
                        outNameImage = '_'.join(outNameImage)
                        nomeOut = False
                else:
                    noJamp = False
            if  noJamp:
                storeOutName.append(outNameImage)
                store += 1
            noJamp = True
            nomeOut = True
        if store == 1:
            return (copy.deepcopy(bandIn),copy.deepcopy(storeOutName))
        else:
            return (None,None)

    # SUPPORT #

    def metaL8(self, path):
        import os
        import copy
        if os.path.isdir(path):
            metaFile = [f for f in os.listdir(path)
                        if ('_MTL.txt' in f) and os.path.isfile(os.path.join(path, f))]
        else:
            return
        with open(path + os.sep + metaFile[0]) as fd:
            result = dict()
            for line in fd.readlines():
                line = line.rstrip('\n')
                if len(line) == 0:
                    continue
                if any(x in line for x in ('DATE_ACQUIRED', 'ACQUISITION_DATE')):
                    result['DATE'] = line.strip().split('=')[1].strip()
                if 'METADATA_FILE_NAME' in line:
                    result['META'] = line.strip().split('=')[1].strip()
                    result['PATH'] = path
                # CLOUD INFO
                if 'CLOUD_COVER' in line:
                    result['CLOUD'] = float(line.strip().split('=')[1].strip())
                if 'CLOUD_COVER_LAND' in line:
                    result['LAND_CLOUD'] = float(line.strip().split('=')[1].strip())
                # IMAGE INFO
                if 'SUN_AZIMUTH' in line:
                    result['AZIMUTH'] = float(line.strip().split('=')[1].strip())
                if 'SUN_ELEVATION' in line:
                    result['ELEVATION'] = float(line.strip().split('=')[1].strip())
                if 'EARTH_SUN_DISTANCE' in line:
                    result['EARTH_SUN_DISTANCE'] = float(line.strip().split('=')[1].strip())
                # IMAGE REFLECTANCE CONSTANT
                if 'REFLECTANCE_MULT_BAND_1' in line:
                    result['MULT_B01'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_MULT_BAND_2' in line:
                    result['MULT_B02'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_MULT_BAND_3' in line:
                    result['MULT_B03'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_MULT_BAND_4' in line:
                    result['MULT_B04'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_MULT_BAND_5' in line:
                    result['MULT_B05'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_MULT_BAND_6' in line:
                    result['MULT_B06'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_MULT_BAND_7' in line:
                    result['MULT_B07'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_MULT_BAND_8' in line:
                    result['MULT_B08'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_MULT_BAND_9' in line:
                    result['MULT_B09'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_1' in line:
                    result['ADD_B01'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_2' in line:
                    result['ADD_B02'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_3' in line:
                    result['ADD_B03'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_4' in line:
                    result['ADD_B04'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_5' in line:
                    result['ADD_B05'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_6' in line:
                    result['ADD_B06'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_7' in line:
                    result['ADD_B07'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_8' in line:
                    result['ADD_B08'] = float(line.strip().split('=')[1].strip())
                if 'REFLECTANCE_ADD_BAND_9' in line:
                    result['ADD_B09'] = float(line.strip().split('=')[1].strip())
        fd.close()
        return copy.deepcopy(result)

    def setMaskBandInDict(self, file, pathIn):
        maskDict = self.inportDictionary[self.satType]['maskList'][0]
        if not maskDict == None:
            selectImageFile = [x for x in file if maskDict['bandNames'].upper() in x.upper() and
                               maskDict['outAttribute'].upper() in x.upper()]
            if not selectImageFile == []:
                return (True, selectImageFile[0], pathIn)
            else:
                return (False, None, None)
        else:
            return (False, None, None)

    def yyyydddToyyyymmdd(self, yyyy, ddd):
        import datetime
        yyyyddd = yyyy + ' ' + ddd
        tmp = datetime.datetime.strptime(yyyyddd, '%Y %j')
        return tmp.strftime('%Y%m%d')

    def extractData(self, imageFile):
        import datetime
        selectImageFile = imageFile

        if self.satType.upper() == 'S2R':
            dataIn = 25
            dataEnd = dataIn + 8
            return selectImageFile[dataIn:dataEnd]
        elif  self.satType.upper() == 'L8R':
            dataIn = 9
            dataEnd = dataIn + 7
            yyyy = selectImageFile[dataIn:dataIn + 4]
            ddd = selectImageFile[dataEnd - 3:dataEnd]
            return self.yyyydddToyyyymmdd(yyyy, ddd)
        else:
            print 'DataType wrong please check'
            return None

    def setTypeSat(self, nomeSensor):
        """
        Set internal satType
        :param nomeSensor: 3 char type str
        :return: True if satType is a key in bandDictionary  or False if it is not a key in bandDictionary
        """
        self.satType = nomeSensor

    def is_imageDict(self, dictionaryImage):
        if isinstance(dictionaryImage, dict):
            if dictionaryImage.has_key('array'):
                if dictionaryImage.has_key('xSize') and dictionaryImage.has_key('ySize'):
                    if dictionaryImage.has_key('gdalTransform') and dictionaryImage.has_key('gdalProjection'):
                        return True
        return False

    def is_sat(self, nomeSensor):
        """
        :param nomeSensor: 3 char type str
        :return: boolean  type True = Sat is in xml
        """
        return self.sensorDictionary.has_key(nomeSensor)

    def reducePixel(self, bandIn, exampleBand):
        from osgeo import gdal, gdalnumeric, gdal_array
        import numpy as np
        import copy
        driver = gdal.GetDriverByName('MEM')
        band = self.setImageDict()
        self.copyImageProperty(band, bandIn)
        band['array'] = copy.deepcopy(bandIn['array'])
        bandOut = self.setImageDict()
        self.copyImageProperty(bandOut, exampleBand)
        xSizeStart = band['xSize']
        ySizeStart = band['ySize']
        typeImageStart = gdal_array.NumericTypeCodeToGDALTypeCode(band['array'].dtype)
        imageIn = driver.Create('', xSizeStart, ySizeStart, 1, typeImageStart)
        bandImageIn = imageIn.GetRasterBand(1)
        bandImageIn.SetNoDataValue(self.noData)
        bandImageIn.WriteArray(band['array'])
        imageIn.SetGeoTransform(band['gdalTransform'])
        imageIn.SetProjection(band['gdalProjection'])
        # Out Image
        xSizeOut = bandOut['xSize']
        ySizeOut = bandOut['ySize']
        imageOut = driver.Create('', xSizeOut, ySizeOut, 1, typeImageStart)
        imageOut.SetGeoTransform(bandOut['gdalTransform'])
        imageOut.SetProjection(bandOut['gdalProjection'])
        gdal.ReprojectImage(imageIn, imageOut, None, None, gdal.GRA_Average)
        bandImageOut = imageOut.GetRasterBand(1)
        bandOut['array'] = np.array(bandImageOut.ReadAsArray())
        imageIn = None
        imageOut = None
        return copy.deepcopy(bandOut)

        # SUPPORT MATRIX FUNCTION #

    def max3Array(self, array1, array2, array3):
        import numpy as np
        import copy
        return copy.deepcopy(np.fmax(np.fmax(array1, array2), array3))

    def min3Array(self, array1, array2, array3):
        import numpy as np
        import copy
        return copy.deepcopy(np.fmin(np.fmin(array1, array2), array3))

    # INDEX #

    def ndvi1(self, nir, red, outPath='', x=0, storeImageDictionary=True):
        import numpy as np
        import copy
        import os
        obj1 = self.noDataToNan(nir, npType='float32')
        obj2 = self.noDataToNan(red, npType='float32')
        ndvi = self.setImageDict()
        self.copyImageProperty(ndvi, obj2)
        ndvi['array'] = np.divide(np.subtract(obj1['array'], obj2['array']), np.add(obj1['array'], obj2['array']))
        if not outPath == '':
            if os.path.isdir(outPath):
                nameFile = self.satType + '_' + self.imageStackedDictionary[x]['yyyymmdd'] + '_NDVI_%d.tif' % x
                self.writeNpBandAsImage(self.nanToNoData(ndvi, npType='float32'), nameFile, outPath, gdalType='',
                                        gdalDrive='GTiff')
        if storeImageDictionary:
            return copy.deepcopy(self.nanToNoData(ndvi, npType='float32'))
        return None

    def ndsi(self, green, swir1):
        import numpy as np
        import copy
        g = self.noDataToNan(green)
        s1 = self.noDataToNan(swir1)
        return copy.deepcopy(np.divide(np.subtract(g['array'], s1['array']), np.add(g['array'], s1['array'])))

    def aweiSh(self, blue, green, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = 2.5
        const2 = -1.5
        const3 = -0.25
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        step11 = np.multiply(const1, g['array'])
        step12 = np.add(b['array'], step11)
        step13 = np.multiply(const2, np.add(n['array'], s1['array']))
        step14 = np.add(step12, step13)
        step15 = np.multiply(const3, s2['array'])
        return np.add(step14, step15)

    def aweiNsh(self, green, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = 4.
        const2 = 0.25
        const3 = 2.75
        g = self.noDataToNan(green)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        step11 = np.subtract(g['array'], s1['array'])
        step12 = np.multiply(const1, step11)
        step13 = np.multiply(const2, n['array'])
        step14 = np.multiply(const3, s2['array'])
        step15 = np.add(step13, step14)
        return np.subtract(step12, step15)

    def wetness(self, blue, green, red, nir, swir1, swir2):
        import numpy as np
        import copy
        constB = 66.96
        constG = 53.55
        constR = 23.61
        constN = 16.72
        constS1 = 194.53
        constS2 = 137.19
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        bx = np.multiply(constB, b['array'])
        gx = np.multiply(constG, g['array'])
        rx = np.multiply(constR, r['array'])
        nx = np.multiply(constN, n['array'])
        s1x = np.multiply(constS1, s1['array'])
        s2x = np.multiply(constS2, s2['array'])
        step11 = np.add(bx, gx)
        step12 = np.add(rx, nx)
        step13 = np.add(step11, step12)
        step21 = np.subtract(s1x, s2x)
        return copy.deepcopy(np.subtract(step13, step21))

        # ANALISYS #

    def saturation(self, green, red, nir):
        import numpy as np
        import copy
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        maxGRN = self.max3Array(g['array'], r['array'], n['array'])
        return copy.deepcopy \
                (
                np.divide
                    (
                    np.subtract(maxGRN, self.min3Array(g['array'], r['array'], n['array'])), maxGRN
                )
            )

    def waterDetection(self,inPath, outPath, wetLimitCon=[]):
        '''
        :param satType: 'S2R' or 'L8R' with R = Raw Data
        :param inPath: image Landsat8 or Sentinel2 path
        :param outPath: path to write out image
        :param wetLimitCon: list [ wetCostant1,wetCostant1,aweiShCostant,aweiNshCostant]
                for L8R = [5,5,0.14,0.14]
                for S2R = [8,8,0.25,0.20]
        :return: tif Image water mask
        '''
        import numpy as np
        import os

        if wetLimitCon == [] or len(wetLimitCon) == 4:
            # Constant definition
            if self.satType == 'S2R' and wetLimitCon == []:
                type = 'S2 Standard Constant'
                wetCostant1 = 8
                wetCostant2 = 8
                aweiShCostant = 0.25
                aweiNshCostant = 0.20
            elif self.satType == 'S2R' and len(wetLimitCon) == 4:
                type = 'S2 User Constant'
                wetCostant1 = wetLimitCon[0]
                wetCostant2 = wetLimitCon[1]
                aweiShCostant = wetLimitCon[2]
                aweiNshCostant = wetLimitCon[3]
            elif self.satType == 'L8R' and wetLimitCon == []:
                type = 'L8 Standard Constant'
                wetCostant1 = 5
                wetCostant2 = 5
                aweiShCostant = 0.14
                aweiNshCostant = 0.14
            elif self.satType == 'L8R' and len(wetLimitCon) == 4:
                type = 'L8 User Constant'
                wetCostant1 = wetLimitCon[0]
                wetCostant2 = wetLimitCon[1]
                aweiShCostant = wetLimitCon[2]
                aweiNshCostant = wetLimitCon[3]
            else:
                return None
        else:
            return None

        if os.path.isfile(inPath):
            print '\n--- START ---\n'

           # with open(inPath) as f:
           #     pathImages = f.readlines()
           #     pathImages = [x.strip() for x in pathImages]
           #     x = 0
           #     numIm = len(pathImages)
        print '\n'+type + ' : '
        print 'wetCostant1 = ', wetCostant1
        print 'wetCostant2 = ', wetCostant2
        print 'aweiShCostant = ', aweiShCostant
        print 'aweiNshCostant = ', aweiNshCostant
        #for pathImage in pathImages:
        pathImage = inPath
	imageIn = self.buildeVirtualStackImage(pathImage)
        image = imageIn[0]
        imageOutput = imageIn[1]
                  
        if not(image is None):
        	x+=1
                outNameImage = imageOutput[0]
                print '\n\nSTART INPORT IMAGE NUM %d ON %d \n\n' % (x,numIm)
                print 'Name outFile : ', outNameImage
                nir = image[3]
                print '\n\nFINE INPORT NIR', nir['array'].shape
                blue = image[0]
                print 'FINE INPORT BLUE', blue['array'].shape
                if not (blue['array'].shape == nir['array'].shape):
                    blue = self.reducePixel(blue, nir)
                    print 'Resize BLUE', blue['array'].shape
                green = image[1]
                print 'FINE INPORT GREEN', green['array'].shape
                if not (green['array'].shape == nir['array'].shape):
                    green = self.reducePixel(green, nir)
                    print 'Resize GREEN', green['array'].shape
                red = image[2]
                print 'FINE INPORT RED', red['array'].shape
                if not (red['array'].shape == nir['array'].shape):
                     red = self.reducePixel(red, nir)
                     print 'Resize RED', red['array'].shape
                swir1 = image[4]
                print 'FINE INPORT SWIR1', swir1['array'].shape
                if not (swir1['array'].shape == nir['array'].shape):
                     swir1 = self.reducePixel(swir1, nir)
                     print 'Resize SWIR1', swir1['array'].shape
                swir2 = image[5]
                print 'FINE INPORT SWIR2', swir2['array'].shape
                if not (swir2['array'].shape == nir['array'].shape):
                     swir2 = self.reducePixel(swir2, nir)
                     print 'Resize SWIR2', swir2['array'].shape
                print 'FINE INPORT IMAGE\n\n'
                print 'NDVI START'
                ndviImage = self.ndvi1(nir, red)
                print 'WATER INDEX START'
                test_aweiSh = self.setImageDict()
                self.copyImageProperty(test_aweiSh, red)
                test_aweiSh['array'] = self.aweiSh(blue, green, nir, swir1, swir2)
                test_aweiNsh = self.setImageDict()
                self.copyImageProperty(test_aweiNsh, red)
                test_aweiNsh['array'] = self.aweiNsh(green, nir, swir1, swir2)
                test_wetness = self.setImageDict()
                self.copyImageProperty(test_wetness, red)
                test_wetness['array'] = self.wetness(blue, green, red, nir, swir1, swir2)
                isWetWithAwai = (test_aweiSh['array'] > aweiShCostant) & (test_aweiNsh['array'] > aweiNshCostant)
                # Set renge pixel for CAT1
                rengeLow = -1
                rengeTop = 0
                # SET  CAT1 arrayt
                print '\nCAT1'
                print 'Step 1'
                cat1 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat1 & self.isSnowshape(blue, green, red, nir, swir1)] = 300
                print 'Step 2'
                cat1 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat1 & self.isDWAT1(blue, green, red, nir, swir1, swir2)] = 100
                print 'Step 3'
                cat1 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat1 & (self.isSWAT1(blue, red, nir, swir1, swir2) & isWetWithAwai)] = 100
                print 'Step 4'
                cat1 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat1 & self.isCl1(blue, green, red, nir)] = 200
                print 'Step 5'
                cat1 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat1 & ((test_wetness['array'] > wetCostant1) & isWetWithAwai)] = 100
                print 'Step 6'
                cat1 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat1] = 500
                print 'End CAT1'
                # Set renge pixel for CAT2
                rengeLow = 0
                rengeTop = 0.45
                # SET  CAT2 array
                print '\nCAT2'
                print 'Step 1'
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2 & self.isSnowshape(blue, green, red, nir, swir1)] = 300
                print 'Step 2'
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2 & self.isCl1(blue, green, red, nir)] = 200
                print 'Step 3'
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2 & (self.isSWAT2(blue, green, red, nir, swir1, swir2) & isWetWithAwai)] = 100
                print 'Step 4'
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2 & self.isCloudshape(blue, green, red, nir, swir1, swir2)] = 200
                print 'Step 5'
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2 & (self.isCl2(blue, green, red, nir) & (ndviImage['array'] > 0.40))] = 200
                print'Step 6'
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2 & (self.isCl3(blue, green, red, nir) & (ndviImage['array'] > 0.35))] = 200
                print'Step 7'
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2 & self.isSadow1(blue, green, red, nir)] = 500
                print'Step 8'
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2 & ((test_wetness['array'] > wetCostant2) & isWetWithAwai)] = 100
                cat2 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat2] = 0
                # Set renge pixel for CAT3
                rengeLow = 0.45
                rengeTop = 1
                cat3 = np.logical_and(ndviImage['array'] >= rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat3] = 0
                # BuidMask
                rengeLow = 100
                rengeTop = 500
                cat3 = np.logical_and(ndviImage['array'] > rengeLow, ndviImage['array'] <= rengeTop)
                ndviImage['array'][cat3] = 0
                self.writeNpBandAsImage(ndviImage, outNameImage, outPath, np.int16)
                nir = None
                blue = None
                green = None
                red = None
                swir1 = None
                swir2 = None
                n = None
                b = None
                g = None
                r = None
                s1 = None
                s2 = None
                print '\n WORKED %d IMAGES ON %d' % (x, numIm)

    # BOOLEANE FUNCTION #

    def isWatershape1(self, blue, green, red, nir, swir1):
        import numpy as np
        import copy
        const1 = -0.2
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        difBG = np.subtract(b['array'], g['array'])
        step11 = difBG > const1
        step12 = step11 & (difBG >= g['array'])
        step13 = step12 & (difBG >= r['array'])
        step14 = step13 & (difBG >= n['array'])
        return copy.deepcopy(step14 & (difBG >= s1['array']))

    def isWatershape2(self, blue, green, red, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = 1.3
        const2 = 0.12
        const3 = 0.039
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        redM = np.multiply(const1, r['array'])
        step11 = b['array'] >= g['array']
        step12 = step11 & (g['array'] >= r['array'])
        step13 = step12 & (r['array'] <= n['array'])
        step14 = step13 & (n['array'] < redM)
        step15 = step14 & (redM < const2)
        step16 = step15 & (const2 > s1['array'])
        step17 = step16 & (s1['array'] > s2['array'])
        step21 = const3 < n['array']
        step22 = step21 & (n['array'] < g['array'])
        return copy.deepcopy(step17 & step22)

    def isSnowshape(self, blue, green, red, nir, swir1):
        import numpy as np
        import copy
        const1 = 0.30
        const2 = 0.65
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        max = self.min3Array(b['array'], g['array'], r['array'])
        max1 = np.fmin(max, n['array']) > const1
        max2 = self.ndsi(green, swir1) > const2
        return copy.deepcopy(max1 & max2)

    def isGrowing14(self, blue, green, red, nir):
        import numpy as np
        import copy
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        step11 = b['array'] < g['array']
        step12 = step11 & (g['array'] < r['array'])
        return copy.deepcopy(step12 & (r['array'] < n['array']))

    def isGrowing15(self, blue, green, red, nir, swir1):
        import numpy as np
        import copy
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        step11 = b['array'] < g['array']
        step12 = step11 & (g['array'] < r['array'])
        step13 = step12 & (r['array'] < n['array'])
        return copy.deepcopy(step13 & (n['array'] < s1['array']))

    def isBrightsoil(self, blue, green, red, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = 0.27
        const2 = 0.038
        b = self.noDataToNan(blue)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        step11 = b['array'] < const1
        step12 = step11 & self.isGrowing15(blue, green, red, nir, swir1)
        step21 = step11
        step22 = step21 & self.isGrowing14(blue, green, red, nir)
        step23 = step22 & (np.subtract(s1['array'], s2['array']) > const2)
        return copy.deepcopy(np.logical_or(step12, step23))

    def isCloudshape1(self, blue, green, red, nir, swir1):
        import numpy as np
        import copy
        const1 = 0.17
        const2 = 0.30
        const3 = 1.3
        const4 = 0.95
        const5 = 0.65
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        step11 = self.min3Array(b['array'], g['array'], r['array']) > const1
        step12 = step11 & (np.fmax(self.max3Array(g['array'], b['array'], r['array']), n['array']) > const2)
        step13 = step12 & (np.divide(n['array'], r['array']) >= const3)
        step14 = step13 & (np.divide(n['array'], g['array']) >= const3)
        step15 = step14 & (np.divide(n['array'], s1['array']) >= const4)
        step16 = step15 & (s1['array'] > self.min3Array(b['array'], g['array'], r['array']))
        return copy.deepcopy(step16 & (self.ndsi(green, swir1) < const5))

    def isCloudshape2(self, blue, green, red, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = 0.47
        const2 = 0.37
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        max1 = self.max3Array(b['array'], g['array'], r['array'])
        max2 = self.max3Array(n['array'], s1['array'], s2['array'])
        step11 = np.fmax(max1, max2) > const1
        step12 = step11 & (np.fmin(self.min3Array(b['array'], g['array'], r['array']), n['array']) > const2)
        step13 = step12 & (self.isSnowshape(blue, green, red, nir, swir1) is False)
        return copy.deepcopy(step13 & (self.isBrightsoil(blue, green, red, nir, swir1, swir2) is False))

    def isCloudshape3(self, blue, green, red, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = 0.21
        const2 = 0.2
        const3 = 0.4
        const4 = 0.35
        const5 = -0.3
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        tmpSaturation = self.saturation(green, red, nir)
        step11 = self.min3Array(b['array'], g['array'], r['array']) > const1
        step12 = step11 & (s1['array'] > self.min3Array(b['array'], g['array'], r['array']))
        step21 = tmpSaturation >= const2
        step22 = step21 & (tmpSaturation <= const3)
        step13 = step12 & step22
        step14 = step13 & (self.max3Array(g['array'], r['array'], n['array']) >= const4)
        step15 = step14 & (self.isSnowshape(blue, green, red, nir, swir1) is False)
        step16 = step15 & (self.ndsi(green, swir1) > const5)
        return copy.deepcopy(step16 & (self.isBrightsoil(blue, green, red, nir, swir1, swir2) is False))

    def isCloudshape(self, blue, green, red, nir, swir1, swir2):
        import copy
        import numpy as np
        return copy.deepcopy(
            np.logical_or
                (
                np.logical_or
                    (
                    self.isCloudshape1(blue, green, red, nir, swir1),
                    self.isCloudshape2(blue, green, red, nir, swir1, swir2)
                ),
                self.isCloudshape3(blue, green, red, nir, swir1, swir2)
            )
        )

    def isDWAT1(self, blue, green, red, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = -0.2
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        difBG = np.subtract(b['array'], g['array'])
        step11 = difBG > const1
        step12 = step11 & (difBG >= g['array'])
        step13 = step12 & (difBG >= r['array'])
        step14 = step13 & (difBG >= n['array'])
        waterShape = step14 & (difBG >= s1['array'])
        # end WaterShape start STEP 1
        const2 = 0.078
        const3 = 0.04
        const4 = 0.12
        step21 = waterShape & (b['array'] >= const2)
        step22 = step21 & (g['array'] >= const3)
        step23 = step22 & (g['array'] <= const4)
        step24 = step23 & (np.fmax(s1['array'], s2['array']) <= const3)
        return copy.deepcopy(step24)

    def isSWAT1(self, blue, red, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = 0.04
        const2 = 0.19
        const3 = 0.078
        b = self.noDataToNan(blue)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        step11 = r['array'] >= self.max3Array(n['array'], s1['array'], s2['array'])
        step12 = step11 & (r['array'] >= const1)
        step13 = step12 & (r['array'] <= const2)
        step14 = step13 & (b['array'] > const3)
        return copy.deepcopy(step14 & (np.fmax(s1['array'], s2['array']) < const1))

    def isSWAT2(self, blue, green, red, nir, swir1, swir2):
        import numpy as np
        import copy
        const1 = 1.3
        const2 = 0.12
        const3 = 0.039
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        s1 = self.noDataToNan(swir1)
        s2 = self.noDataToNan(swir2)
        redM = np.multiply(const1, r['array'])
        step11 = b['array'] >= g['array']
        step12 = step11 & (g['array'] >= r['array'])
        step13 = step12 & (r['array'] <= n['array'])
        step14 = step13 & (n['array'] < redM)
        step15 = step14 & (redM < const2)
        step16 = step15 & (const2 > s1['array'])
        step17 = step16 & (s1['array'] > s2['array'])
        step21 = const3 < n['array']
        step22 = step21 & (n['array'] < g['array'])
        isWatershape2 = step17 & step22
        const4 = 0.078
        const5 = 0.058
        step31 = b['array'] > const4
        step32 = step31 & (np.fmax(s1['array'], s2['array']) < const5)
        return copy.deepcopy(isWatershape2 & step32)

    def isCl1(self, blue, green, red, nir):
        import numpy as np
        import copy
        const1 = 0.94
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        step11 = b['array'] > const1
        step12 = step11 & (g['array'] > const1)
        step13 = step12 & (r['array'] > const1)
        return copy.deepcopy(step13 & (n['array'] > const1))

    def isCl2(self, blue, green, red, nir):
        import numpy as np
        import copy
        const1 = 0.254
        const2 = 0.165
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        step11 = b['array'] > g['array']
        step12 = step11 & (b['array'] > r['array'])
        step13 = step12 & (n['array'] > const1)
        return copy.deepcopy(step13 & (b['array'] > const2))

    def isCl3(self, blue, green, red, nir):
        import numpy as np
        import copy
        const1 = 0.27
        const2 = 0.21
        const3 = 0.1
        const4 = 0.35
        b = self.noDataToNan(blue)
        g = self.noDataToNan(green)
        r = self.noDataToNan(red)
        n = self.noDataToNan(nir)
        step11 = b['array'] > g['array']
        step12 = step11 & (b['array'] > const1)
        step13 = step12 & (g['array'] > const2)
        step14 = step13 & (np.abs(np.subtract(r['array'], g['array'])) <= const3)
        return copy.deepcopy(step14 & (n['array'] > const4))

    def isSadow1(self, blue, green, red, nir):
        import numpy as np
        import copy
        const1 = 0.13
        const2 = 0.05
        const3 = -0.04
        b = self.noDataToNan(blue)
        r = self.noDataToNan(red)
        g = self.noDataToNan(green)
        n = self.noDataToNan(nir)
        step11 = b['array'] < const1
        step12 = step11 & (b['array'] > g['array'])
        step13 = step12 & (g['array'] > r['array'])
        step14 = step13 & (r['array'] < const2)
        return copy.deepcopy(step14 & (np.subtract(b['array'], n['array']) < const3))





def water_OpticalSat_detection_body(image_list=None, type_sat=None, window=None, outdir=None, smallest_flood_pixels=None, proc_param=None):
    
    cc=window.split()
    if not cc==[]:
        window=[(int(cc[0]),int(cc[1])),(int(cc[2]),int(cc[3]))]
    else:
        window=cc

    cc=proc_param.split()
    if not cc==[]:
        proc_param=[float(cc[0]),float(cc[1]),float(cc[2]),float(cc[3])]
    else:
        proc_param=cc

    print type_sat
    satData = ImageWork(type_sat)
    satData.waterDetection(image_list, outdir, proc_param)





def water_OpticalSat_detection_main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--image_list", default="", type=str, help="optical image list")
    parser.add_argument("-t", "--type_sat", default="S2R", type=str, help="Sensor type: it can be S2R for Sentinel2 and L8R for Landsat-8")
    parser.add_argument("-w", "--window", default="", type=str, help="AOI for analysis")
    parser.add_argument("-o", "--outdir", default="", type=str, help="output directory")
    parser.add_argument("-s", "--smallest_flood_pixels", default="9", type=int, help="minimum size for flood areas")
    parser.add_argument("-p", "--proc_param", default="", type=str, help="minimum size for flood areas")


    kwargs = vars(parser.parse_args())
    water_OpticalSat_detection_body(**kwargs)







if __name__ == '__main__':
    water_OpticalSat_detection_main()


##riga di lancio:  water_OpticalSat_detection.py --image_list lista_immagini.txt --type_sat 'S2R' --window 'xmin ymin xdim ydim' --outdir=./ --proc_param='8 8 0.20 0.25'

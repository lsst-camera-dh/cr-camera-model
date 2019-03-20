"""
Makes a raft's worth of image files based on a user-specified single-sensor
image file, the raft designator, and the sensor type (ITL or E2V, which
probably could be guessable from the single-sensor image).  The test type
probably also needs to be specified (for the file naming)
The data in the nine output images are the same as in the input image.  What
is new is the raft coordinate-related header keywords.  These are derived
from the specification in LCA-13501.
"""

from __future__ import print_function, absolute_import, division

import sys
import os
import fnmatch
import glob

import numpy as np
import astropy.io.fits as fits

import siteUtils
from datacat.error import DcClientException
from DataCatalog import DataCatalog
import camera_components

# Set the Data Catalog folder.
ROOT_FOLDER = camera_components.ROOT_FOLDER

def make_datacat_path(**kwargs):
    """ Build the data catalog path for a particular test on a particular sensor

    Looking at LSSTTD-690 it appear that with run numbers the best we can do is
    <root_folder>/<sensor_type>/<sensor_id>.
    """
    return os.path.join(kwargs.get('root_folder', ROOT_FOLDER),
                        kwargs['sensor_type'],
                        kwargs['sensor_id'])


def sort_unique(filelist):
    """ Remove duplicate files from re-running particular scripts

    This just keeps the file with the highest JOB ID
    """
    # Build a dictionary of dictionaries mapping base_id : job_id : filename
    sort_dict = {}
    for filename in filelist:
        fbase = os.path.splitext(os.path.basename(filename))[0]
        job_id = fbase.split('_')[-1]
        fid = fbase.replace('_%s'%job_id, '')
        try:
            sort_dict[fid].update({job_id:filename})
        except KeyError:
            sort_dict[fid] = {job_id:filename}

    # For each dictionary pull out just the filename associated to the latest job_id
    ret_list = []
    for value in sort_dict.values():
        keys2 = value.keys()
        keys2.sort()
        ret_list += [value[keys2[-1]]]

    return ret_list


def make_outfile_path(outpath='.', slot_name=None, file_string=None,
                      jobname=None):
    """
    Build the staged output path for a file produced by a harnessed
    job.  This means that the path will be relative to the cwd.  For
    raft testing of sensors, the output file path will be of the form

    ./<slot_name>/<file_name>

    where file_name is the basename of file_string with the jobname
    inserted before the timestamp, thereby making the job name part of
    the original run number.

    If slot_name is None, then the output file path is

    ./<file_name>
    """
    if file_string is None:
        raise RuntimeError("fake_raft.make_outfile_path: file_string must not be None")
    if jobname is None:
        raise RuntimeError("fake_raft.make_outfile_path: jobname must not be None")
    fname, ext = os.path.splitext(file_string)
    tokens = fname.split('_')
    tokens.insert(-1, jobname)
    fname = '_'.join(tokens)
    fname += ext
    if slot_name is not None:
        outfile_path = os.path.join(outpath, slot_name, fname)
    else:
        outfile_path = os.path.join(outpath, fname)
    return outfile_path


def get_file_suffix(filepath):
    """ Get the last part of the file path, basically everything after the sensor_id.
    e.g., '_mean_bias_25.fits' or '_eotest_results.fits'
    """
    basename = os.path.basename(filepath)
    tokens = basename.split('_')
    retval = ""
    for tok in tokens[1:]:
        retval += '_'
        retval += tok
    return retval


def get_template_files(root_folder, sensor_type, sensor_id, process_name, **kwargs):
    """ Get files that serve as templates for a particular process_type

    Parameters
    ----------
    root_folder : str
        Top level data catalog folder for search
    sensor_type : str
        Type of sensor, e.g., 'E2V-CCD'
    sensor_id : str
        Name of the sensor, e.g., 'E2V-CCD250-104'
    process_name : str
        Name of the eTraveler process, e.g., 'fe55_acq' or 'dark_acq'

    Keyword arguments
    -----------
    image_type : str or None
        Type of images to find
    test_type : str or None
        Type of test to find images for
    pattern : str
        Regular expression specifying which files to get / copy
    site : str
        Specifies data catalog database to access
    sort : bool
        Sort the file names before returning them
    test_version : str
        Version of the test process to search for

    Returns
    ----------
    file_list : list
        List of file names for files that match the process_name and sensor_id
    """
    pattern = kwargs.get('pattern', '*.fits')
    image_type = kwargs.get('image_type', None)
    test_type = kwargs.get('test_type', None)

    try:
        folder = os.environ['LCATR_DATACATALOG_FOLDER']
    except KeyError:
        folder = make_datacat_path(root_folder=root_folder, sensor_type=sensor_type,
                                   sensor_id=sensor_id)

    datacat = DataCatalog(folder=folder, site=kwargs.get('site', 'slac.lca.archive'))
    query = '&&'.join(('LSST_NUM=="%s"' % sensor_id,
                       'ProcessName=="%s"' % process_name))
    if image_type is not None:
        query += '&& IMGTYPE == "%s"' % image_type
    if test_type is not None:
        query += '&& TESTTYPE == "%s"' % test_type

    file_list = []
    try:
        datasets = datacat.find_datasets(query)
    except DcClientException as eobj:
        # Make the error message a bit more useful for debbuging
        msg = eobj.message + (":\nFolder = %s\n" % folder)
        msg += "Query = %s\n" % query
        sys.stderr.write(msg)
        return file_list
        #raise DcClientException(msg)

    for item in datasets.full_paths():
        if fnmatch.fnmatch(os.path.basename(item), pattern):
            file_list.append(item)

    file_list = sort_unique(file_list)
    if kwargs.get('sort', False):
        file_list = sorted(file_list)

    return file_list


def get_crosstalk_files(root_folder, sensor_id):
    """ Get files that serve as templates for a particular process_type

    Parameters
    ----------
    root_folder : str
        Top level data catalog folder for search
    sensor_id : str
        Name of the sensor, e.g., 'E2V-CCD250-104'

    Returns (spot_files, dark_files)
    ----------
    spot_files : list
        List of files from spot exposures for this sensor

    dark_files : list
        List of files from darks for this sensor
    """
    folder = os.path.join(root_folder, sensor_id.replace('-', '_'), 'spot')
    globstring_spot = '%s/*/%s_spot_spot_multi*.fits'%(folder, sensor_id.replace('-', '_'))
    globstring_dark = '%s/*/%s_dark_dark*.fits'%(folder, sensor_id.replace('-', '_'))

    spot_files = glob.glob(globstring_spot)
    dark_files = glob.glob(globstring_dark)

    return (spot_files, dark_files)


class RaftImages(object):
    '''
    Writes a raft's worth of images based on a user-supplied single-sensor
    image.

    Parameters
    ----------
    raft_id : str
        Name of the raft, e.g., 'RAFT_000'.  This is used to evaluate coordinate
        keywords for focal plane coordinates
    process_name : str
        Name of the 'test' being applied.  This probably could be extracted
        from single_sensor_file.  No analysis is implied; the string is just
        used for assigning names to output files.
    sensor_type: str
        ITL-CCD or E2V-CCD
    output_path : str
        Path to prepended to the output file names

    Attributes
    ----------
    ccd_image : astropy.io.fits.HDUList
        This is an Astropy HDUList that contains all the headers and image
        extensions of the single_sensor_file
    '''

    def __init__(self, raft_id, process_name, sensor_type, output_path):
        """
        Class constructor.
        """
        self.raft_id = raft_id
        self.process_name = process_name
        self.sensor_type = sensor_type
        self.output_path = output_path
        try:
            os.mkdir(self.output_path)
        except OSError:
            pass

    @staticmethod
    def update_primary_header(slot_name, hdu):
        """
        Update the primary header

        Parameters
        ----------
        slot_name : str
            Name of the slot within the raft
        hdu : fits.Image
            FITS image whose header is being updated
        """
        pass

    @staticmethod
    def update_image_header(slot_name, hdu):
        """
        Update the image header for one of the readout segments.  Adds raft-level
        coordinates (one set in Camera coordinates and one set rotated so the CCD
        orientation has serial direction horizontal).  Adds amplifier and rotated
        CCD coordinates as well.  See LCA-13501.

        Parameters
        ----------
        slot_name : str
            Name of the slot within the raft
        hdu : fits.Image
            FITS image whose header is being updated
        """
        # The coordinate keyword values depend on the type of CCD.
        # Kind of awkward, but below the CCD type is identified by the assumed
        # values of DETSIZE for each type.  The image extension headers do not
        # include the sensor type explicitly

        hdu.header['SLOT'] = slot_name

        if hdu.header['DETSIZE'] == '[1:4096,1:4004]':
            # pixel parameters for e2v sensors
            dimv = 2002
            dimh = 512
            ccdax = 4004
            ccday = 4096
            ccdpx = 4197
            ccdpy = 4200
            gap_inx = 28
            gap_iny = 25
            gap_outx = 26.5
            gap_outy = 25
            preh = 10
        elif hdu.header['DETSIZE'] == '[1:4072,1:4000]':
            # pixel parameters for ITL sensors
            dimv = 2000
            dimh = 509
            ccdax = 4000
            ccday = 4072
            ccdpx = 4198
            ccdpy = 4198
            gap_inx = 27
            gap_iny = 27
            gap_outx = 26.0
            gap_outy = 26
            preh = 3
        else:
            raise RuntimeError("Sensor DETSIZE not recognized")

        # get the segment 'coordinates' from the extension name
        extname = hdu.header['EXTNAME']
        sx = int(extname[-2:-1])
        sy = int(extname[-1:])

        # For convenience of notation in LCA-13501 these are also defined as 'serial'
        # and 'parallel' indices, with Segment = Sp*10 + Ss
        sp = sx
        ss = sy

        # Extract the x and y location indexes in the raft from the slot name. Also
        # define the serial and parallel versions
        cx = int(slot_name[1:2])
        cy = int(slot_name[2:3])
        cp = cx
        cs = cy

        # Define the WCS and Mosaic keywords
        wcsnamea = 'AMPLIFIER'
        ctype1a = 'Seg_X   '
        ctype2a = 'Seg_Y   '
        wcsnamec = 'CCD     '
        ctype1c = 'CCD_X   '
        ctype2c = 'CCD_Y   '
        wcsnamer = 'RAFT    '
        ctype1r = 'RAFT_X  '
        ctype2r = 'RAFT_Y  '
        wcsnamef = 'FOCAL_PLANE'
        wcsnameb = 'CCD_SERPAR'
        ctype1b = 'CCD_S   '
        ctype2b = 'CCD_P   '
        wcsnameq = 'RAFT_SERPAR'
        ctype1q = 'RAFT_S  '
        ctype2q = 'RAFT_P  '

        if hdu.header['DETSIZE'] == '[1:4072,1:4000]':
            # header coordinate parameters for ITL sensors
            pc1_1a = 0
            pc1_2a = 1 - 2*sx
            pc2_1a = -1
            pc2_2a = 0
            crpix1a = 0
            crpix2a = 0
            crval1a = sx*(dimv + 1)
            crval2a = dimh + 1 - preh
            pc1_1c = 0
            pc1_2c = 1 - 2*sx
            pc2_1c = -1
            pc2_2c = 0
            crpix1c = 0
            crpix2c = 0
            crval1c = sx*(2*dimv + 1)
            crval2c = dimh + 1 + sy*dimh - preh
            pc1_1r = 0
            pc1_2r = 1 - 2*sx
            pc2_1r = -1
            pc2_2r = 0
            crpix1r = 0
            crpix2r = 0
            crval1r = (sx*(2*dimv + 1) + gap_outx + (ccdpx - ccdax)/2.
                       + cx*(2*dimv + gap_inx + ccdpx - ccdax))
            crval2r = (dimh + 1 + sy*dimh + gap_outy + (ccdpy - ccday)/2.
                       + cy*(8*dimh + gap_iny + ccdpy - ccday) - preh)
            pc1_1b = -1
            pc1_2b = 0
            pc2_1b = 0
            pc2_2b = 1 - 2*sp
            cdelt1b = 1
            cdelt2b = 1
            crpix1b = 0
            crpix2b = 0
            crval1b = (ss + 1)*dimh + 1 - preh
            crval2b = sp*(2*dimv + 1)
            pc1_1q = -1
            pc1_2q = 0
            pc2_1q = 0
            pc2_2q = 1 - 2*sp
            cdelt1q = 1
            cdelt2q = 1
            crpix1q = 0
            crpix2q = 0
            crval1q = (gap_outy + (ccdpy - ccday)/2. +\
                           cs*(8*dimh + gap_iny + ccdpy - ccday) + (ss + 1)*dimh + 1 - preh)
            crval2q = (sp*(2*dimv + 1) + gap_outx + (ccdpx - ccdax)/2.
                       + cp*(2*dimv + gap_inx + ccdpx - ccdax))
            dtm1_1 = -1
            dtm1_2 = 0
            dtm2_1 = 0
            dtm2_2 = 2*sx - 1
            dtv1 = (dimh + 1) + sy*dimh + preh
            dtv2 = (2*dimv + 1)*(1 - sx)
        elif hdu.header['DETSIZE'] == '[1:4096,1:4004]':
            # header coordinate parameters for e2v sensors
            pc1_1a = 0
            pc1_2a = 1 - 2*sx
            pc2_1a = 1 - 2*sx
            pc2_2a = 0
            crpix1a = 0
            crpix2a = 0
            crval1a = sx*(dimv + 1)
            crval2a = sx*(dimh + 1) + (2*sx - 1)*preh
            pc1_1c = 0
            pc1_2c = 1 - 2*sx
            pc2_1c = 1 - 2*sx
            pc2_2c = 0
            crpix1c = 0
            crpix2c = 0
            crval1c = sx*(2*dimv+1)
            crval2c = sx*(dimh+1) + sy*dimh + (2*sx - 1)*preh
            pc1_1r = 0
            pc1_2r = 1 - 2*sx
            pc2_1r = 1 - 2*sx
            pc2_2r = 0
            cdelt1r = 1
            cdelt2r = 1
            crpix1r = 0
            crpix2r = 0
            crval1r = (sx*(2*dimv + 1) + gap_outx + (ccdpx - ccdax)/2. +\
                           cx*(2*dimv + gap_inx + ccdpx - ccdax))
            crval2r = (sx*(dimh + 1) + sy*dimh + gap_outy +\
                           (ccdpy - ccday)/2. +\
                       cy*(8*dimh + gap_iny + ccdpy - ccday) + (2*sx - 1)*preh)
            pc1_1b = 1 - 2*sp
            pc1_2b = 0
            pc2_1b = 0
            pc2_2b = 1 - 2*sp
            cdelt1b = 1
            cdelt2b = 1
            crpix1b = 0
            crpix2b = 0
            crval1b = sp*(dimh + 1) + ss*dimh + (2*sp-1)*preh
            crval2b = sp*(2*dimv+1)
            pc1_1q = 1 - 2*sp
            pc1_2q = 0
            pc2_1q = 0
            pc2_2q = 1 - 2*sp
            cdelt1q = 1
            cdelt2q = 1
            crpix1q = 0
            crpix2q = 0
            crval1q = (gap_outy + (ccdpy - ccday)/2. +\
                       cs*(8*dimh + gap_iny + ccdpy - ccday) +\
                           sp*(dimh + 1) + ss*dimh) + (2*sp - 1)*preh
            crval2q = (sp*(2*dimv + 1) + gap_outx + (ccdpx - ccdax)/2. +\
                           cp*(2*dimv + gap_inx + ccdpx - ccdax))
            dtm1_1 = 1 - 2*sx
            dtm1_2 = 0
            dtm2_1 = 0
            dtm2_2 = 2*sx - 1
            dtv1 = (dimh+1 + 2*preh)*sx + sy*dimh - preh
            dtv2 = (2*dimv + 1)*(1 - sx)
        else:
            raise RuntimeError("Sensor DETSIZE not recognized")

        hdu.header['DTM1_1'] = dtm1_1
        hdu.header['DTM1_2'] = dtm1_2
        hdu.header['DTM2_1'] = dtm2_1
        hdu.header['DTM2_2'] = dtm2_2
        hdu.header['DTV1'] = dtv1
        hdu.header['DTV2'] = dtv2

        hdu.header['WCSNAMEA'] = wcsnamea
        hdu.header['CTYPE1A'] = ctype1a
        hdu.header['CTYPE2A'] = ctype2a
        hdu.header['CRVAL1A'] = crval1a
        hdu.header['CRVAL2A'] = crval2a
        hdu.header['PC1_1A'] = pc1_1a
        hdu.header['PC1_2A'] = pc1_2a
        hdu.header['PC2_1A'] = pc2_1a
        hdu.header['PC2_2A'] = pc2_2a
        hdu.header['CDELT1A'] = 1
        hdu.header['CDELT2A'] = 1
        hdu.header['CRPIX1A'] = crpix1a
        hdu.header['CRPIX2A'] = crpix2a
        hdu.header['WCSNAMEC'] = wcsnamec
        hdu.header['CTYPE1C'] = ctype1c
        hdu.header['CTYPE2C'] = ctype2c
        hdu.header['CRVAL1C'] = crval1c
        hdu.header['CRVAL2C'] = crval2c
        hdu.header['PC1_1C'] = pc1_1c
        hdu.header['PC1_2C'] = pc1_2c
        hdu.header['PC2_1C'] = pc2_1c
        hdu.header['PC2_2C'] = pc2_2c
        hdu.header['CDELT1C'] = 1
        hdu.header['CDELT2C'] = 1
        hdu.header['CRPIX1C'] = crpix1c
        hdu.header['CRPIX2C'] = crpix2c
        hdu.header['WCSNAMER'] = wcsnamer
        hdu.header['CTYPE1R'] = ctype1r
        hdu.header['CTYPE2R'] = ctype2r
        hdu.header['CRVAL1R'] = crval1r
        hdu.header['CRVAL2R'] = crval2r
        hdu.header['PC1_1R'] = pc1_1r
        hdu.header['PC1_2R'] = pc1_2r
        hdu.header['PC2_1R'] = pc2_1r
        hdu.header['PC2_2R'] = pc2_2r
        hdu.header['CDELT1R'] = cdelt1r
        hdu.header['CDELT2R'] = cdelt2r
        hdu.header['CRPIX1R'] = crpix1r
        hdu.header['CRPIX2R'] = crpix2r
        hdu.header['WCSNAMEF'] = wcsnamef
        hdu.header['WCSNAMEB'] = wcsnameb
        hdu.header['CTYPE1B'] = ctype1b
        hdu.header['CTYPE2B'] = ctype2b
        hdu.header['CRVAL1B'] = crval1b
        hdu.header['CRVAL2B'] = crval2b
        hdu.header['PC1_1B'] = pc1_1b
        hdu.header['PC1_2B'] = pc1_2b
        hdu.header['PC2_1B'] = pc2_1b
        hdu.header['PC2_2B'] = pc2_2b
        hdu.header['CDELT1B'] = cdelt1b
        hdu.header['CDELT2B'] = cdelt2b
        hdu.header['CRPIX1B'] = crpix1b
        hdu.header['CRPIX2B'] = crpix2b
        hdu.header['WCSNAMEQ'] = wcsnameq
        hdu.header['CTYPE1Q'] = ctype1q
        hdu.header['CTYPE2Q'] = ctype2q
        hdu.header['CRVAL1Q'] = crval1q
        hdu.header['CRVAL2Q'] = crval2q
        hdu.header['PC1_1Q'] = pc1_1q
        hdu.header['PC1_2Q'] = pc1_2q
        hdu.header['PC2_1Q'] = pc2_1q
        hdu.header['PC2_2Q'] = pc2_2q
        hdu.header['CDELT1Q'] = cdelt1q
        hdu.header['CDELT2Q'] = cdelt2q
        hdu.header['CRPIX1Q'] = crpix1q
        hdu.header['CRPIX2Q'] = crpix2q

    def write_sensor_image(self, single_sensor_file, slot_name, sensor_id, **kwargs):
        """
        Write a FITS image with pixel coordinates matching the specified
        sensor_id and raft_id.  The output file name is constructed from
        the test type, sensor type, sensor_id and raft_id.

        Parameters
        ----------
        single_sensor_file : str
            Name of the file to be copied
        slot_name:  str
            Name of the slot this sensor occupies
        sensor_id:  str
            Name of the sensor, e.g., 'E2V-CCD250-104'

        Keyword arguments
        -----------------
        clobber : bool, optional
            Flag indicating whether to overwrite an existing output file
        dry_run : bool, optional
            If true, just print output file names, but do not copy files
        job_id : str, optional
            Used to construct the output file name.  Defaults to siteUtils.getJobName()
        frame_suffix : str, optional
            Used in case we are making multiple frames for the same input.
            Defaults to None
        """
        file_suffix = get_file_suffix(single_sensor_file)

        clobber = kwargs.get('clobber', True)
        dry_run = kwargs.get('dry_run', False)
        job_id = kwargs.get('job_id', siteUtils.getJobName())
        frame_suffix = kwargs.get('frame_suffix', None)
        if frame_suffix is None:
            basename = "%s%s" % (sensor_id, file_suffix)
        else:
            basename = "%s%s%s" % (sensor_id, frame_suffix, file_suffix)
        outfilename = make_outfile_path(outpath=self.output_path,
                                        slot_name=slot_name,
                                        file_string=basename,
                                        jobname=jobname)
        outdir = os.path.dirname(outfilename)
        try:
            os.makedirs(outdir)
        except OSError:
            pass

        if dry_run:
            os.system("touch %s"% outfilename)
            return
        output = fits.open(single_sensor_file)

        RaftImages.update_primary_header(slot_name, output[0])

        for ext_num in range(1, 16):
            RaftImages.update_image_header(slot_name, output[ext_num])


        output.writeto(outfilename, clobber=clobber)
        output.close()

    def write_xtalk_sensor_dark(self,
                                single_sensor_file,
                                victim_slot_name,
                                aggresor_slot_name,
                                sensor_id,
                                **kwargs):
        """
        Write a FITS image with pixel coordinates matching the specified
        sensor_id and raft_id emulating a dark from a crosstalk test.
        The output file name is constructed from the test type, sensor type, sensor_id,
        raft_id and aggresor_slot_names.

        Parameters
        ----------
        single_sensor_file : str
            Name of the file to be copied
        victim_slot_name:  str
            Name of the slot this sensor occupies
        aggresor_slot_name:  str
            Name of the slot this sensor occupies
        sensor_id:  str
            Name of the sensor, e.g., 'E2V-CCD250-104'

        Keyword arguments
        -----------------
        clobber : bool, optional
            Flag indicating whether to overwrite an existing output file
        dry_run : bool, optional
            If true, just print output file names, but do not copy files
        job_id : str, optional
            Used to construct the output file name.  Defaults to siteUtils.getJobName()
        """
        self.write_sensor_image(single_sensor_file,
                                slot_name=victim_slot_name,
                                sensor_id=sensor_id,
                                frame_suffix="_%s"%(aggresor_slot_name),
                                **kwargs)

    def copy_aggressor_spot(self,
                            single_sensor_file,
                            victim_slot_name,
                            aggresor_slot_name,
                            sensor_id,
                            **kwargs):
        """
        Write a FITS image with pixel coordinates matching the specified
        sensor_id and raft_id emulating a multiaggressor spot from a crosstalk test.
        The output file name is constructed from the test type, sensor type, sensor_id,
        raft_id and aggresor_slot_names.

        Parameters
        ----------
        single_sensor_file : str
            Name of the file to be copied
        victim_slot_name:  str
            Name of the slot this sensor occupies
        aggresor_slot_name:  str
            Name of the slot this sensor occupies
        sensor_id:  str
            Name of the sensor, e.g., 'E2V-CCD250-104'

        Keyword arguments
        -----------------
        clobber : bool, optional
            Flag indicating whether to overwrite an existing output file
        dry_run : bool, optional
            If true, just print output file names, but do not copy files
        job_id : str, optional
            Used to construct the output file name.  Defaults to siteUtils.getJobName()
        """
        self.write_sensor_image(single_sensor_file,
                                slot_name=victim_slot_name,
                                sensor_id=sensor_id,
                                frame_suffix="_%s"%(aggresor_slot_name),
                                **kwargs)

    def make_victim_sensor_image(self,
                                 aggressor_spot_file,
                                 aggressor_dark_file,
                                 victim_dark_file,
                                 crosstalk,
                                 victim_slot_name,
                                 aggresor_slot_name,
                                 sensor_id,
                                 **kwargs):
        """
        Write a FITS image with pixel coordinates matching the specified
        sensor_id and raft_id emulating a victim sensor from a crosstalk test.
        The output file name is constructed from the test type, sensor type, sensor_id,
        raft_id and aggresor_slot_names.

        The output image is: crosstalk * (aggressor_spot - aggressor_dark ) +  victim_dark

        Parameters
        ----------
        aggressor_spot_file : str
            Name of the file with the spot data for the aggressor sensor
        aggressor_dark_file : str
            Name of the file with the dark data for the aggressor sensor
        victim_dark_file : str
            Name of the file with the dark data for the victim sensor
        crosstalk : float
            Crosstalk value.
        victim_slot_name:  str
            Name of the slot this sensor occupies
        aggresor_slot_name:  str
            Name of the slot this sensor occupies
        sensor_id:  str
            Name of the sensor, e.g., 'E2V-CCD250-104'

        Keyword arguments
        -----------------
        clobber : bool, optional
            Flag indicating whether to overwrite an existing output file
        dry_run : bool, optional
            If true, just print output file names, but do not copy files
        job_id : str, optional
            Used to construct the output file name.  Defaults to siteUtils.getJobName()
        """
        file_suffix = get_file_suffix(victim_dark_file).replace('sim_dark_dark',
                                                                'sim_spot_spot_multi')

        clobber = kwargs.get('clobber', True)
        dry_run = kwargs.get('dry_run', False)
        job_id = kwargs.get('job_id', siteUtils.getJobName())
        basename = "%s_%s%s" % (sensor_id, aggresor_slot_name, file_suffix)

        outfilename = make_outfile_path(outpath=self.output_path,
                                        slot_name=victim_slot_name,
                                        file_string=basename,
                                        job_id=job_id)
        outdir = os.path.dirname(outfilename)
        try:
            os.makedirs(outdir)
        except OSError:
            pass
        if dry_run:
            os.system("touch %s"% outfilename)
            return

        aggressor_spot = fits.open(aggressor_spot_file)
        aggressor_dark = fits.open(aggressor_dark_file)
        victim_dark = fits.open(victim_dark_file)

        output = aggressor_spot
        for amp in range(1, 17):
            victimimage = aggressor_spot[amp].data.astype(float)
            victimimage -= aggressor_dark[amp].data
            victimimage *= crosstalk
            victimimage += victim_dark[amp].data
            output[amp].data = victimimage.astype(int)

        RaftImages.update_primary_header(victim_slot_name, output[0])

        for ext_num in range(1, 17):
            RaftImages.update_image_header(victim_slot_name, output[ext_num])

        output.writeto(outfilename, clobber=clobber)
        output.close()


    def make_victim_raft_images(self,
                                aggressor_spot_files,
                                dark_files,
                                crosstalk_matrix,
                                raft,
                                **kwargs):
        """
        Write rafts worth of FITS images emulating a crosstalk test.

        This will loop over the input file, and either copy or combine them to
        emulate the crosstalk data aquisition.

        Parameters
        ----------
        aggressor_spot_files : list of str
            Names of the files with the aggressor spot data for the sensor
        dark_file : list of str
            Names of the files with the dark data for the sensors
        crosstalk_matrix : np.ndarray(9,9)
            Matrix of inter-sensor crosstalk values, expressed as a
            fraction of the intra-sensor values
        raft : `camera_components.raft` object
            Object that encapsulates the raft-to-sensor relationships
        """
        slot_names = raft.slot_names
        for i_agg, (aggressor_spot_file, aggressor_slot_name, aggressor_dark_file)\
                in enumerate(zip(aggressor_spot_files, slot_names, dark_files)):
            sys.stdout.write("Working on aggressor %s: "%slot_names[i_agg])
            sys.stdout.flush()
            for i_vic, (victim_dark_file, victim_slot_name, sensor_id)\
                    in enumerate(zip(dark_files, slot_names, raft.sensor_names)):
                sys.stdout.write('.')
                sys.stdout.flush()
                crosstalk_value = crosstalk_matrix[i_agg, i_vic]
                self.write_xtalk_sensor_dark(victim_dark_file,
                                             victim_slot_name,
                                             aggressor_slot_name,
                                             sensor_id,
                                             **kwargs)
                if i_agg == i_vic:
                    self.copy_aggressor_spot(aggressor_spot_file,
                                             victim_slot_name,
                                             aggressor_slot_name,
                                             sensor_id,
                                             **kwargs)
                else:
                    self.make_victim_sensor_image(aggressor_spot_file,
                                                  aggressor_dark_file,
                                                  victim_dark_file,
                                                  crosstalk_value,
                                                  victim_slot_name,
                                                  aggressor_slot_name,
                                                  sensor_id,
                                                  **kwargs)

            sys.stdout.write('!\n')


def copy_single_sensor_data(raft, process_name, output_path, **kwargs):
    """ Copy a single input file to the correct output location for
        each sensor in this raft.  Possibly updating FITS header
        infomation along the way.

    Parameters
    ----------
    raft : `camera_components.raft` object
        Object that encapsulates the raft-to-sensor relationships
    process_name : str
        The name of the assoicated eTraveler process, used in making
        the output file name
    output_path : str
        The prefix for the output file paths

    Keyword Arguments
    ----------
    test_type : str, optional
        Test type to copy files for
    image_type : str, optional
        Types of images copy
    pattern : str, optional
        Regular expression specifying which files to copy
    site : str
        Specifies data catalog database to access
    test_version : str
        Version of the test process to search for
    root_folder : str, defaults to 'LSST/mirror/BNL-prod/prod
        Allow overriding the top-level folder for the template file search
    process_name_out : str
        The name of the output eTraveler process, if it differs from
        process_name
    clobber : bool, optional
        Allow overwriting existing files
    dry_run : bool, optional
        If true, just print output file names, but do not copy files
    """
    kwargs = kwargs.copy()
    root_folder = kwargs.pop('root_folder', ROOT_FOLDER)
    kwargs_write = dict(clobber=kwargs.pop('clobber', False),
                        dry_run=kwargs.pop('dry_run', False),
                        process_name_out=kwargs.pop('process_name_out',
                                                    process_name))

    writer = RaftImages(raft.raft_id, process_name, raft.sensor_type,
                        output_path)

    for slot_name, sensor_id in raft.items():
        template_files = get_template_files(root_folder, raft.sensor_type,
                                            sensor_id=sensor_id.replace('_sim', ''),
                                            process_name=process_name, **kwargs)
        for fname in template_files:
            writer.write_sensor_image(fname, slot_name, sensor_id,
                                      **kwargs_write)


def simulate_crosstalk_data(raft, output_path, crosstalk_matrix, **kwargs):
    """ Simulate the data from a raft-level crosstalk test.
    This takes simulated sensor-level crosstalk data as inputs.

    Parameters
    ----------
    raft : `camera_components.raft` object
        Object that encapsulates the raft-to-sensor relationships
    output_path : str
        The prefix for the output file paths

    Keyword Arguments
    ----------
    root_folder : str, defaults to 'LSST/mirror/BNL-prod/prod
        Allow overriding the top-level folder for the template file search
    clobber : bool, optional
        Allow overwriting existing files
    dry_run : bool, optional
        If true, just print output file names, but do not copy files
    """

    kwargs = kwargs.copy()
    root_folder = kwargs.pop('root_folder', ROOT_FOLDER)
    kwargs_write = dict(clobber=kwargs.pop('clobber', False),
                        dry_run=kwargs.pop('dry_run', False),
                        process_name_out='crosstalk')

    writer = RaftImages(raft.raft_id, 'crosstalk', raft.sensor_type,
                        output_path)

    spotfiles = []
    darkfiles = []
    use_db = False

    for sensor_id in raft.sensor_names:
        if use_db:
            template_files_spot = get_template_files(root_folder, raft.sensor_type,
                                                     sensor_id=sensor_id.replace('_sim', ''),
                                                     process_name='spot', **kwargs)
            template_files_dark = get_template_files(root_folder, raft.sensor_type,
                                                     sensor_id=sensor_id.replace('_sim', ''),
                                                     process_name='dark', **kwargs)
        else:
            template_files_spot, template_files_dark = get_crosstalk_files('sensorData', sensor_id)
        spotfiles.append(template_files_spot[-1])
        darkfiles.append(template_files_dark[-1])

    writer.make_victim_raft_images(spotfiles, darkfiles, crosstalk_matrix, raft, **kwargs_write)



if __name__ == '__main__':

    USER = os.environ['USER']
    ETRAV_DB = 'Dev'

    # These are in caps to keep pylint happy
    TESTTYPE = 'FE55'
    IMGTYPE = 'BIAS'
    PROCESS_NAME_IN = 'vendorIngest'
    PROCESS_NAME_OUT = 'fe55_acq'
    PATTERN = '*.fits'
    OUTPATH = 'sensorData/raft_xtalk'
    RAFT_ID = 'LCA-10753-RSA_sim-0000'

    RAFT = camera_components.Raft.create_from_yaml("test_raft.yaml")
    #RAFT = camera_components.Raft.create_from_etrav(RAFT_ID, user=USER,
    #                                                db_name=ETRAV_DB)

    CROSSTALK = 0.1
    CROSSTALK_MATRIX = np.ones((9, 9), float)
    CROSSTALK_MATRIX *= CROSSTALK
    for i in range(9):
        CROSSTALK_MATRIX[i, i] = 1.0

    os.environ['LCATR_JOB'] = '00000'

    #copy_single_sensor_data(RAFT, PROCESS_NAME_IN, OUTPATH,
    #                        root_folder=ROOT_FOLDER, dry_run=True,
    #                        test_type=TESTTYPE, image_type=IMGTYPE,
    #                        pattern=PATTERN)

    simulate_crosstalk_data(RAFT, OUTPATH, CROSSTALK_MATRIX,
                            root_folder=ROOT_FOLDER, dry_run=True)

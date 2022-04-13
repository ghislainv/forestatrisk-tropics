#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
***************************************************************************
Download Globbiomass AGB or GSV maps			2018-09-18
***************************************************************************
"""
import os,sys
if (sys.version_info > (3, 0)):
    # Python 3 code 
    import urllib.request
else:
    # Python 2 code 
    import urllib

def myargsparse():
    import argparse

    thisprog=os.path.basename(sys.argv[0])
    
    epilog=\
    """****************************************************************
    \n*             Download Globbiomass AGB or GSV maps              *
    \n*****************************************************************
    \nUsage:
    \n{thisprog} -m AGB -o /home/user/Download
    \n{thisprog} -m GSV -o /home/user/Download
    """.format(thisprog=thisprog)
        
    help_m=\
    '''AGB or GSV'''
    help_o=\
    '''Output directory'''
  
    p = argparse.ArgumentParser(usage=None,description=epilog,prog=thisprog,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-m",required=True,help=help_m,action='store',dest='maps',choices=['AGB','GSV'],default=None)
    p.add_argument("-o",required=True,help=help_o,action='store',dest='outpath',default=None)
    args=p.parse_args()
    return args

def map_download (args):
    maps       = args.maps              
    outdir    = args.outpath               
   
    if maps == 'AGB':    
        flist={'http://globbiomass.org/wp-content/uploads/GB_Maps/N00E020_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00E060_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00E100_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00E140_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00W020_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00W060_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00W100_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N00W180_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40E020_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40E060_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40E100_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N40E140_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W020_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W060_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W100_agb.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W140_agb.zip'}
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W180_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80E020_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80E060_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80E100_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80E140_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W020_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W060_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W100_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W140_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W180_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/S40E140_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/S40W060_agb.zip',
               #'http://globbiomass.org/wp-content/uploads/GB_Maps/S40W100_agb.zip'}
    else:
        flist={'http://globbiomass.org/wp-content/uploads/GB_Maps/N00E020_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00E060_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00E100_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00E140_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00W020_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00W060_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00W100_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N00W180_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40E020_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40E060_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40E100_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40E140_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W020_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W060_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W100_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W140_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N40W180_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80E020_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80E060_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80E100_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80E140_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W020_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W060_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W100_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W140_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/N80W180_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/S40E140_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/S40W060_gsv.zip',
               'http://globbiomass.org/wp-content/uploads/GB_Maps/S40W100_gsv.zip'}

    for f in flist:
        print('Download ' + f)
        bname=os.path.basename(f)
        oname= outdir + '/' + bname
        if (sys.version_info > (3, 0)):
            urllib.request.urlretrieve(f, oname)
        else:
            urllib.urlretrieve(f, oname)

#########################################################################
def main():
    args=myargsparse()
    map_download(args)

if __name__ == '__main__':
    main()

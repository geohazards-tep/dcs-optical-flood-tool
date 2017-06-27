#!/opt/anaconda/bin/python
#Classe runSnap legge da una singola cartella i file in essa contenuti
#li ordina in modo decrescente per data e crea
#le coppie per lo start di SNAP
#infine crea il file name da associare all'output di SNAP




import subprocess
import os,sys
import cioppy
import string
sys.path.append('./util/')

import water_OpticalSat_detection
ciop = cioppy.Cioppy()




# define the exit codes - need to be better assessed
SUCCESS = 0
ERR_FAILED = 134

# add a trap to exit gracefully
def clean_exit(exit_code):
    log_level = 'INFO'
    if exit_code != SUCCESS:
        log_level = 'ERROR'

    msg = { SUCCESS: 'Download successfully concluded',
           ERR_FAILED: 'Unable to complete the download'}

    ciop.log(log_level, msg[exit_code])



def main():
    outdir=ciop.tmp_dir
    input = sys.stdin.readlines()
    print "input to the new node:",  input
    print "tmpdir: ", outdir
    try:
    	input_file = input[0][string.find(input[0], "'")+1:string.rfind(input[0],"'")].strip()    
	local_file = ciop.copy(input_file, outdir)
	print 'result of ciop.copy: ', local_file
	print "input file: ", input_file
        print "outdir: ", outdir
	#create dir with proper name
        filename = os.path.basename(local_file)
	print "filename: ", filename
        if 'tar.gz' in filename:
	    print 'tar.gz'
	    tmpdir=os.path.splitext(os.path.splitext(filename)[0])[0]
	    print tmpdir
	    extract_dir=outdir+os.sep+tmpdir
	    print "extract_dir: ", extract_dir
            subprocess.call(["ls","-l",local_file])
	    os.mkdir(extract_dir)
	    subprocess.call(["tar","xzf",local_file,"-C",extract_dir])
	elif 'zip' in filename:
	    print 'zip'
	    tmpdir=os.path.splitext(filename)[0]
	    print tmpdir
	    extract_dir=outdir+os.sep+tmpdir
	    print "extract_dir: ", extract_dir
            subprocess.call(["ls","-l",local_file])
	    os.mkdir(extract_dir)
	    subprocess.call(["unzip",local_file,"-d",extract_dir])
        elif 'tar.bz' in filename:
	    print 'tar.bz'
            tmpdir=os.path.splitext(os.path.splitext(filename)[0])[0]
            print tmpdir
	    extract_dir=outdir+os.sep+tmpdir
            print "extract_dir: ", extract_dir
	    subprocess.call(["ls","-l",local_file])
	    os.mkdir(extract_dir)
            subprocess.call(["tar","xjf",local_file,"-C",extract_dir])
    except:
 	print "flood_extraction: unexpected error...", sys.exc_info()[0]
        return 15     
    

    subprocess.call(["ls","-l",extract_dir])
    
    flood_file_result = water_OpticalSat_detection_body(image_folder=extract_dir, type_sat=None, outdir=extract_dir, smallest_flood_pixels=9, proc_param='8 8 0.20 0.25')
    #water_OpticalSat_detection --image_folder lista_immagini.txt --type_sat 'S2R' --window 'xmin ymin xdim ydim' --outdir=./ --proc_param='8 8 0.20 0.25'
    
    #ciop.publish(flood_file_result
    print "flood_file_result: ", flood_file_result
   
   
    
    
    #estrazione nome directory dove estrarre il file
    #etsrarlo
    
    #print "sys.stdin ", input
    #for input in sys.stdin:
    #print "sys.stdin ", input
    #process=subprocess.Popen(['opensearch-client',input_file,'enclosure'], stdout=subprocess.PIPE)
    #out, err=process.communicate()
    #res=ciop.copy(out,outdir, extract=False)
    #print res
    #output_file = ciop.publish(res, mode='silent', metalink=True)
    #print "output: ", output_file




try:
    main()
except SystemExit as e:
    if e.args[0]:
        clean_exit(e.args[0])
    raise
#else:
#    atexit.register(clean_exit, 0)




#ciop.publish(outdef, metalink = true)        


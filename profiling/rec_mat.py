import subprocess
import traceback

def call_matlab(fname, handler):
    process = subprocess.Popen(["""matlab -nojvm -nodisplay -nosplash -r "try_exit('{}')" """.format(fname)], shell=True,
                           stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE)
    try: 
        while True:
            output = process.stdout.read(1)

            if output == '' and process.poll() != None:
                break

            if output != '':
                handler(output)
    except Exception as e:
        process.kill()
        print traceback.print_tb(e.__traceback__) 




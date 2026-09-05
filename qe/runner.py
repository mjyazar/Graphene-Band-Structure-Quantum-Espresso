import subprocess
import time

NPROC = 4

def run(process, input_path, output_path):
    """
    Atomate the process of manually running pw.x << input_file.pwi >> output_file.pwo 
    in the terminal after creating an input file.
    """
    
    with open(input_path, mode="r") as input_file, open(output_path, mode="w") as output_file:
        
        start = time.perf_counter()
        
        calculation = subprocess.Popen(["mpirun", "-np", str(NPROC), process], stdin=input_file, stdout=output_file)

        while calculation.poll() is None:
            elapsed = int(time.perf_counter() - start)

            hours = elapsed // 3600
            minutes = (elapsed % 3600) // 60
            seconds = elapsed % 60

            print(f"\rElapsed: {hours:02d}:{minutes:02d}:{seconds:02d}", end="", flush=True)

            time.sleep(1)

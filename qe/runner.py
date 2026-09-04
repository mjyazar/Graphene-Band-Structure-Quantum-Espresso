import subprocess   


def run(process, input_path, output_path):
    """
    Atomate the process of manually running pw.x << input_file.pwi >> output_file.pwo 
    in the terminal after creating an input file.
    """
    
    with open(input_path, mode="r") as input_file, open(output_path, mode="w") as output_file:
        subprocess.run([process], stdin=input_file, stdout=output_file, check=True)

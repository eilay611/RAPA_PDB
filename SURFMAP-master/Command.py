import subprocess
import sys

class Command:

	@staticmethod	
	def do_command(command,allowed_return_codes=[],verbose=True):
		if verbose:
			print ("Running : " + command)
		try:
			result = subprocess.check_output(command,stderr=subprocess.STDOUT,shell=True)
                        #subprocess.run(command)
		except subprocess.CalledProcessError as e:
			print ("Returned with error(%d) : %s" % (e.returncode,e.output))
	
			if (e.returncode in allowed_return_codes):
				return
			else:
				raise e

	
	

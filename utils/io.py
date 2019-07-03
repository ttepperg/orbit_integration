'''
Author:	Thorsten Tepper Garcia
Date:	02/07/2019
'''

def get_input():
	import sys
	if len(sys.argv) < 2:
		print("\nUSAGE:")
		print("{} <input parameter file>\n".format(sys.argv[0]))
		exit()
	else:
		ics_file = sys.argv[1]
		ext = ".py"
		if ext in ics_file:
			print("\nRemoved unnecessary extension {} in input parameter file.".format(ext))
			return ics_file.replace(ext,'') 	# remove file extension if present
		else:
			return ics_file


# def write_table(file = None):
# 	'''dumps code results to ascii file'''
# 	if file is None:
# 		raise ValueError("output file name not provided in write_table")
# 	else:
# 		f = open(file, 'wt')
# 		f.write("test\n")
# 		f.close()

	
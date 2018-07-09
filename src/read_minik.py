#!usr/bin/env python3

#import pandas as pd
import logging



# 34 columns come after the desired data
# minik_columns = ["event_num", "timestamp", "???", "azimuth", "zenith", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34"]

# df = pd.DataFrame()
# df_events = pd.DataFrame(columns = minik_columns)

# df = pd.read_csv(in_file, header = None, sep = "\t", names = minik_columns)

"""
# Allows the creation of several loggers
file_formatter = logging.Formatter("%(asctime)s: %(name)s: %(levelname)-8s %(message)s")
console_formatter = logging.Formatter("%(name)s: %(levelname)-8s %(message)s")
def setup_logger(name, log_file, consol, level = logging.INFO):

    handler = logging.FileHandler(directory + "/analysis/" + log_file)        
    handler.setFormatter(file_formatter)
    
    logger = logging.getLogger(name)    
    logger.setLevel(level)
    logger.addHandler(handler)

    if consol == True:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)    
        console.setFormatter(console_formatter)
        logger.addHandler(console)    
    
    return logger


minik_logger = setup_logger("minik_logger", "minik_log.txt", consol = True)
"""


# Used to piece together multiple lines which makeup an event
#def read_event():
#	return

in_file  = open("/home/user/Desktop/rise/Antenna-Analysis/databases/Measurement_20180605/minik_20180605.txt", "r")


def find_mk_event(ant_timestamp):


	# print("Searching MiniK data...")
	# print(df)
	# df.iloc[0]["timestamp"]
	# print(tim)

	# print("... Completed.")

	for line in in_file:
		data_line = line.split("\t")#.strip("\n").strip("\t")
		
		if data_line[-1] == "\n":
			data_line.pop(-1)
			data_line += in_file.readline().split("\t")

		event_num = int(data_line[0])
		timestamp = int(data_line[1])
		
		if timestamp >= (ant_timestamp - 3) and timestamp <= (ant_timestamp + 3):
			# print(data_line[3])
			# print(data_line[4])
			azimuth = round(float(data_line[3]), 4)
			zenith = round(float(data_line[4]), 4)
			
			return event_num, timestamp, azimuth, zenith
			#minik_logger.info("    Shared Event:")
			#minik_logger.info("        Antenna:\tEvent " + str(ant_event_num) + "\tTime: "str(time) + " ns")
			#minik_logger.info("        MiniK:\tEvent " + str(event_num) + "\tTime: "str(timestamp) + " ns")

'''
if __name__  == "__main__":
	ant_event_time = 1528217308
	event_num, timestamp, azimuth, zenith = find_event(ant_event_time)
	print("MiniK Timestamp: " + str(timestamp))
	print("Antenna Timestamp: " + str(ant_event_time))
	print(round(azimuth - 194.5, 4))
	print(round(zenith, 4))
'''

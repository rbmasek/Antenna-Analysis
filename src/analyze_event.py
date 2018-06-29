import logging
import argparse
import sys
import operator
import numpy as np
import matplotlib.pyplot as plt
from numpy import array, sign, zeros
from scipy.interpolate import interp1d
from read_minik import *

# Define constants
SPEED_OF_LIGHT = 0.299792458 # In meters per nanosecond
NUMBER_OF_CHANNELS = 4

# Define the constants for the x axis
X_MIN = -2000
X_MAX = 100

# Define directories
directory = "/home/reed/Desktop/rise/Antenna-Analysis"
input_file  = directory + "/databases/Measurement_20180605/WaveDump_20180605_124246.db"

# Define all plot objects
plot_file = directory + "/analysis/radio.pdf"
fig1 = plt.figure(figsize = (12, 6))
plot1 = fig1.add_subplot(3, 1, 1)
plot2 = fig1.add_subplot(3, 1, 2)
plot3 = fig1.add_subplot(3, 1, 3)
fig2 = plt.figure(figsize = (12, 10))
plt.subplots_adjust(hspace = 1)
cmap = plt.cm.get_cmap("gist_rainbow", NUMBER_OF_CHANNELS) # Automatically assigns a color to each channel

# Define global variables
time_list = [] # Stores the time of each channel's peak frequency value. The index corresponds to the channel (i.e. time_list[0] is the time for ch0).
amplitude_list = [] # Stores the peak amplitude of each channel. The index corresponds to the channel.
event_list = [] # Stores all potential radio events across all channels
sorted_channel_list = [] # Stores the channels in the order in which they peaked
#envelope_list = [] # Stores the envelopes of each channel
#mean_list = [] # Stores the mean value of each channel. The index corresponds to the channel.

# +y axis aligned with north and +x axis aligned with east
A0 = np.array([4.80, 2.15, 0.001])
A1 = np.array([7.00, 4.82, 0.001])
A2 = np.array([0, 0, 0])

# Allows the creation of several loggers
file_formatter = logging.Formatter("%(asctime)s: %(name)s: %(levelname)-8s %(message)s")
console_formatter = logging.Formatter("%(name)s: %(levelname)-8s %(message)s")
def setup_logger(name, log_file, consol, level = logging.INFO):

    handler = logging.FileHandler(directory + "/analysis/" + log_file)        
    handler.setFormatter(file_formatter)
    
    logger = logging.getLogger(name)    
    logger.setLevel(level)
    logger.addHandler(handler)
    logging.getLogger("matplotlib").setLevel(logging.WARNING) # Suppresses matplotlib debug

    if consol == True:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)    
        console.setFormatter(console_formatter)
        logger.addHandler(console)    
    
    return logger


event_logger = setup_logger("event_logger", "events.log", consol = True)
coincidence_logger = setup_logger("coincidence_logger", "coincidences.log", consol = False)
cosmic_ray_logger = setup_logger("cosmic_ray_logger", "cosmic.log", consol = False) # Logs all events which have a coinciding signal within -1000 ns and 0 ns

#Retrieved from https://stackoverflow.com/questions/34235530/python-how-to-get-high-and-low-envelope-of-a-signal
#Creates an envelope then plots it
def create_and_plot_envelope(time, chan_num, adcValues):
    s = adcValues #This is your noisy vector of values.

    q_u = zeros(s.shape)

    #Prepend the first value of (s) to the interpolating values. This forces the model to use the same starting point for both the upper and lower envelope models.
    u_x = [0,]
    u_y = [s[0],]

    #Detect peaks and troughs and mark their location in u_x,u_y,l_x,l_y respectively.
    for k in range(1,len(s)-1):
        if (sign(s[k]-s[k-1])==1) and (sign(s[k]-s[k+1])==1):
            u_x.append(k)
            u_y.append(s[k])

    #Append the last value of (s) to the interpolating values. This forces the model to use the same ending point for both the upper and lower envelope models.
    u_x.append(len(s)-1)
    u_y.append(s[-1])

    #Fit suitable models to the data. Here I am using cubic splines, similarly to the MATLAB example given in the question.
    u_p = interp1d(u_x,u_y, kind = 'cubic',bounds_error = False, fill_value=0.0)

    #Evaluate each model over the domain of (s)
    for k in range(0,len(s)):
        q_u[k] = np.real(u_p(k))

    #Plot everything
    chan_name = "ch" + str(chan_num)
    
    col = cmap(chan_num)
    plot3.plot(time, q_u, color = col, linewidth = 3, label = chan_name + " upper envelope")

    # Find peak coordinates
    index = np.where(q_u == np.max(q_u))[0][0] # Envelope index of peak ampltiude
    ind_y = q_u[index]

    # Build string of peak coordinates and list of times
    coords = "({:6.3f}".format(time[index]) + ", {:6.3f}".format(ind_y) + ")"
    time_list.append(time[index]) # Adds the time at which the channel's signal began
    amplitude_list.append(float(ind_y))
    event_logger.info("        Envelope peak coordinate: " + coords + "")

    return q_u


def find_channel_mean(cut):
    chan_mean = np.real(np.mean(cut))
    event_logger.info("        Mean value: {:6.3f}".format(chan_mean) + " mV")

    return chan_mean


def find_signals(row, time, env_list, mean_list, bin_range, timestamp):
    
    global time_list
    if max(time_list) - min(time_list) <= 100:
        signal_found = False
        event_info = []
        event_list = []

        for chan_num, env in enumerate(env_list):
            chan_mean = np.mean(env)
            chan_max = np.max(env)
            mean_max_diff = chan_max - chan_mean

            for i in range(0, len(time) - bin_range):
                
                temp_mean = np.mean(np.real(env[i : (i + bin_range)]))
                temp_diff = temp_mean - chan_mean

                if (temp_diff > mean_max_diff * 0.5) and signal_found is False:
                    signal_found = True
                    signal_begin = time[i]

                else:
                    if ((temp_diff < mean_max_diff * 0.5) or (i == len(time) - bin_range - 1)) and signal_found is True:
                        signal_found = False
                        signal_end = time[i + bin_range]

                        if (signal_begin != time[0] and signal_end != time[-1]) and ((signal_end - signal_begin) < 300):
                            event_info.append(chan_num)
                            event_info.append(signal_begin)
                            event_info.append(signal_end)
                            event_list.append(event_info)
                            event_info = []
                            
        #Sort event_list by signal_begin
        sorted_event_list = sorted(event_list, key = operator.itemgetter(1))

        # Check for coincidence amongst the events
        event_logger.info("    Coinciding signals:")
        coincidence = False
        log_string = ""
        chans_with_c_event = [] # Will store channel numbers with coinciding event
        coincidence_list = [] # Will store information about the coinciding event

        for i, comp_event in enumerate(sorted_event_list):
            coincidence_list.append(comp_event)
            comp_event_chan = comp_event[0]
            comp_event_begin = comp_event[1]
            comp_event_end = comp_event[2]

            # Add initial event
            if len(chans_with_c_event) == 0:
                chans_with_c_event.append(comp_event_chan)
                event_begin = comp_event_begin # Marks the beginning of the event as the start of the first event

            if float(len(chans_with_c_event)) >= 0.75 * NUMBER_OF_CHANNELS:
                coincidence = True
                # Check if all channels coincide
                if len(chans_with_c_event) == NUMBER_OF_CHANNELS:
                    event_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))
                    coincidence_logger.info("    Event " + str(row) + ":")
                    coincidence_logger.info("        Event Timestamp (sec): " + str(timestamp))
                    coincidence_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    if comp_event_begin >= -1000 and comp_event_end <= 0:
                        cosmic_ray_logger.info("    Event " + str(row) + ":")
                        cosmic_ray_logger.info("        Event Timestamp (sec): " + str(timestamp))
                        cosmic_ray_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    reconstructed = find_direction(coincidence_list, timestamp)
                    if True:#reconstructed == True:
                        make_histogram(time, coincidence_list, bin_range)
                        make_heatmap(time, coincidence_list)

                    chans_with_c_event = []
                    coincidence_list = []
                    continue

                # Check if the remaining channels coincide
                if i != len(sorted_event_list) - 1:
                    next_event = sorted_event_list[i + 1]
                    next_event_chan = next_event[0]
                    next_event_begin = next_event[1]    
                
                    if next_event_chan not in chans_with_c_event and np.abs(next_event_begin - comp_event_begin) < 100:
                        chans_with_c_event.append(next_event_chan)
                        continue

                    else:
                        event_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                        coincidence_logger.info("    Event " + str(row) + ":")
                        coincidence_logger.info("        Event Timestamp (sec): " + str(timestamp))
                        coincidence_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                        if comp_event_begin >= -1000 and comp_event_end <= 0:
                            cosmic_ray_logger.info("    Event " + str(row) + ":")
                            cosmic_ray_logger.info("        Event Timestamp (sec): " + str(timestamp))
                            cosmic_ray_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                        reconstructed = find_direction(coincidence_list, timestamp)
                        if True:#reconstructed == True:
                            make_histogram(time, coincidence_list, bin_range)
                            make_heatmap(time, coincidence_list)
                        chans_with_c_event = []
                        coincidence_list = []
                        continue
                else:
                    event_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    coincidence_logger.info("    Event " + str(row) + ":")
                    coincidence_logger.info("        Event Timestamp (sec): " + str(timestamp))
                    coincidence_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    if comp_event_begin >= -1000 and comp_event_end <= 0:
                        cosmic_ray_logger.info("    Event " + str(row) + ":")
                        cosmic_ray_logger.info("        Event Timestamp (sec): " + str(timestamp))
                        cosmic_ray_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    reconstructed = find_direction(coincidence_list, timestamp)
                    if True:#reconstructed == True:
                        make_histogram(time, coincidence_list, bin_range)
                        make_heatmap(time, coincidence_list)
                    chans_with_c_event = [] 
                    coincidence_list = []
                    continue                        
                
            if i != len(sorted_event_list) - 1:
                next_event = sorted_event_list[i + 1]
                next_event_chan = next_event[0]
                next_event_begin = next_event[1]    
                
                if next_event_chan not in chans_with_c_event and np.abs(next_event_begin - comp_event_begin) < 100:
                    chans_with_c_event.append(next_event_chan)
                    continue

                else:
                    chans_with_c_event = []
                    coincidence_list = []
                    continue

        # Reset each list after use
        amplitude_list = []
        event_list = []
        time_list = []

        return coincidence
    else:
        # Reset each list after use
        amplitude_list = []
        event_list = []
        time_list = []

        return False

# Creates and prints a sorted list of the peak amplitude and time difference for each channel
def sort_channels(cut_list):
    global time_list
    global amplitude_list
    global sorted_channel_list

    difference_list = [0] * len(cut_list) # Stores the sequencial differences in TOA of the peak. The index corresponds to the channel.
    sorted_channel_list = sorted(range(len(time_list)), key = lambda k: time_list[k]) # Stores the channels in order of the signal's arrival time

    initial_time = min(time_list) # The first time at which a channel peaked
    event_logger.info("    The first envelope peaked at: " + "{:07.3f}".format(initial_time) + " ns in ch" + str(sorted_channel_list[0])) 
    event_logger.info("    The channels peaked in the following order: ")

    for channel in sorted_channel_list:
        time_difference = time_list[channel] - initial_time
        difference_list[channel] = int(time_difference)
        event_logger.info("        ch" + str(channel) + "    t: " + "{:+8.0f}".format(time_list[channel]) + " ns ({:+04.0f} ns)".format(time_difference) + " with an amplitude of " + "{:10.3f}".format(amplitude_list[channel]) + " mV")


def make_histogram(time, coincidence_list, bin_range):

    # Creates histo_list if this is the first iteration of the program
    if (time[1] - time[0] != 1):
        time = np.arange(time[0], time[-1] + 1, step = 1) # Must fix time so that it has tep of 1

    if "histo_list" not in globals():
        global histo_list
        histo_list = np.zeros((NUMBER_OF_CHANNELS, len(time)))

    if coincidence_list != None:
        for signal in coincidence_list:
            chan_num = signal[0]
            event_time = signal[1]
            event_index = np.where(time == event_time)[0][0]
            diff = event_index % bin_range
            event_index -= diff
            for b in range(0, bin_range):
                if (event_index + b) < len(time) - 1:
                    histo_list[chan_num][event_index + b] += 1
    for chan_num in range(0, NUMBER_OF_CHANNELS): 
        plot4 = fig2.add_subplot(NUMBER_OF_CHANNELS, 2, chan_num * 2 + 1)
        plot4.clear()

        col = cmap(chan_num)
        plot4.plot(time, histo_list[chan_num][:], label = "ch" + str(chan_num) + " histo", color = col)
        plot4.fill_between(time[:], histo_list[chan_num][:], y2 = 0, color = col)
        plot4.set_xlabel("Time (ns)")
        plot4.set_ylabel("Event Counts (Bin = " + str(bin_range) + " ns)")
        plot4.set_title("Channel " + str(chan_num) + " Histogram")
        plot4.set_xlim(time[0], time[len(time) - 1])
        plot4.set_ylim(0, np.max(histo_list) + 1)
        plot4.grid(1)


def make_heatmap(time, coincidence_list):

    if (time[1] - time[0] != 1):
        time = np.arange(time[0], time[-1] + 1, step = 1) # Must fix time so that it has a step of 1

    # Creates heat_list if this is the first iteration of the program
    if "heat_list" not in globals():
        global heat_list
        heat_list = np.zeros((NUMBER_OF_CHANNELS, len(time)))

    if coincidence_list != None:
        for signal in coincidence_list:
            chan_num = signal[0]
            event_time_begin = signal[1]
            event_time_end = signal[2]
            event_index = np.where(time == event_time_begin)[0][0]
            diff = int(event_time_end - event_time_begin)
            for d in range(0, diff):
                heat_list[chan_num][event_index + d] += 1

    for chan_num in range(0, NUMBER_OF_CHANNELS):
        plot5 = fig2.add_subplot(NUMBER_OF_CHANNELS, 2, chan_num * 2 + 2)
        plot5.clear()

        col = cmap(chan_num)
        heat_max = np.max(heat_list)
        heat_sum = np.sum(heat_list[chan_num][:])
        if heat_sum == 0:
            heat_sum = 1
                
        plot5.plot(time, heat_list[chan_num][:] / heat_sum, label = "ch" + str(chan_num) + " heatmap", color = col)
        plot5.fill_between(time[:], heat_list[chan_num][:] / heat_sum * 100, y2 = 0, color = col)
        plot5.set_xlabel("Time (ns)")
        plot5.set_ylabel("Relative Frequency (% of total)")
        plot5.set_title("Channel " + str(chan_num) + " Heatmap")
        plot5.set_xlim(time[0], time[len(time) - 1])
        plot5.set_ylim(0, (heat_max / heat_sum * 100) * 1.1)
        plot5.grid(1)
                 
# Method derived from arXiv:1702.04902
def find_direction(coincidence_list, timestamp):
    time_list = [] # Stores the signal begin times where the index corresponds to the channel number
    event_list = [] # Stores the event by channel number (event_list[0] corresponds to ch0)
    new_coincidence_list = []
    reconstructed = False # Returned value to determine if direction reconsruction was successful.

    # Removes the channel without matching polariation from sorted list
    for i, signal in enumerate(coincidence_list):
        chan = signal[0]
        if int(chan) != 2:
            new_coincidence_list.append(signal)

    event_list = sorted(new_coincidence_list, key = operator.itemgetter(0))

    for event in event_list:
        time_list.append(event[1])

    # Checks that all three antennas received a signal
    if len(new_coincidence_list) != 3:
        return

    # Choose the antennas that will be used in the direction reconstruction
    antenna_list = [A0, A1, A2]

    # D is a unit vector
    D = ( np.cross( (antenna_list[1] - antenna_list[0]), (antenna_list[2] - antenna_list[0]) ) ) / ( np.linalg.norm( np.cross( (antenna_list[1] - antenna_list[0]), (antenna_list[2] - antenna_list[0]) ) ) )
    D_x = D[0]
    D_y = D[1]
    D_z = D[2]

    A_tilde = np.zeros((3,3))
    xy_sqrt = np.sqrt( np.square(D_x) + np.square(D_y) )

    A_tilde[0][0] = (D_x * D_z) / xy_sqrt
    A_tilde[0][1] = (D_y * D_z) / xy_sqrt
    A_tilde[0][2] = np.negative( xy_sqrt )
    A_tilde[1][0] = np.negative(D_y) / xy_sqrt
    A_tilde[1][1] = D_x / xy_sqrt
    A_tilde[2][0] = D_x
    A_tilde[2][1] = D_y
    A_tilde[2][2] = D_z
    
    r_prime = []
    for r in antenna_list:
        r_prime.append( A_tilde.dot(r) )

    d_x_prime = SPEED_OF_LIGHT * ( ( (time_list[0] - time_list[2]) * (r_prime[1][1] - r_prime[0][1]) ) - ( (time_list[0] - time_list[1]) * (r_prime[2][1] - r_prime[0][1]) ) ) / ( ( (r_prime[2][0] - r_prime[0][0]) * (r_prime[1][1] - r_prime[0][1]) ) - ( (r_prime[1][0] - r_prime[0][0]) * (r_prime[2][1] - r_prime[0][1]) ) ) 
    d_y_prime = SPEED_OF_LIGHT * ( ( (time_list[0] - time_list[1]) * (r_prime[2][0] - r_prime[0][0]) ) - ( (time_list[0] - time_list[2]) * (r_prime[1][0] - r_prime[0][0]) ) ) / ( ( (r_prime[2][0] - r_prime[0][0]) * (r_prime[1][1] - r_prime[0][1]) ) - ( (r_prime[1][0] - r_prime[0][0]) * (r_prime[2][1] - r_prime[0][1]) ) ) 
    d_z_prime = np.sqrt( 1 - np.square(d_x_prime) - np.square(d_y_prime) )

    #d_prime is a unit vector
    d_prime = np.array((d_x_prime, d_y_prime, d_z_prime))

    if np.linalg.norm(d_prime) >= 0.999 and np.linalg.norm(d_prime) < 1.0001:
        reconstructed = True
        
        # d is a unit vector
        d = (np.linalg.inv(A_tilde)).dot(d_prime)

        zenith = np.arccos(d[2]) * 57.2958
        azimuth = np.arccos(d[0] / np.sin(zenith)) * 57.2958

        event_logger.info("        Antenna Zenith: {:5.2f}  degrees".format(zenith))
        event_logger.info("        Antenna Azimuth: {:5.2f}  degrees".format(azimuth))

        coincidence_logger.info("        Antenna Zenith: {:5.2f}  degrees".format(zenith))
        coincidence_logger.info("        Antenna Azimuth: {:5.2f}  degrees".format(azimuth))

        # Checks that the event falls within the time range in which cosmic ray events would be found
        if new_coincidence_list[0][1] >= -1000 and new_coincidence_list[0][2] <= 0:
            cosmic_ray_logger.info("        Antenna Zenith: {:5.2f}  degrees".format(zenith))
            cosmic_ray_logger.info("        Antenna Azimuth: {:5.2f}  degrees".format(azimuth))

            # MikiK data is not always structured properly. This prevents errors from being processed.
            try:
                mk_event_num, mk_timestamp, mk_azimuth, mk_zenith = find_mk_event(timestamp)

            except Exception as e:
                cosmic_ray_logger.warning("        Error when reading MiniK data: " + str(e) + ". Skipping...")
                return

            # Convert mk_azimuth so that it is with respect to North
            mk_azimuth += 194.5
            if mk_azimuth >= 360:
                mk_azimuth %= 360

            # Compares the minik angles to the antenna angles
            if True:#azimuth >= mk_azimuth - 5 and azimuth <= mk_azimuth + 5:
                event_logger.info("        MiniK Zenith: {:5.2f}  degrees".format(mk_zenith))
                event_logger.info("        MiniK Azimuth: {:5.2f}  degrees".format(mk_azimuth))
                coincidence_logger.info("        MiniK Zenith: {:5.2f}  degrees".format(mk_zenith))
                coincidence_logger.info("        MiniK Azimuth: {:5.2f}  degrees".format(mk_azimuth))
                cosmic_ray_logger.info("        MiniK Zenith: {:5.2f}  degrees".format(mk_zenith))
                cosmic_ray_logger.info("        MiniK Azimuth: {:5.2f}  degrees".format(mk_azimuth))

    return reconstructed



# Searches for coincidence amongst channels
def analyze_channels(row, time, cut_list, bin_range, timestamp):

    channel_envelopes = [] # Stores the channel envelopes. Index corresponds to channel number
    channel_means = [] # Stores the channel means. Index corresponds to channel number

    time_cut = np.round(time) # Rounds each value to an int (each value is originally a float)
    x_min_index = np.where(time_cut == X_MIN)[0][0] # Returns time index of X_MIN
    x_max_index = np.where(time_cut == X_MAX)[0][0] + 1 # Returns time index of X_MAX
    time_cut = time_cut[x_min_index : x_max_index] # Cuts the time array down to only what will be analyzed
    
    # Analyze each channel separately
    for channel_number, cut in enumerate(cut_list):
        event_logger.info("    Channel " + str(channel_number) + ":")
        channel_means.append(find_channel_mean(cut)) # Finds mean of the channel cut
        channel_envelopes.append(create_and_plot_envelope(time_cut, channel_number, cut[x_min_index : x_max_index]))

    sort_channels(cut_list)
    coincidence = find_signals(row, time_cut, channel_envelopes, channel_means, bin_range, timestamp)
    
    return coincidence

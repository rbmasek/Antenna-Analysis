import numpy as np
import logging
import operator
from numpy import array, sign, zeros
from scipy.interpolate import interp1d

# Define constants
SPEED_OF_LIGHT = 0.299792458 # In meters per nanosecond
NUMBER_OF_CHANNELS = 4

directory = "/home/user/Desktop/rise/signal_analysis"

time_list = [] # Stores the time of each channel's peak frequency value. The index corresponds to the channel (i.e. time_list[0] is the time for ch0).
amplitude_list = [] # Stores the peak amplitude of each channel. The index corresponds to the channel.
event_list = [] # Stores all potential radio events across all channels
sorted_channel_list = [] # Stores the channels in the order in which they peaked
#envelope_list = [] # Stores the envelopes of each channel
#mean_list = [] # Stores the mean value of each channel. The index corresponds to the channel.

# +y axis aligned with north and +x axis aligned with east
A0 = np.array([4.80, 2.15, 0.001])#.reshape((3,1)) # ch0
A1 = np.array([7.00, 4.82, 0.001])#.reshape((3,1)) # ch1
A2 = np.array([0, 0, 0])#.reshape((3,1)) # ch3

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


event_logger = setup_logger("event_logger", "event_log.txt", consol = True)
coincidence_logger = setup_logger("coincidence_logger", "coincidence_report.txt", consol = False)
cosmic_ray_logger = setup_logger("cosmic_ray_logger", "possible_cosmic_rays.txt", consol = False) # Logs all events which have a coinciding signal within -1000 ns and 0 ns

#Retrieved from https://stackoverflow.com/questions/34235530/python-how-to-get-high-and-low-envelope-of-a-signal
#Creates an envelope then plots it
def create_and_plot_envelope(plt, time, chan_num, adcValues):
    s = adcValues#[x_min_index : x_max_index] #This is your noisy vector of values.

    q_u = zeros(s.shape)
    #q_l = zeros(s.shape)

    #Prepend the first value of (s) to the interpolating values. This forces the model to use the same starting point for both the upper and lower envelope models.
    u_x = [0,]
    u_y = [s[0],]

    #l_x = [0,]
    #l_y = [s[0],]

    #Detect peaks and troughs and mark their location in u_x,u_y,l_x,l_y respectively.
    for k in range(1,len(s)-1):
        if (sign(s[k]-s[k-1])==1) and (sign(s[k]-s[k+1])==1):
            u_x.append(k)
            u_y.append(s[k])

        #if (sign(s[k]-s[k-1])==-1) and ((sign(s[k]-s[k+1]))==-1):
        #    l_x.append(k)
        #    l_y.append(s[k])

    #Append the last value of (s) to the interpolating values. This forces the model to use the same ending point for both the upper and lower envelope models.
    u_x.append(len(s)-1)
    u_y.append(s[-1])

    #l_x.append(len(s)-1)
    #l_y.append(s[-1])

    #Fit suitable models to the data. Here I am using cubic splines, similarly to the MATLAB example given in the question.
    u_p = interp1d(u_x,u_y, kind = 'cubic',bounds_error = False, fill_value=0.0)
    #l_p = interp1d(l_x,l_y,kind = 'cubic',bounds_error = False, fill_value=0.0)

    #Evaluate each model over the domain of (s)
    for k in range(0,len(s)):
        q_u[k] = np.real(u_p(k))
        #q_l[k] = np.real(l_p(k))

    #Plot everything
    chan_name = "ch" + str(chan_num)

    col = ""

    if chan_num == 0:
        col = "blue"
    if chan_num == 1:
        col = "red"
    if chan_num == 2:
        col = "green"
    if chan_num == 3:
        col = "orange"
    if chan_num == 4:
        col = "purple"
    if chan_num == 5:
        col = "brown"

    plt.plot(time, q_u, color = col, linewidth = 3, label = chan_name + " upper envelope")
    #plt.plot(time, q_l,'g', label = chan_name + " lower envelope")

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


def find_signals(row, fig, time, env_list, mean_list, bin_range):
    
    global time_list
    if max(time_list) - min(time_list) <= 100:
        signal_found = False
        event_info = []
        event_list = []

        for chan_num, env in enumerate(env_list):
            chan_mean = np.mean(env)
            chan_max = np.max(env)
            mean_max_diff = chan_max - chan_mean

            # env = env - chan_mean
            #env = env - np.min(env) # Sets the envelope minimum to y = 0
            #chan_max = np.max(env) # Finds the new maximum
            #chan_mean = np.mean(env)

            # print("Chan: " + str(chan_num) + " Mean: " + str(chan_mean) + " Max: " + str(chan_max)+ " Diff: " + str(mean_max_diff))

            for i in range(0, len(time) - bin_range):
                
                temp_mean = np.mean(np.real(env[i : (i + bin_range)]))# - chan_mean
                temp_diff = temp_mean - chan_mean
                # print(temp_diff)

                if (temp_diff > mean_max_diff * 0.5) and signal_found is False:
                    signal_found = True
                    signal_begin = time[i] #+ (bin_range / 2)
                    #log_string = "            A signal may begin at or around t {:.0f} ns".format(signal_begin)# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))

                else:
                    if ((temp_diff < mean_max_diff * 0.5) or (i == len(time) - bin_range - 1)) and signal_found is True:
                        signal_found = False
                        signal_end = time[i + bin_range]
                        #log_string = log_string + " and end at or around t {:.0f} ns".format(signal_end)
                        #event_logger.info(log_string)

                        if (signal_begin != time[0] and signal_end != time[-1]) and ((signal_end - signal_begin) < 300):
                            event_info.append(chan_num)
                            event_info.append(signal_begin)
                            event_info.append(signal_end)
                            event_list.append(event_info)
                            event_info = []
                            
        #Sort event_list by signal_begin
        sorted_event_list = sorted(event_list, key = operator.itemgetter(1))
        number_of_channels = len(env_list) # Should eventually make this value a global constant defined at the beginning of the execution

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

            if float(len(chans_with_c_event)) >= 0.75 * number_of_channels:
                coincidence = True
                # Check if all channels coincide
                if len(chans_with_c_event) == number_of_channels:
                    event_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))
                    coincidence_logger.info("    Event " + str(row) + ":")
                    coincidence_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    if comp_event_begin >= -1000 and comp_event_end <= 0:
                        cosmic_ray_logger.info("    Event " + str(row) + ":")
                        cosmic_ray_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    make_histogram(fig, time, coincidence_list, bin_range)
                    make_heatmap(fig, time, coincidence_list)
                    find_direction(coincidence_list)   
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
                        coincidence_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                        if comp_event_begin >= -1000 and comp_event_end <= 0:
                            cosmic_ray_logger.info("    Event " + str(row) + ":")
                            cosmic_ray_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                        print(coincidence_list) 
                        make_histogram(fig, time, coincidence_list, bin_range)
                        make_heatmap(fig, time, coincidence_list)
                        find_direction(coincidence_list)    
                        chans_with_c_event = []
                        coincidence_list = []
                        continue
                else:
                    event_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    coincidence_logger.info("    Event " + str(row) + ":")
                    coincidence_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    if comp_event_begin >= -1000 and comp_event_end <= 0:
                        cosmic_ray_logger.info("    Event " + str(row) + ":")
                        cosmic_ray_logger.info("        A coinciding signal was detected in channels " + str(chans_with_c_event) + " and begins at or around t {:.0f} ns".format(event_begin) + " and ends at or around t {:.0f} ns".format(comp_event_end))# + " ns ( of {:+07.3f}".format(cut[i + bin_range] + " mV)"))
                    make_histogram(fig, time, coincidence_list, bin_range)
                    make_heatmap(fig, time, coincidence_list)
                    find_direction(coincidence_list)    
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

        # Ensures that the histogram and heatmap is updated when a coinciding signal is found
        #for channel_number in range(0, number_of_channels):
        #    make_histogram(fig, time, None, bin_range)
        #    make_heatmap(fig, time, None)

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

    # Each list must be emptied after each event
    #time_list = []
    
    #mean_list = []
    #envelope_list = []

    #return sorted_channel_list


def make_histogram(fig, time, coincidence_list, bin_range):

    # Creates histo_list if this is the first iteration of the program
    if (time[1] - time[0] != 1):
        time = np.arange(time[0], time[-1] + 1, step = 1)# - time[0]) # Must fix time so that it has tep of 1

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
        plot4 = fig.add_subplot(NUMBER_OF_CHANNELS, 2, chan_num * 2 + 1)
        plot4.clear()
        if chan_num == 0:
            col = "blue"
        if chan_num == 1:
            col = "red"
        if chan_num == 2:
            col = "green"
        if chan_num == 3:
            col = "orange"
        if chan_num == 4:
            col = "purple"
        if chan_num == 5:
            col = "brown"

        plot4.plot(time, histo_list[chan_num][:], label = "ch" + str(chan_num) + " histo", color = col)
        plot4.fill_between(time[:], histo_list[chan_num][:], y2 = 0, color = col)
        plot4.set_xlabel("Time (ns)")
        plot4.set_ylabel("Event Counts (Bin = " + str(bin_range) + " ns)")
        plot4.set_title("Channel " + str(chan_num) + " Histogram")
        plot4.set_xlim(time[0], time[len(time) - 1])
        plot4.set_ylim(0, np.max(histo_list) + 1)
        plot4.grid(1)


def make_heatmap(fig, time, coincidence_list):

    if (time[1] - time[0] != 1):
        time = np.arange(time[0], time[-1] + 1, step = 1)# - time[0]) # Must fix time so that it has tep of 1

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
        plot5 = fig.add_subplot(NUMBER_OF_CHANNELS, 2, chan_num * 2 + 2)
        plot5.clear()
        if chan_num == 0:
            col = "blue"
        if chan_num == 1:
            col = "red"
        if chan_num == 2:
            col = "green"
        if chan_num == 3:
            col = "orange"
        if chan_num == 4:
            col = "purple"
        if chan_num == 5:
            col = "brown"

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
        plot5.set_ylim(0, (heat_max / heat_sum * 100) * 1.1)# + (heat_max * 0.2))
        plot5.grid(1)
                 

def find_direction(coincidence_list):
    # print("helo direction")
    time_list = [] # Stores the begin times in the order at which the signals began

    # Remove ch2 from sorted list
    new_coincidence_list = []
    for i, signal in enumerate(coincidence_list):
        chan = signal[0]
        time_begin = signal[1]
        if int(chan) != 2:
            new_coincidence_list.append(signal)
            time_list.append(time_begin)

    # Checks that all three antennas received a signal
    if len(new_coincidence_list) != 3:
        return

    #print(coincidence_list)
    #print(time_list)
    #print(new_coincidence_list)

    # Sort antenna coordinates by time of signal begin
    antenna_list = [A0, A1, A2]
    sorted_antenna_list = []
    for signal in new_coincidence_list:
        channel = signal[0]
        if int(channel) == 3:
            sorted_antenna_list.append(antenna_list[channel - 1])
        else:
            sorted_antenna_list.append(antenna_list[channel])

    #print(sorted_antenna_list)
    #print(sorted_antenna_list[2][1])
    #print(time_list[new_coincidence_list[0]])

    # D is a unit vector
    D = ( np.cross( (sorted_antenna_list[1] - sorted_antenna_list[0]), (sorted_antenna_list[2] - sorted_antenna_list[0]) ) ) / ( np.linalg.norm( np.cross( (sorted_antenna_list[1] - sorted_antenna_list[0]), (sorted_antenna_list[2] - sorted_antenna_list[0]) ) ) )
    D_x = D[0]
    D_y = D[1]
    D_z = D[2]
    # print(D)
    #print("D length: " + str(np.linalg.norm(D)))



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
    #print(A_tilde)
    

    r_prime = []
    for r in sorted_antenna_list:
        r_prime.append( A_tilde.dot(r) )
    print(r_prime)

    d_x_prime = SPEED_OF_LIGHT * ( ( (time_list[0] - time_list[2]) * (r_prime[1][1] - r_prime[0][1]) ) - ( (time_list[0] - time_list[1]) * (r_prime[2][1] - r_prime[0][1]) ) ) / ( ( (r_prime[2][0] - r_prime[0][0]) * (r_prime[1][1] - r_prime[0][1]) ) - ( (r_prime[1][0] - r_prime[0][0]) * (r_prime[2][1] - r_prime[0][1]) ) ) 
    d_y_prime = SPEED_OF_LIGHT * ( ( (time_list[0] - time_list[1]) * (r_prime[2][0] - r_prime[0][0]) ) - ( (time_list[0] - time_list[2]) * (r_prime[1][0] - r_prime[0][0]) ) ) / ( ( (r_prime[2][0] - r_prime[0][0]) * (r_prime[1][1] - r_prime[0][1]) ) - ( (r_prime[1][0] - r_prime[0][0]) * (r_prime[2][1] - r_prime[0][1]) ) ) 
    d_z_prime = np.sqrt( 1 - np.square(d_x_prime) - np.square(d_y_prime) )
    print(d_x_prime)
    print(d_y_prime)
    print(d_z_prime)
    #print(1 - np.square(d_x_prime) - np.square(d_y_prime))

    #d_prime is a unit vector
    d_prime = np.array((d_x_prime, d_y_prime, d_z_prime))
    if np.linalg.norm(d_prime) == 1.0:
        print("D PRIME LENGHT EQUALS ONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("d_prime length: " + str(np.linalg.norm(d_prime)))

        d = (np.linalg.inv(A_tilde)).dot(d_prime)
        print(d)

        zenith = np.arccos(d[2])
        azimuth = np.arcsin(d[1] / np.sin(zenith))
        print(zenith)
        print(azimuth)
        event_logger.info("        Zenith: {:5.2f}  degrees".format(zenith * 57.2958))
        event_logger.info("        Azimuth: {:5.2f}  degrees".format(azimuth * 57.2958))

        coincidence_logger.info("        Zenith: {:5.2f}  degrees".format(zenith * 57.2958))
        coincidence_logger.info("        Azimuth: {:5.2f}  degrees".format(azimuth * 57.2958))

        if new_coincidence_list[0][1] >= -1000 and new_coincidence_list[0][2] <= 0:
            cosmic_ray_logger.info("        Zenith: {:5.2f}  degrees".format(zenith * 57.2958))
            cosmic_ray_logger.info("        Azimuth: {:5.2f}  degrees".format(azimuth * 57.2958))

    #return zenith, azimuth



#def analyze_channels(fig, plt, time, cut_list, bin_range, X_MIN, X_MAX):
#    global envelope_list
#    global time_list 
#    global amplitude_list
#    global mean_list
#    global event_list
#    global sorted_channel_list

#    sorted_channel_list = sort_channels(cut_list)

#    return sorted_channel_list

# Plots a close up of the peak events. 
def make_peak_plot(plt, time):
    begin_times = [] # Stores all event begin times. Used to find smallest value
    end_times = [] # Stores all event end times. Used to find largest value

    for event in event_list:
        event = event.split(",")
        chan_num = int(event[0])
        signal_begin = int(event[1])
        signal_end = int(event[2])

        if time_list[chan_num] > signal_begin and time_list[chan_num] < signal_end:
            begin_times.append(signal_begin)
            end_times.append(signal_end)

            if chan_num == 0:
                col = "blue"
            if chan_num == 1:
                col = "red"
            if chan_num == 2:
                col = "green"
            if chan_num == 3:
                col = "orange"
            plt.plot(time, envelope_list[chan_num], col, linewidth = 3, label = "ch" + str(chan_num) + " peak event")

    plt.set_xlabel("Time (ns)")
    plt.set_ylabel("Amplitude (mV)")
    plt.set_title("Envelopes of Peak Events")
    plt.set_xlim(min(begin_times), max(end_times))
    plt.legend()
    plt.grid()



    
#Old stuff

# Hilbert Method
#nfft = len(ch0FourierCut)
#analytic_signal = hilbert(np.real(adcValuesCh0Cut), N = nfft)
#envelope = np.abs(np.imag(analytic_signal))


# Multitaper Method
#N = len(adcValuesCh0Cut)
#nfft = np.power( 2, int(np.ceil(np.log2(N))) )
#NW = N / 2
#fm = int(np.round(float(200) * nfft / N))
#(dpss, eigs) = nt_alg.dpss_windows(N, NW, 2 * NW)
#xk = nt_alg.tapered_spectra(adcValuesCh0Cut, dpss, nfft)
#w, n = nt_ut.adaptive_weights(xk, eigs, sides='onesided')
#mtm_bband = np.sum(2 * (xk[:, fm] * np.sqrt(eigs))[:,None] * dpss, axis = 0)# - np.abs(np.real(adcValuesCh0Cut[0]))
#mtm_bband = butter_lowpass_filter(mtm_bband, np.max(adcValuesCh0Cut), samplingFreq)


#def save_cut(t, env):
#    combined_array = np.vstack((t, env)).T
#    np.savetxt("cut_data.txt", combined_array)
#    print("Cut data saved to \"cut_data.txt\"")
#save_cut(time, adcValuesCh0Cut)

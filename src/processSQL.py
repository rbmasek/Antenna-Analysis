#!usr/bin/env python3


import sqlite3
import datetime
import matplotlib as mpl
import scipy.fftpack as sfft
from time import sleep
from analyze_event import *
# Cutting the frequencies is done in the frequencyCut.py
from frequencyCut import fourierCut


k = 1 / 12.57


def convertCh0(adcValue):
    return 0.0311 * adcValue + 2.17


def convertCh1(adcValue):
    return 0.0311 * adcValue + 2.15


def convertCh2(adcValue):
    return 0.0311 * adcValue + 2.66


def convertCh3(adcValue):
    return 0.0311 * adcValue + 2.04


def convertCh4(adcValue):
    return 0.0311 * adcValue + 2.23


def convertCh5(adcValue):
    return 0.0311 * adcValue + 2.15


def convertCh6(adcValue):
    return 0.0311 * adcValue + 2.29


def convertCh7(adcValue):
    return 0.0311 * adcValue + 2.15


'''
if len(sys.argv) == 1:
    fileLoc = input_file
else:
    fileLoc = sys.argv[1]
'''

fileLoc = input_file
conn = sqlite3.connect(fileLoc)
conn.row_factory = sqlite3.Row
print("\nConnection established\n")

c = conn.cursor()
s = conn.cursor()
ss = conn.cursor()

# Get sampling frequency
c.execute("SELECT frequency FROM digitizer")
samplingFreq = c.fetchone()[0]
samplingStepTime = 1 / samplingFreq
# print("samplingFreq: " + str(samplingFreq) \
#        + ", samplingSteps: " + str(samplingStepTime) )

c.execute("SELECT * FROM events")
rowEvents = c.fetchall()
c.execute("SELECT * FROM samples")
rowSamples = c.fetchall()
c.execute("SELECT * FROM settings_root")
rowSettings = c.fetchall()
c.execute("SELECT * FROM settings_dcoffsets")
rowDCOffsets = c.fetchall()

# Finds the number of channels in the data
c.execute("SELECT DISTINCT channel FROM samples")
NUMBER_OF_CHANNELS = len(c.fetchall())


def search_events(args):
    for row in range(args.SEARCH_MIN, args.SEARCH_MAX):

        plot1.clear()
        plot2.clear()
        plot3.clear()

        event_logger.info("################################################")
        event_logger.info("Analyzing Event " + str(row) + "...")

        event_id = rowEvents[row]['id']
        event_timestamp = rowEvents[row]['time_stamp']
        settings_id = 0
        record_length = c.execute(
            "SELECT record_length FROM settings_root").fetchone()[0]
        post_trigger = c.execute(
            "SELECT post_trigger FROM settings_root").fetchone()[0]

        dcOffsetCH0 = rowDCOffsets[0]['offset']
        dcOffsetCH0 = int((2**(14) - 1) * (dcOffsetCH0 / 65536.0))
        #event_logger.info("dcOffsetCh0: " + str(dcOffsetCH0))
        dcOffsetCH1 = rowDCOffsets[1]['offset']
        dcOffsetCH1 = int((2**(14) - 1) * (dcOffsetCH1 / 65536.0))
        #event_logger.info("dcOffsetCh1: " + str(dcOffsetCH1))
        dcOffsetCH2 = rowDCOffsets[2]['offset']
        dcOffsetCH2 = int((2**(14) - 1) * (dcOffsetCH2 / 65536.0))
        #event_logger.info("dcOffsetCh2: " + str(dcOffsetCH2))
        dcOffsetCH3 = rowDCOffsets[3]['offset']
        dcOffsetCH3 = int((2**(14) - 1) * (dcOffsetCH3 / 65536.0))
        #event_logger.info("dcOffsetCh3: " + str(dcOffsetCH3))
        dcOffsetCH4 = rowDCOffsets[4]['offset']
        dcOffsetCH4 = int((2**(14) - 1) * (dcOffsetCH4 / 65536.0))

        dcOffsetCH5 = rowDCOffsets[5]['offset']
        dcOffsetCH5 = int((2**(14) - 1) * (dcOffsetCH5 / 65536.0))

        dcOffsetCH6 = rowDCOffsets[6]['offset']
        dcOffsetCH6 = int((2**(14) - 1) * (dcOffsetCH6 / 65536.0))

        dcOffsetCH7 = rowDCOffsets[7]['offset']
        dcOffsetCH7 = int((2**(14) - 1) * (dcOffsetCH7 / 65536.0))

        event_logger.info("    Event Timestamp (sec): " + str(event_timestamp))
        # event_logger.info("Event timestamp (UTC): " + time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(event_timestamp)))

        # Shift windows_size to negative by amount of post_trigger -
        # percentage of windows_size to get time axis
        binShift = int(int(record_length) * float(100 - post_trigger) / 100)
        time = []
        for i in range(int(record_length)):
            time.append((i - binShift) * samplingStepTime * 10**9)

        # Get ADC values
        c.execute("SELECT samples FROM samples WHERE event_id = " + str(event_id))
        adcValuesCh0 = c.fetchone()[0]
        adcValuesCh1 = c.fetchone()[0]
        adcValuesCh2 = c.fetchone()[0]
        adcValuesCh3 = c.fetchone()[0]
        adcValuesCh4 = c.fetchone()[0]
        adcValuesCh5 = c.fetchone()[0]
        adcValuesCh6 = c.fetchone()[0]
        adcValuesCh7 = c.fetchone()[0]

        adcValuesCh0 = adcValuesCh0.split()
        adcValuesCh0 = [convertCh0((float(i) - dcOffsetCH0)) for i in adcValuesCh0]

        adcValuesCh1 = adcValuesCh1.split()
        adcValuesCh1 = [convertCh1((float(i) - dcOffsetCH1)) for i in adcValuesCh1]

        adcValuesCh2 = adcValuesCh2.split()
        adcValuesCh2 = [convertCh2((float(i) - dcOffsetCH2)) for i in adcValuesCh2]

        adcValuesCh3 = adcValuesCh3.split()
        adcValuesCh3 = [convertCh3((float(i) - dcOffsetCH3)) for i in adcValuesCh3]

        adcValuesCh4 = adcValuesCh4.split()
        adcValuesCh4 = [convertCh4((float(i) - dcOffsetCH4)) for i in adcValuesCh4]

        adcValuesCh5 = adcValuesCh5.split()
        adcValuesCh5 = [convertCh5((float(i) - dcOffsetCH5)) for i in adcValuesCh5]

        adcValuesCh6 = adcValuesCh6.split()
        adcValuesCh6 = [convertCh6((float(i) - dcOffsetCH6)) for i in adcValuesCh6]

        adcValuesCh7 = adcValuesCh7.split()
        adcValuesCh7 = [convertCh7((float(i) - dcOffsetCH7)) for i in adcValuesCh7]

        exp = np.log(len(adcValuesCh1)) / np.log(2)
        expDif = len(adcValuesCh1) - 2 ** int(exp)
        if expDif % 2 == 0:
            adcValuesCh0 = adcValuesCh0[int(expDif / 2):-int(expDif / 2)]
            adcValuesCh1 = adcValuesCh1[int(expDif / 2):-int(expDif / 2)]
            adcValuesCh2 = adcValuesCh2[int(expDif / 2):-int(expDif / 2)]
            adcValuesCh3 = adcValuesCh3[int(expDif / 2):-int(expDif / 2)]
            adcValuesCh4 = adcValuesCh4[int(expDif / 2):-int(expDif / 2)]
            adcValuesCh5 = adcValuesCh5[int(expDif / 2):-int(expDif / 2)]
            adcValuesCh6 = adcValuesCh6[int(expDif / 2):-int(expDif / 2)]
            adcValuesCh7 = adcValuesCh7[int(expDif / 2):-int(expDif / 2)]

            time = time[int(expDif / 2):-int(expDif / 2)]

        else:
            adcValuesCh0 = adcValuesCh0[int(expDif / 2) - 1:-int(expDif / 2)]
            adcValuesCh1 = adcValuesCh1[int(expDif / 2) - 1:-int(expDif / 2)]
            adcValuesCh2 = adcValuesCh2[int(expDif / 2) - 1:-int(expDif / 2)]
            adcValuesCh3 = adcValuesCh3[int(expDif / 2) - 1:-int(expDif / 2)]
            adcValuesCh4 = adcValuesCh4[int(expDif / 2) - 1:-int(expDif / 2)]
            adcValuesCh5 = adcValuesCh5[int(expDif / 2) - 1:-int(expDif / 2)]
            adcValuesCh6 = adcValuesCh6[int(expDif / 2) - 1:-int(expDif / 2)]
            adcValuesCh7 = adcValuesCh7[int(expDif / 2) - 1:-int(expDif / 2)]

            time = time[int(expDif / 2) - 1:-int(expDif / 2)]

        # event_logger.info(len(adcValuesCh0))
        event_logger.info("    Exp: {:6.3f}".format(exp))  # + "\n")

        ch0Fourier = np.complex()
        ch1Fourier = np.complex()
        ch2Fourier = np.complex()
        ch3Fourier = np.complex()
        ch4Fourier = np.complex()
        ch5Fourier = np.complex()
        ch6Fourier = np.complex()
        ch7Fourier = np.complex()

        ch0Fourier = sfft.fft(adcValuesCh0, len(adcValuesCh0))
        ch1Fourier = sfft.fft(adcValuesCh1, len(adcValuesCh1))
        ch2Fourier = sfft.fft(adcValuesCh2, len(adcValuesCh2))
        ch3Fourier = sfft.fft(adcValuesCh3, len(adcValuesCh3))
        ch4Fourier = sfft.fft(adcValuesCh4, len(adcValuesCh4))
        ch5Fourier = sfft.fft(adcValuesCh5, len(adcValuesCh5))
        ch6Fourier = sfft.fft(adcValuesCh6, len(adcValuesCh6))
        ch7Fourier = sfft.fft(adcValuesCh7, len(adcValuesCh7))

        freq = sfft.fftfreq(len(adcValuesCh3), samplingStepTime)
        for i in range(len(freq)):
            freq[i] = freq[i] * 10**(-6)

        ch0FourierCut = np.complex()
        ch1FourierCut = np.complex()
        ch2FourierCut = np.complex()
        ch3FourierCut = np.complex()
        ch4FourierCut = np.complex()
        ch5FourierCut = np.complex()
        ch6FourierCut = np.complex()
        ch7FourierCut = np.complex()

        ch0FourierCut = fourierCut(ch0Fourier, freq)
        ch1FourierCut = fourierCut(ch1Fourier, freq)
        ch2FourierCut = fourierCut(ch2Fourier, freq)
        ch3FourierCut = fourierCut(ch3Fourier, freq)
        ch4FourierCut = fourierCut(ch4Fourier, freq)
        ch5FourierCut = fourierCut(ch5Fourier, freq)
        ch6FourierCut = fourierCut(ch6Fourier, freq)
        ch7FourierCut = fourierCut(ch7Fourier, freq)

        adcValuesCh0Cut = sfft.ifft(ch0FourierCut, len(ch0FourierCut))
        adcValuesCh1Cut = sfft.ifft(ch1FourierCut, len(ch1FourierCut))
        adcValuesCh2Cut = sfft.ifft(ch2FourierCut, len(ch2FourierCut))
        adcValuesCh3Cut = sfft.ifft(ch3FourierCut, len(ch3FourierCut))
        adcValuesCh4Cut = sfft.ifft(ch4FourierCut, len(ch4FourierCut))
        adcValuesCh5Cut = sfft.ifft(ch5FourierCut, len(ch5FourierCut))
        adcValuesCh6Cut = sfft.ifft(ch6FourierCut, len(ch6FourierCut))
        adcValuesCh7Cut = sfft.ifft(ch7FourierCut, len(ch7FourierCut))

        # , adcValuesCh4Cut[x_min_index : x_max_index], adcValuesCh5Cut[x_min_index : x_max_index]] #Stores each channel in a list so that the commands can be iterated instead of hard coded.
        cut_list = [adcValuesCh0Cut, adcValuesCh1Cut, adcValuesCh2Cut, adcValuesCh3Cut, adcValuesCh4Cut, adcValuesCh5Cut, adcValuesCh6Cut, adcValuesCh7Cut]

        coincidence = False  # Reset boolean each iteration
        coincidence = analyze_channels(row, time, cut_list, args.BIN_RANGE, event_timestamp)

        '''
        try:
            mk_event_num, mk_timestamp, mk_a, mk_z = find_mk_event(row)
            print(event_timestamp)
            print(mk_timestamp)

            make_time_plot(event_timestamp, mk_timestamp)

        except Exception as e:
            event_logger.warning(
                "        Error when reading MiniK data: " + str(e))
		'''

        if coincidence == True:

            plot1.plot(time, adcValuesCh0, label="CH0", color=cmap(0))
            plot1.plot(time, adcValuesCh1, label="CH1", color=cmap(1))
            plot1.plot(time, adcValuesCh2, label="CH2", color=cmap(2))
            plot1.plot(time, adcValuesCh3, label="CH3", color=cmap(3))
            plot1.plot(time, adcValuesCh4, label="CH4", color=cmap(4))
            plot1.plot(time, adcValuesCh5, label="CH5", color=cmap(5))
            plot1.plot(time, adcValuesCh6, label="CH6", color=cmap(6))
            plot1.plot(time, adcValuesCh7, label="CH7", color=cmap(7))

            plot1.set_xlabel("Time (ns)")
            plot1.set_ylabel("Amplitude (mV)")
            plot1.set_title("Recorded waveforms")
            plot1.set_xlim(X_MIN, X_MAX)
            plot1.legend()
            plot1.grid(1)

            plot2.plot(freq[1:int(len(ch0Fourier) / 2)], 20 * np.log10(abs(ch0Fourier)[1:int(len(ch0Fourier) / 2)]),
                       label="Fourier Trafo CH0", lw=2, color='blue')
            plot2.plot(freq[1:int(len(ch1Fourier) / 2)], 20 * np.log10(abs(ch1Fourier)[1:int(len(ch1Fourier) / 2)]),
                       label="Fourier Trafo CH1", lw=2, color='red')
            plot2.plot(freq[1:int(len(ch2Fourier) / 2)], 20 * np.log10(abs(ch2Fourier)[1:int(len(ch2Fourier) / 2)]),
                       label="Fourier Trafo CH2", lw=2, color='green')
            plot2.plot(freq[1:int(len(ch3Fourier) / 2)], 20 * np.log10(abs(ch3Fourier)[1:int(len(ch3Fourier) / 2)]),
                       label="Fourier Trafo CH3", lw=2, color='orange')
            plot2.plot(freq[1:int(len(ch0Fourier) / 2)], 20 * np.log10(abs(ch0FourierCut)[1:int(len(ch0Fourier) / 2)]),
                       label="Fourier Trafo Cut", color='black')
            plot2.set_xlabel("Frequency (MHz)")
            plot2.set_ylabel("Power (dB)")
            plot2.set_title("Fourier transformation")
            plot2.set_ylim((0, 20 * np.log10(abs(ch0Fourier[1:]).max()) + 3))
            plot2.legend()
            plot2.grid(1)

            plot3.plot(time, adcValuesCh0Cut, label="ch0Cut", color=cmap(0))
            plot3.plot(time, adcValuesCh1Cut, label="ch1Cut", color=cmap(1))
            plot3.plot(time, adcValuesCh2Cut, label="ch2Cut", color=cmap(2))
            plot3.plot(time, adcValuesCh3Cut, label="ch3Cut", color=cmap(3))
            plot3.plot(time, adcValuesCh4Cut, label="ch4Cut", color=cmap(4))
            plot3.plot(time, adcValuesCh5Cut, label="ch5Cut", color=cmap(5))
            plot3.plot(time, adcValuesCh6Cut, label="ch6Cut", color=cmap(6))
            plot3.plot(time, adcValuesCh7Cut, label="ch7Cut", color=cmap(7))

            plot3.set_xlabel("Time (ns)")
            plot3.set_ylabel("Amplitude (mV)")
            plot3.set_title("Waveforms with frequency cuts")
            plot3.set_xlim(X_MIN, X_MAX)
            plot3.grid(1)

            plt.tight_layout()

            fig1.canvas.draw()
            fig1.show()
            fig1.canvas.flush_events()

            fig2.canvas.draw()
            fig2.show()
            fig2.canvas.flush_events()

            fig3.canvas.draw()
            fig3.show()
            fig3.canvas.flush_events()

        event_logger.info("################################################\n")


if __name__ == "__main__":

    try:
        # Define some constants
        parser = argparse.ArgumentParser(description="Define constants for the minimum event number, the maximum event number, and the bin range.")
        parser.add_argument("-min", type=int, dest="SEARCH_MIN", default=0,
                            help="integer value for the minimum event number")
        parser.add_argument("-max", type=int, dest="SEARCH_MAX", default=len(rowEvents),
                            help="integer value for the maximum event number")
        parser.add_argument("-bin", type=int, dest="BIN_RANGE", default=20,
                            help="integer value for the bin range")

        args = parser.parse_args()

        # Ensure that the specified conditions are acceptable
        if args.SEARCH_MIN > args.SEARCH_MAX or args.SEARCH_MIN == args.SEARCH_MAX:
            event_logger.error("Argument [SEARCH_MIN] (" + str(args.SEARCH_MIN) +
                               ") cannot be greater than or equal to argument [SEARCH_MAX] (" + str(args.SEARCH_MAX) + "). Exitting...")
            sys.exit()

        if args.SEARCH_MAX > len(rowEvents):
            event_logger.error("Argument [SEARCH_MAX] (" + str(args.SEARCH_MAX) +
                               ") cannot be larger than the total number of events (" + str(len(rowEvents)) + "). Exitting...")
            sys.exit()

        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info("BEGIN LOG")
        event_logger.info("Input file: " + input_file.split("/")[-1])
        event_logger.info("Number of antenna channels: " +
                          str(NUMBER_OF_CHANNELS))
        event_logger.info("Sampling Frequency: " + str(samplingFreq) + " Hz")
        event_logger.info("Sampling Steps: " + str(samplingStepTime))
        event_logger.info(
            "Event range: [" + str(args.SEARCH_MIN) + ", " + str(args.SEARCH_MAX) + ")")
        event_logger.info("Bin range: " + str(args.BIN_RANGE) + " ns")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info(
            "------------------------------------------------------------\n")

        coincidence_logger.info(
            "------------------------------------------------------------")
        coincidence_logger.info("BEGIN COINCIDENCE REPORT")
        coincidence_logger.info("Input file: " + input_file.split("/")[-1])
        coincidence_logger.info(
            "------------------------------------------------------------")

        cosmic_ray_logger.info(
            "------------------------------------------------------------")
        cosmic_ray_logger.info("BEGIN COSMIC_RAY REPORT")
        cosmic_ray_logger.info("Input file: " + input_file.split("/")[-1])
        cosmic_ray_logger.info(
            "------------------------------------------------------------")

        # Scans each database event
        search_events(args)

        # Saves the final figure to a PDF
        plt.tight_layout()
        fig1.canvas.draw()
        fig2.canvas.draw()
        fig3.canvas.draw()
        fig1.savefig(plot_file + "waveforms.pdf", bbox_inches='tight')
        fig2.savefig(plot_file + "histo.pdf", bbox_inches='tight')
        fig3.savefig(plot_file + "polar.pdf", bbox_inches='tight')

        event_logger.info("Figures saved to \"" + plot_file + "\"\n")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info("END LOG")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info(
            "------------------------------------------------------------\n\n")

        coincidence_logger.info(
            "------------------------------------------------------------")
        coincidence_logger.info("END COINCIDENCE REPORT")
        coincidence_logger.info(
            "------------------------------------------------------------\n\n")

        cosmic_ray_logger.info(
            "------------------------------------------------------------")
        cosmic_ray_logger.info("END COSMIC_RAY REPORT")
        cosmic_ray_logger.info(
            "------------------------------------------------------------\n\n")

    except KeyboardInterrupt:

        event_logger.info("Keyboard Interrupt. Exitting...\n")

        # Saves the existing figure to a PDF
        fig1.canvas.draw()
        fig2.canvas.draw()
        fig3.canvas.draw()
        fig1.savefig(plot_file + "waveforms.pdf", bbox_inches='tight')
        fig2.savefig(plot_file + "histo.pdf", bbox_inches='tight')
        fig3.savefig(plot_file + "polar.pdf", bbox_inches='tight')

        event_logger.info("Figures saved to \"" + plot_file + "\"\n")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info("END LOG")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info(
            "------------------------------------------------------------\n\n")

        coincidence_logger.info(
            "------------------------------------------------------------")
        coincidence_logger.info("END COINCIDENCE REPORT")
        coincidence_logger.info(
            "------------------------------------------------------------\n\n")

        cosmic_ray_logger.info(
            "------------------------------------------------------------")
        cosmic_ray_logger.info("END COSMIC_RAY REPORT")
        cosmic_ray_logger.info(
            "------------------------------------------------------------\n\n")

    except Exception as e:

        event_logger.info(
            "------------------------------------------------------------")
        event_logger.error("ERROR")
        event_logger.error(e)

        # Saves the existing figure to a PDF
        fig1.canvas.draw()
        fig2.canvas.draw()
        fig3.canvas.draw()
        fig1.savefig(plot_file + "waveforms.pdf", bbox_inches='tight')
        fig2.savefig(plot_file + "histo.pdf", bbox_inches='tight')
        fig3.savefig(plot_file + "polar.pdf", bbox_inches='tight')

        event_logger.info("Figures saved to \"" + plot_file + "\"\n")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info("END LOG")
        event_logger.info(
            "------------------------------------------------------------")
        event_logger.info(
            "------------------------------------------------------------\n\n")

        coincidence_logger.info(
            "------------------------------------------------------------")
        coincidence_logger.info("END COINCIDENCE REPORT")
        coincidence_logger.info(
            "------------------------------------------------------------\n\n")

        cosmic_ray_logger.info(
            "------------------------------------------------------------")
        cosmic_ray_logger.info("END COSMIC_RAY REPORT")
        cosmic_ray_logger.info(
            "------------------------------------------------------------\n\n")


conn.close()
print("\nConnection closed\n")

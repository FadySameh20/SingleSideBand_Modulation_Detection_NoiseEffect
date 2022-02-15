# SingleSideBand_Modulation_Detection_NoiseEffect
Single Side Band (SSB) modulation for both suppressed carrier and transmitted carrier by filtering the double side band (DSB) signal using ideal filter or 4th order butterworth filter. SSB detection using coherent detector in addition to showing the effect of noise on the signal with different SNR (signal to noise ratio) values.

* Procedures:
1. Use Matlab to read the attached audio file which has a sampling frequency Fs= 48 KHz. Find the spectrum of this signal.
2. Using an ideal Filter, remove all frequencies greater than 4 KHz.
3. Obtain the filtered signal in time domain, this is a band limited signal of BW=4KHz. You could play the sound back, to make sure only small distortion was introduced.
4. Generate a DSB-SC modulated signal and plot its spectrum. Choose the carrier frequency to be 100kHz. Remember to set the sampling frequency to Five times the carrier frequency Fs = 5𝐹𝑐.
5. Obtain the SSB by filtering out the USB (we need to get LSB only) of the DSB-SC modulated signal using an ideal filter then Plot the spectrum again.
6. Use coherent detection with no noise interference to get the received signal (to demodulate the SSB-SC) and play the file back also sketch the received waveform and spectrum.
7. Repeat steps 5 and 6, only this time. Use a practical 4th order Butterworth filter. [butter, filter]
8. For the ideal filter case, get the received signal again but when noise is added to SSB-SC with SNR = 0, 10, and 30 also play the received sound and sketch the received waveform/spectrum in each case. [awgn]
9. For the ideal filter case, generate a SSB-TC. As with experiment one, set the DC bias to twice the maximum of the message. Use envelope detector to demodulate the message (without noise). Play back the received message and sketch its waveform.

import pyaudio
import numpy as np


def frequency(eingenvalue, speed_sound = 300): 

    frequency  = speed_sound*800*(eingenvalue)**(1/2)/(np.pi*2)

    return frequency


def play_sound(volume, sampling_rate, duration, freq) : 

    p = pyaudio.PyAudio()

    samples = (np.sin(2*np.pi*np.arange(sampling_rate*duration)*freq/sampling_rate)).astype(np.float32)

    stream = p.open(format = pyaudio.paFloat32, channels = 1, rate = sampling_rate, output = True)

    stream.write(volume*samples)

    stream.stop_stream()
    stream.close()




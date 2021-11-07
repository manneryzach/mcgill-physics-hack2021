import pyaudio
import numpy as np
from pygame import mixer


def frequency(eingenvalue, speed_sound = 300): 

    frequency  = speed_sound*800*(eingenvalue)**(1/2)/(np.pi*2)

    return frequency


def play_sound(volume, sampling_rate, duration, freq) : 

    if (2000 < freq < 2100) : 

        file = 'therock.mp3'
        mixer.init()
        mixer.music.load(file)
        mixer.music.set_volume(0.6)
        mixer.music.play()
    
    
    else : 
        p = pyaudio.PyAudio()

        samples = (np.sin(2*np.pi*np.arange(sampling_rate*duration)*freq/sampling_rate)).astype(np.float32)

        stream = p.open(format = pyaudio.paFloat32, channels = 1, rate = sampling_rate, output = True)

        stream.write(volume*samples)

        stream.stop_stream()
        stream.close()

play_sound(0.5 ,44100, 1, 100 )



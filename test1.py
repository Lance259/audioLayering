#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.io.wavfile import write, read
from scipy.signal import stft, istft


class waveLayer:
    
    def __init__(self, samplerate=44100, length_s=10):
        self.length = length_s
        self.samplerate = samplerate
        self.audio_layers = []
        self.numLayers = 0
        self.blend_modes = []
        self.names = {}
        self.out_render = np.zeros(samplerate*length_s) # can be resized later
    
    """ 
    contains: 
        addLayer(path_to_file, name/index, blend_mode)
        
        render() # applies all blend modes and compiles all layers into one array
        play() # plays the current layer? 
        export() # exports the current layers with blend modes to wav
        
    """
    
    def add_layer(self, source, name='', gain=1.0, blend_mode=None):
        #todo: maybe check for exceptions
        if(isinstance(source, str)):
            try:
                rate, data = read(source)
                if(rate != self.samplerate):
                    print('WARNING: sample rate of imported file does not match global samplerate')
            except:
                print('problem reading wav file ' + str(source))
                return -1
        elif(isinstance(source, np.ndarray)):
            data = source
        else:
            print('invalid input')
            return -2
        self.audio_layers.append(np.array(data)*gain)
        self.blend_modes.append(blend_mode)
        self.names[name] = self.numLayers
        self.numLayers += 1
        return 0
        
    def export(self, path_to_file='export.wav'):
        # scale to 16 bit wave resolution, export
        scaled_16Bit = np.int16(self.out_render  / np.max(np.abs(self.out_render)) * 32767)
        #return scaled_16Bit
        write(path_to_file, self.samplerate, scaled_16Bit)

    def render(self):
        # goes through all layers one by one, and applies blend modes
        # todo: choose specific layers
        if len(self.audio_layers) == 1:
            # nothing to render if only one layer
            self.out_render = np.array(self.audio_layers[0])
            
        else: 
            intermediate_render = np.array(self.audio_layers[0])
            for i in range (1, len(self.audio_layers)):
                try:
                    intermediate_render = self.blend_modes[i](lower_layer = intermediate_render, higher_layer = self.audio_layers[i])
                except:    
                    print('no valid blend mode for layer ' + str(self.names[i]))
            self.out_render = intermediate_render
        return self.out_render
     
    def shift_layer(self, name, sample_shift=0):
        sample_shift = int(sample_shift)
        if(sample_shift > 0):
            self.audio_layers[self.names[name]] = np.pad(self.audio_layers[self.names[name]],((sample_shift,0),(0,0)), mode='constant')
        elif(sample_shift < 0):
            self.audio_layers[self.names[name]]  = self.audio_layers[self.names[name]][-sample_shift:,:]
    
    def get_layer(self, name):
        return self.audio_layers[self.names[name]]
    
    #def pan(name, pan) #pan between -1 (left) and 1 (right)
    
    

###################### BLEND MODES #######################


#defining blend modes
def normalize(arr):
    return arr / np.max(arr)

def add(lower_layer, higher_layer):
    # find the max array length and pad the arrays to same length
    length = np.max([lower_layer.shape[0], higher_layer.shape[0]])
    lower_layer = np.pad(lower_layer, ((0,length - lower_layer.shape[0]),(0,0)), 'constant')
    higher_layer = np.pad(higher_layer, ((0,length - higher_layer.shape[0]),(0,0)), 'constant')
    return lower_layer + higher_layer

def mul(lower_layer, higher_layer):
    # find the max array length and pad the arrays to same length
    length = np.max([lower_layer.shape[0], higher_layer.shape[0]])
    lower_layer = np.pad(lower_layer, ((0,length - lower_layer.shape[0]),(0,0)), 'constant')
    higher_layer = np.pad(higher_layer, ((0,length - higher_layer.shape[0]),(0,0)), 'constant')
    return higher_layer * lower_layer

def crosssynth(lower_layer, higher_layer): #lower_layer = carrier, higher_layer = modulator
    # compute stft of lower_layer (carrier)
    # and higher_layer (modulator)
    
    # spectral envelope? 
    # -> find peaks in fft, polynomial interpolation between those peaks
    
    # flatten spectrum of carrier (divide by own spectral envelope)?
    
    # multiply both spectra? 
    
    # apply inverse stft
    
    n_per_seg = 512
    lower_layer = np.array(waveLayer2.get_layer('cello_mid'))
    higher_layer = np.array(waveLayer2.get_layer('drums'))
    
    length = np.max([lower_layer.shape[0], higher_layer.shape[0]])
    lower_layer = np.pad(lower_layer, ((0,length - lower_layer.shape[0]),(0,0)), 'constant')
    higher_layer = np.pad(higher_layer, ((0,length - higher_layer.shape[0]),(0,0)), 'constant')
    
    carrier = lower_layer.T
    modulator = higher_layer.T
    
    carrier[0] = normalize(carrier[0])
    carrier[1] = normalize(carrier[1])
    Carrier_f, Carrier_t, Carrier_Zxx = stft(carrier, samplerate, window='hann', nperseg=n_per_seg)
    
    modulator[0] = normalize(modulator[0])
    modulator[1] = normalize(modulator[1])
    Modulator_f, Modulator_t, Modulator_Zxx = stft(modulator, samplerate, window='hann', nperseg=n_per_seg)
    
    Modulated_stft = np.zeros((Carrier_Zxx.shape), dtype=np.complex256)
    for channel in range(2):
        for timeslice in range(Carrier_Zxx.shape[1]):        
            Modulated_stft[channel][timeslice] = Carrier_Zxx[channel][timeslice] * Modulator_Zxx[channel][timeslice]
        
    plt.pcolormesh(Carrier_t, Carrier_f, np.abs(Carrier_Zxx[0]), vmin=0, vmax=0.01, shading='gouraud')
    
    
    t, x = istft(Modulated_stft, fs=samplerate, window='hann')
    
    return x.T

    

#################### MAIN PROGRAM ####################



    
audio_files = ['/home/lhinz/Documents/python_tools/audioLayering/audio/0_Drum1.wav',
               '/home/lhinz/Documents/python_tools/audioLayering/audio/Pali_nochNieErlebt.wav',
               '/home/lhinz/Documents/python_tools/audioLayering/audio/Pali_stickOderBag.wav']    
samplerate = 44100

       
waveLayer1 = waveLayer(samplerate=samplerate) 
tempWaveLayer = waveLayer(samplerate=samplerate) 
tempWaveLayer.add_layer(audio_files[2], 'Pali_stickOderBag2', add) 
tempWaveLayer.add_layer(np.ones((5*len(tempWaveLayer.get_layer('Pali_stickOderBag2')),2)), 'const1', add) 
tempWaveLayer.shift_layer('Pali_stickOderBag2', 2*len(tempWaveLayer.get_layer('Pali_stickOderBag2')))


waveLayer1.add_layer(audio_files[2], 'Pali_stickOderBag', add) 
waveLayer1.add_layer(audio_files[1], 'Pali_nochNieErlebt', add)
waveLayer1.add_layer(audio_files[0], 'drum', add)
waveLayer1.add_layer(tempWaveLayer.render(), 'test', mul)




waveLayer1.shift_layer('Pali_stickOderBag', 1*samplerate)
#waveLayer1.shift_layer('test', 3*samplerate)
print(waveLayer1.get_layer('Pali_stickOderBag').shape)
print(waveLayer1.get_layer('Pali_stickOderBag'))

waveLayer1.render()
waveLayer1.export()

#%%

# cello mit drums multiplizieren? 
audio_files = ['/home/lhinz/Documents/python_tools/audioLayering/audio/0_Drum1.wav',
               '/home/lhinz/Documents/python_tools/audioLayering/audio/Pali_nochNieErlebt.wav',
               '/home/lhinz/Documents/python_tools/audioLayering/audio/Pali_stickOderBag.wav',
               '/home/lhinz/Documents/python_tools/audioLayering/audio/01_C3_bowed_long_vibrato.wav',
               '/home/lhinz/Music/SuperCollider/samples/samples/samples-cello-bowed/07_C5_bowed_long.wav']    
samplerate = 44100

waveLayer2 = waveLayer(samplerate=samplerate)

waveLayer2.add_layer(audio_files[1], 'cello_mid', 0.5, add)
waveLayer2.add_layer(audio_files[4], 'cello_high', 0.5, add)

waveLayer2.shift_layer('cello_high', len(waveLayer2.get_layer('cello_mid'))-2*samplerate)
waveLayer2.add_layer(audio_files[0], 'drums', 1.0, mul)
waveLayer2.render()
waveLayer2.export('export2.wav')

#%%


# cello mit drums multiplizieren? 
audio_files = ['/home/lhinz/Documents/python_tools/audioLayering/audio/0_Drum1.wav',
               '/home/lhinz/Documents/python_tools/audioLayering/audio/Pali_nochNieErlebt.wav',
               '/home/lhinz/Documents/python_tools/audioLayering/audio/Pali_stickOderBag.wav',
               '/home/lhinz/Documents/python_tools/audioLayering/audio/01_C3_bowed_long_vibrato.wav',
               '/home/lhinz/Music/SuperCollider/samples/samples/samples-cello-bowed/07_C5_bowed_long.wav']    
samplerate = 44100


waveLayer3 = waveLayer(samplerate=samplerate)
waveLayer3.add_layer(audio_files[0], 'drum', 0.5, add)
waveLayer3.add_layer(audio_files[1], 'nieErlebt', 0.5, crosssynth)
waveLayer3.render()
waveLayer3.export('export3.wav')
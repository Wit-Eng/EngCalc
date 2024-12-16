# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:18:15 2024

@author: pimwi
"""

import numpy as np
from matplotlib import pyplot as plt
from math import pi, sqrt, log
import scipy.stats
import streamlit as st

st.title("Engineering Calculator")

# description = ''' This is a demo front end for the Spiral simulator

# '''
# st.markdown(description)

tab1, tab2 = st.tabs(["cantilever dynamics","cantilever bending"])

with tab1:
    st.header("Cantilever dynamics")
    col1, col2, col3 = st.columns([2,1,4])
    col1.write("Input")
    col3.write("Output")
    
    with col1:
        # define units
        # m = 1
        # mm = 1E-3
        # um = 1E-6
        # nm = 1E-9
        # ms = 1E-3
        # Hz = 1
        # dB = 1
        # kg = 1
        # sec = 1
        # pi = 3.141592653589793
        
        # Frequency range
        f_low = 1                                  # Lower limit of input frequency range [Hz]
        f_high = 2000                                 # Upper limit of input frequency range [Hz]
        df = 0.01                                   # Increment of input frequency range [Hz]
        f = np.arange(f_low, f_high, df)                # Input frequencies [Hz]
        w = f*2*pi                                      # Input frequencies [rad/s]
        s = complex('j')*w                              # Laplace domain variable s = j*omega
    
        # System
        L = st.number_input("Cantilever length [mm]: ",None,None,100)/1000
        w = st.number_input("Cantilever width [mm]: ",None,None,10)/1000
        t = st.number_input("Cantilever thickness [mm]: ",None,None,1)/1000
        I = (1/12)*w*t**3
        V = L*w*t
        rho = 7850
        M = rho*V
        Meff = (33/140)*M
        E = 2e11
        k = (3*E*I)/(L**3)
        wn = sqrt(k/Meff)/(2*pi)
    
        # print(wn,wn*2*pi)
    
        # zeta = 0.05                                    # stiffness
        zeta = st.number_input("damping ratio: []", None,None,0.05)
        z_P = 2*(zeta)*sqrt(k*Meff)                     # Damping ratio [-]
    
        G = 1/(Meff*s**2 + z_P*s + k)                      #response to external input
        H_P = (z_P*s + k)/(1*s**2 + z_P*s + k)
        # H_P = (z_P*s + k)/(Meff*s**2 + z_P*s + k)
        # Plot results
        xlim = (f_low, f_high)
    
    with col3:
        st.write("Stiffness [N/m]: " , k)
        st.write("resonance [Hz]: ", wn)
        
        fig, ax = plt.subplots()
        plt.loglog(f,(np.abs(G)))
        plt.grid(which='both')
        ax.set_title('System frequency response')
        ax.set_ylabel('Magnitude [m/N]')
        ax.set_xlabel("Frequency [Hz]")
        st.pyplot(fig)
        
        fig, ax = plt.subplots()
        plt.semilogx(f, 20*np.log10(np.abs(G)), label='System', color='red')
        plt.grid(which='both')
        ax.set_title('System frequency response')
        ax.set_ylabel('Magnitude [dB]')
        ax.set_xlabel("Frequency [Hz]")
        st.pyplot(fig)
        
        fig, ax = plt.subplots()
        plt.semilogx(f, (180/pi)*np.unwrap(np.angle(G)), label='System', color='red')
        ax.set_xlabel('Frequency [Hz]', fontsize=12)
        ax.set_ylabel(r'Phase $\left[\degree\right]$', fontsize=12)
        ax.set_xlim(xlim)
        # ax.set_yticks(np.arange(-180, 360+90, 90))
        ax.grid(which='both')
        st.pyplot(fig)
        
        # ax, fig = plt.subplots()
        # # ax1 = plt.subplot(211)
        # ax.semilogx(f, 20*np.log10(np.abs(G)), label='System', color='red')
        # # ax1.semilogx(f, 20*np.log10(np.abs(G_CL)), label='Closed loop', color='blue')
        # ax.set_title('System frequency response', fontsize=12)
        # ax.set_ylabel('Magnitude [dB]', fontsize=12)
        # ax.set_xlim(xlim)
        # # ax.set_yticks(np.arange(-40, 180+20, 20))
        # ax.grid(which='both')
        # ax.legend()
    
        # ax, fig = plt.subplots()
        # ax = plt.subplot(212)
        # ax.semilogx(f, (180/pi)*np.unwrap(np.angle(G)), label='System', color='red')
        # ax.set_xlabel('Frequency [Hz]', fontsize=12)
        # ax.set_ylabel(r'Phase $\left[\degree\right]$', fontsize=12)
        # ax.set_xlim(xlim)
        # # ax.set_yticks(np.arange(-180, 360+90, 90))
        # ax.grid(which='both')
        
        
        
        
    
with tab2:
    st.header("Cantilever dynamics")
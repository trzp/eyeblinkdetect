#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/11/15 22:24
# @Version : 1.0
# @File    : eyeblink.py
# @Author  : Jingsheng Tang
# @Version : 1.0
# @Contact : mrtang@nudt.edu.cn   mrtang_cs@163.com
# @License : (C) All Rights Reserved

from __future__ import division
from __future__ import print_function
from scipy import signal as scipy_signal
import numpy as np
#from matplotlib import pyplot as plt


# 基本原理
# 1-5hz滤波
# 降采样
# 峰值检测（波峰波谷）
# 对波峰-波谷（可表征一次眨眼）的幅值差和宽度进行判断，没有使用绝对的幅值，可靠性更高

class EyeBlink():
    def __init__(self,srate=250,notch=False,th_w = [60,250],th_h = 150):
        self.srate = 250
        self.notch = notch

        # notch filter
        fs = self.srate / 2.
        Fo = 50
        Q = 15
        w0 = Fo / (fs)
        self.notchB, self.notchA = scipy_signal.iirnotch(w0=w0, Q=Q)

        # bandpass filter
        Wp = np.array([1 / fs, 5 / fs])
        Ws = np.array([0.5 / fs, 10 / fs])
        N, Wn = scipy_signal.cheb1ord(Wp, Ws, 3, 40)
        self.bpB, self.bpA = scipy_signal.cheby1(N, 0.5, Wn, 'bandpass')

        # 均值降采样
        targetSrate = 50
        self.meanNum = int(self.srate / targetSrate) # 将信号降采样到50Hz附近
        self.meanNumf = np.ones(self.meanNum)/self.meanNum

        tu = self.meanNum * 1000 / self.srate # 降采样后采样点间隔时间
        self.th_w = [int(i/tu) for i in th_w]
        self.th_h = th_h

        self.signal = np.array([0])

    def detect(self,sig):
        if self.notch:
            sig = scipy_signal.filtfilt(self.notchB, self.notchA, sig)
        sig = scipy_signal.filtfilt(self.bpB, self.bpA, sig)
        sig = scipy_signal.lfilter(self.meanNumf,1,sig)
        sig = sig[1+int(0.5*self.meanNum)::self.meanNum]    #均值降采样
        locs = scipy_signal.find_peaks(sig)[0]
        pks = sig[locs]
        locs1 = scipy_signal.find_peaks(-1*sig)[0]
        pks1 = sig[locs1]
        tem = np.vstack([locs,pks,np.ones(locs.size)])
        tem1 = np.vstack([locs1,pks1,-1*np.ones(locs1.size)])
        ppks = np.hstack([tem,tem1])
        ppks = ppks[:,np.argsort(ppks[0,:])]
        if ppks[2,0]==-1:
            ppks = ppks[:,1:]
        dppks = np.diff(ppks,axis=-1)
        ind = np.where((dppks[0,:]>self.th_w[0]) & (dppks[0,:]<self.th_w[1]) & (dppks[1,:]<-self.th_h))[0]
        blink = ppks[:2,ind]
        return blink

if __name__ == '__main__':
    eog = np.load('eog.npy')
    signal = eog[:,0].transpose()
    eye = EyeBlink(srate=250,notch=True)
    eye.detect(signal)






import os, time, sys, threading, math
from os import listdir
from os.path import abspath, isfile, join
from pathlib import Path
import datetime
import csv
import copy
import warnings

import numpy as np
from astropy.io import fits
from astropy.io.fits.verify import VerifyError
from astropy import wcs
import astropy.units as u
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.modeling import models, fitting
import PIL.Image as PILimage

from PyQt5.QtWidgets import QDesktopWidget

from ginga import Bindings
from ginga.misc import log
from ginga.qtw.QtHelp import QtGui, QtCore
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.util import iqcalc
from ginga.util.loader import load_data
from ginga.AstroImage import AstroImage

import ktl

class ScannerSignals(QtCore.QObject):
    load = QtCore.Signal(object)

class Scanner(QtCore.QRunnable):
    '''
    Scanner thread
    '''
    def __init__(self, fn, *args, **kwargs):
        super(Scanner, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = ScannerSignals()
        self.kwargs['file_callback'] = self.signals.load

        # Add the callback to our kwargs
    @QtCore.Slot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        self.fn(*self.args, **self.kwargs)

class VideoSignals(QtCore.QObject):
    load = QtCore.Signal(object, object, object)

class Video(QtCore.QRunnable):
    '''
    Video thread
    '''
    def __init__(self, fn, *args, **kwargs):
        super(Video, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = VideoSignals()
        self.kwargs['file_callback'] = self.signals.load

        # Add the callback to our kwargs
    @QtCore.Slot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        self.fn(*self.args, **self.kwargs)

class UpdateControlWindowSignals(QtCore.QObject):
    load = QtCore.Signal()

class UpdateControlWindow(QtCore.QRunnable):
    '''
    Control Window thread
    '''
    def __init__(self, fn, *args, **kwargs):
        super(UpdateControlWindow, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = UpdateControlWindowSignals()
        self.kwargs['file_callback'] = self.signals.load

        # Add the callback to our kwargs
    @QtCore.Slot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        self.fn(*self.args, **self.kwargs)



class Util:
    def __init__(self):
        super().__init__()

        #KTL keywords
        self.trickxpos = ktl.cache('tds', 'TRKRO1X')
        self.trickypos = ktl.cache('tds', 'TRKRO1Y')
        self.trickxsize = ktl.cache('tds', 'TRKRO1SX')
        self.trickysize = ktl.cache('tds', 'TRKRO1SY')
        self.stopex = ktl.cache('tds', 'STOPEX')
        self.timfile = ktl.cache('tds', 'TIMFILE')
        self.init = ktl.cache('tds', 'INIT')
        self.sampmode = ktl.cache('tds', 'SAMPMODE')
        self.cdsmode = ktl.cache('tds', 'CDSMODE')
        self.readmode = ktl.cache('tds', 'READMODE')
        self.itime = ktl.cache('tds', 'ITIME')
        self.getkw = ktl.cache('tds', 'GETKW')
        self.getdcskw = ktl.cache('tds', 'GETDCSKW')
        self.getaokw = ktl.cache('tds', 'GETAOKW')
        self.go = ktl.cache('tds', 'GO')
        self.trkro1px = ktl.cache('ao', 'TRKRO1PX')
        self.trkro1px.monitor()
        self.trkro1ff = ktl.cache('ao', 'TRKRO1FF')
        self.trkro1ff.monitor()
        self.trkro1bg = ktl.cache('ao', 'TRKRO1BG')
        self.trkro1bg.monitor()

        self.ops = "MGAO" #set mode to MGAO by default, then checks for missing keywords
        try:
            self.trkenapx = ktl.cache('ao', 'trkenapx')
            self.trkfpspx = ktl.cache('ao', 'trkfpspx')
        except KeyError:
            self.ops = "RTC"

        self.tkenrup = ktl.cache('ao', 'tkenrup')
        self.tkcrxs = ktl.cache('ao','tkcrxs')
        self.tkcrys = ktl.cache('ao','tkcrys')
        self.tkcrevxp = ktl.cache('ao','tkcrevxp')
        self.tkcrevyp = ktl.cache('ao','tkcrevyp')
        self.tkcrevxo = ktl.cache('ao','tkcrevxo')
        self.tkcrevyo = ktl.cache('ao','tkcrevyo')
        self.tksrtrg = ktl.cache('ao','tksrtrg')
        self.trkrordy = ktl.cache('ao','trkrordy')
        self.tkcxim = ktl.cache('ao','tkcxim')
        self.tkcyim = ktl.cache('ao','tkcyim')
        self.trkstop = ktl.cache('trick', 'trkstop')
        self.trkstsx = ktl.cache('trick', 'trkstsx')
        self.trkrocpr = ktl.cache('trick', 'trkrocpr')
        self.cyclespr = ktl.cache('tds', 'cyclespr') #TODO figure out differnce in cpr
        self.trknmad1 = ktl.cache('ao', 'trknmad1')
        self.targname = ktl.cache('tfs', 'TARGNAME')
        self.targname.monitor()
        self.tfshome = ktl.cache('tfs', 'home')
        self.tfsinit = ktl.cache('tfs', 'init')
        self.tfsstatus = ktl.cache('tfs', 'status')
        self.tfsposname = ktl.cache('tfs', 'posname')
        self.trkstat = ktl.cache('trick', 'trkstat')
        self.trkstat.monitor()
        self.dtlp = ktl.cache('ao', 'dtlp')
        self.dttmastr = ktl.cache('ao', 'dttmastr')
        self.dtsensor = ktl.cache('ao', 'dtsensor')
        self.dtsensor.monitor()

    def trk_set_cpr(self, cpr):
        cpr_now = int(self.trkrocpr.read())
        if cpr == cpr_now:
            print('Requested CPR already in place')
            return

        print('Setting CPR to %d' % cpr)
        trkstat_init = self.trkstat.read()
        if trkstat_init == '2nd-channel video in progress':
            # if the tip-tilt loop is closed on TRICK, open it briefly then reclose it
            dtlp = self.dtlp.read()
            dttmastr = self.dttmastr.read()
            if (dtlp == 'CLOSE') and (dttmastr == 'ROI1'):
                print('Opening the tip-tilt loop')
                ktl.cache('ao', 'dtlp').write('OPEN')

            print('Stopping exposure while we change CPR')
            self.trkstop.write('1')

        status = self.trkrocpr.write(cpr)

        if trkstat_init == '2nd-channel video in progress':
            print('Restarting video')
            status = self.trkstsx.write(1) # restart the video
            time.sleep(0.02)
            dtlp = self.dtlp.write(dtlp)      # return the loop to original state

    def trk_set_coadds(self, ncoadds):

        # error checking on the number of coadds; must be an integer >= 1
        ncoadds = int(ncoadds)
        ncoadds = np.clip(1,ncoadds,100) # not sure what the upper limit is!

        # check that this is not already the value in place; if so, don't set
        ncoadds_now = int(self.trknmad1.read())
        if ncoadds == ncoadds_now:
            print('Requested number of coadds already in place')
            return

        # check that the CPR is an integer number of ncoadds
        cpr_now = int(self.trkrocpr.read())
        if cpr_now % ncoadds != 0:
            print('CPR is not an integer multiple of requested ncoadds')
            # we will have to make it an integer number of ncoadds for *both* the current and new number of coadds
            # calculate the lowest common multiple of the current and new number of coadds
            lcm = np.lcm(ncoadds_now,ncoadds)
            cpr_new = int(np.ceil(cpr_now/lcm)*lcm)
            trk_set_cpr(cpr_new)

        status = self.trknmad1.write(ncoadds)

    def desaturate(self):
        print("Desaturate button pressed")
        cpr = int(self.trkrocpr.read())
        coadd = int(self.trknmad1.read())
        xsatu = int(128*cpr/(69700*coadd^2)*100) #info.maxroi = 128
        if xsatu >= 40:
            newcpr = int(cpr*40/xsatu)
            if newcpr < 2: #cannot be less than 2
                newcpr = 2

            if cpr != newcpr:
                print(f'Resetting CPR to {newcpr}')

                status = self.trkstop.write(1)
                status = self.trkrocpr.write(newcpr)
                status = self.trkstsx.write(1)
            else:
               print('Unable to change CPR to desaturate detector')
        else:
            print('The detector is not near saturation')

    def trk_distortion_model(self, x0, y0, p):
        # Translation of DistortionModel into Python
        x0T = x0-p[0]
        y0T = y0-p[10]

        x1 = p[0]+p[1]*y0T+p[2]*x0T+p[3]*y0T**2+p[4]*y0T*x0T+p[5]*x0T**2+\
            p[6]*y0T**3+p[7]*y0T**2*x0T+p[8]*y0T*x0T**2+p[9]*x0T**3

        y1 = p[10]+p[11]*y0T+p[12]*x0T+p[13]*y0T**2+p[14]*y0T*x0T+p[15]*x0T**2+\
            p[16]*y0T**3+p[17]*y0T**2*x0T+p[18]*y0T*x0T**2+p[19]*x0T**3

        return (x1,y1)

    def trk_undo_distmodel(self, x1, y1, DistC):
        eps = 0.5 # error in pixels
        maxiter = 20

        k = 0
        x0 = copy.copy(x1)
        y0 = copy.copy(y1)
        (x1p,y1p) = self.trk_distortion_model(x0,y0,DistC)

        while (np.sqrt((x1p-x1)**2+(y1p-y1)**2) > eps):
            if k >= maxiter:
                print('No convergence in '+str(iter)+' iterations')
                break

            x0 -= (x1p-x1)
            y0 -= (y1p-y1)

            (x1p,y1p) = self.trk_distortion_model(x0,y0,DistC)
            k += 1

        return (x0,y0)


    def trk_putxy_spoc(self, xroi, yroi, distcoeff, roisz=None):
        # TODO: error handling

        if roisz not in [None,2, 4, 8, 16]:
            raise('Invalid value for roisz')
            return

        sctrack = self.tkenrup.read() # what is the error handling here?
        if sctrack:
            # turn sctracking off
            status = self.tkenrup.write(0)

        if roisz is not None:
            status = self.tkcrxs.write(roisz)
            status = self.tkcrys.write(roisz)

        (xp,yp) = self.trk_undo_distmodel(xroi, yroi, distcoeff)
        status = self.tkcrevxp.write(xp)
        status = self.tkcrevyp.write(yp)

        xim = self.tkcrevxo.read()
        yim = self.tkcrevyo.read()
        status = self.tkcxim.write(xim)
        status = self.tkcyim.write(yim)

        # need to tigger the ROI calculator TWICE (as told to Bruno Femenia by Paul Stomski)
        status = self.tksrtrg.write(1)
        time.sleep(0.1)    # wait 0.1 s
        status = self.tksrtrg.write(1)

        # need to set to propagate SPOC-Camera-WFC
        status = self.tkenrup.write(1)
        if sctrack == 0:
            status = self.tkenrup.write(0)

        # Setting trkrordy syncs SPOC and camera values with WFC and triggers actions by WFC: bad pixel map, etc
        status = self.trkrordy.write(1)
        return

    def change_filter(self, button):
        target = button.text()

        if target != str(self.targname):
        #     if self.tiptilt == "closed": #todo figure out keyword here
        #         self.tiptilt_popup()
        #         return
            if target == 'DISMISS':
                return
            if target == 'Block':
                target = 'blocking'
            if target == 'home' or target == 'open':
                target = target.lower()
            self.targname.write(target)
        return

    def init_filter(self):
        print("Initing filter")
        #if self.tiptilt == "closed": #todo figure out keyword here
        #         self.tiptilt_popup()
        #         return
        targname = str(self.targname)
        print('Sending TRICK filter wheel home')
        status = self.tfshome.write(1)
        time.sleep(3)
        print('Initializing TRICK filter wheel')
        status = self.tfsinit.write(1)
        tfsstatus = self.tfsstatus.read()
        while tfsstatus != 'OK':
            print('tfsstatus')
            time.sleep(1)
            tfsstatus = self.tfsstatus.read()
        posname = self.tfsposname.read()
        if posname != 'home':
            print('Waiting for TFS to go home, this could take 30s')
        while posname != 'home':
            time.sleep(0.5)
            posname = self.tfsposname.read()
        print(f'Sending TRICK filter wheel to {targname}')
        status = self.targname.write(targname)



class ControlWindow(QtGui.QWidget):
    """
    This "window" is a QWidget. If it has no parent, it
    will appear as a free-floating window as we want.
    """
    def __init__(self):
        super().__init__()

        self.util = Util()

        vbox = QtGui.QVBoxLayout()
        title_hbox = QtGui.QHBoxLayout()
        title = QtGui.QLabel("Trick Control")
        title.setAlignment(QtCore.Qt.AlignHCenter)
        title_hbox.addWidget(title)
        vbox.addLayout(title_hbox)
        roisz_hbox = QtGui.QHBoxLayout()
        self.roisz_label = QtGui.QLabel("ROI Size: ")
        roisz_hbox.addWidget(self.roisz_label)
        self.wroisz = QtGui.QComboBox()
        sizes = ['2','4','8','16']
        for size in sizes:
            self.wroisz.addItem(size)
        item = str(self.util.trickxsize.read())
        self.wroisz.setCurrentText(item)
        roisz_hbox.addWidget(self.wroisz)
        vbox.addLayout(roisz_hbox)
        roixy_hbox = QtGui.QHBoxLayout()
        self.roi_label = QtGui.QLabel("ROI")
        roixy_hbox.addWidget(self.roi_label)
        self.roix = QtGui.QLineEdit()
        roixy_hbox.addWidget(self.roix)
        self.roiy = QtGui.QLineEdit()
        roixy_hbox.addWidget(self.roiy)
        vbox.addLayout(roixy_hbox)
        cpr_hbox = QtGui.QHBoxLayout()
        self.cpr_label = QtGui.QLabel("CPR: ")
        cpr_hbox.addWidget(self.cpr_label)
        self.cpr = QtGui.QLineEdit()
        cpr_hbox.addWidget(self.cpr)
        self.coadd_label = QtGui.QLabel("Coadd: ")
        cpr_hbox.addWidget(self.coadd_label)
        self.coadd = QtGui.QLineEdit()
        cpr_hbox.addWidget(self.coadd)
        vbox.addLayout(cpr_hbox)
        apply_hbox = QtGui.QHBoxLayout()
        wapply = QtGui.QPushButton("Apply")
        wapply.clicked.connect(self.apply)
        apply_hbox.addWidget(wapply)
        wdismiss = QtGui.QPushButton("Dismiss")
        wdismiss.clicked.connect(self.dismiss)
        apply_hbox.addWidget(wdismiss)
        vbox.addLayout(apply_hbox)
        self.setLayout(vbox)

        screen = QDesktopWidget().screenGeometry()
        widget = self.geometry()
        x = screen.width() - widget.width()
        y = screen.height() - widget.height()
        self.move(x, y)
        self.resize(250, 0)

        self.threadpool = QtCore.QThreadPool()

        self.start_updating()

    def start_updating(self):
        self.updating = True
        updater = UpdateControlWindow(self.update)
        updater.signals.load.connect(self.update_gui)
        self.threadpool.start(updater)

    def update(self, file_callback):
        while self.updating:
            file_callback.emit()
            time.sleep(1)

    def stop_updating(self):
        self.updating = False

    def update_gui(self):
        roix = int(self.util.trickxsize.read())
        roiy = int(self.util.trickysize.read())
        roisz = str(self.util.trickxsize.read())
        centerx = int(int(self.util.trickxpos.read()) + roix/2)
        centery = int(int(self.util.trickypos.read()) + roiy/2)
        coadd = self.util.trknmad1.read()
        cpr = self.util.trkrocpr.read()
        self.roi_label.setText(f"ROI {centerx} {centery}")
        self.roisz_label.setText(f"ROI Size: {roisz}")
        self.coadd_label.setText(f"Coadd: {coadd}")
        self.cpr_label.setText(f"CPR: {cpr}")

    def apply(self):
        roix = int(self.util.trickxsize.read())
        roiy = int(self.util.trickysize.read())
        roisz = str(self.util.trickxsize.read())
        centerx = int(int(self.util.trickxpos.read()) + roix/2)
        centery = int(int(self.util.trickypos.read()) + roiy/2)
        coadd = int(self.util.trknmad1.read())
        cpr = int(self.util.trkrocpr.read())
        if self.roix.text() == '':
            set_roix = centerx
        else:
            set_roix = int(self.roix.text())
        if self.roiy.text() == '':
            set_roiy = centery
        else:
            set_roiy = int(self.roiy.text())
        if self.coadd.text() == '':
            set_coadd = coadd
        else:
            set_coadd = int(self.coadd.text())
        if self.cpr.text() == '':
            set_cpr = cpr
        else:
            set_cpr = int(self.cpr.text())
        set_roisz = None
        if self.wroisz.currentText() != roisz:
            set_roisz = int(self.wroisz.currentText())
        if (set_roix != centerx) or (set_roiy != centery) or (set_roisz != None):
            print("ROI change")
            xroi = set_roix
            yroi = set_roiy
            self.util.trickxpos.write(xroi)
            self.util.trickypos.write(yroi)
            distcoeff = np.zeros(20)
            rows = csv.reader(open('/usr/local/qfix/data/Trick/setup_files/TRICK_DistCoeff.dat','r'))
            for idx,row in enumerate(rows):
                distcoeff[idx] = float(row[0][5:])
            self.util.trk_putxy_spoc(xroi, yroi, distcoeff, roisz=set_roisz)
        if set_coadd != coadd:
            status = self.util.trk_set_coadds(set_coadd)
        if set_cpr != cpr:
            print("CPR change")
            status = self.util.trk_set_cpr(set_cpr)
        self.update_gui()


    def dismiss(self):
        self.close()




class FitsViewer(QtGui.QMainWindow):

    def __init__(self, logger):
        super(FitsViewer, self).__init__()
        self.logger = logger

        self.cachedFiles = None
        self.video = False

        self.util = Util()

        self.rawfile = ''
        self.mode = ''

        self.img = AstroImage()

        self.threadpool = QtCore.QThreadPool()

        self.iqcalc = iqcalc.IQCalc(self.logger)

        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params('zscale')
        fi.enable_autozoom('on')
        # fi.set_callback('drag-drop', self.drop_file)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        self.fitsimage = fi

        # enable some user interaction
        self.bd = fi.get_bindings()
        self.bd.enable_all(True)
        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(0, 0, 0, 0)
        vbox.setObjectName("vbox")
        status_hbox = QtGui.QHBoxLayout()
        status_hbox.setSpacing(1)
        status_hbox.setContentsMargins(QtCore.QMargins(5, 5, 5, 5))
        self.rmlabel = QtGui.QLabel("Read Mode: ")
        status_hbox.addWidget(self.rmlabel)
        self.rmstatus = QtGui.QLabel("")
        self.rmstatus.setStyleSheet("border: 1px solid black;")
        status_hbox.addWidget(self.rmstatus)
        self.tslabel  = QtGui.QLabel(" Status: ")
        status_hbox.addWidget(self.tslabel)
        self.trickstatus  = QtGui.QLabel("")
        self.trickstatus.setStyleSheet("border: 1px solid black;")
        status_hbox.addWidget(self.trickstatus)
        self.dtlabel = QtGui.QLabel(" DT: ")
        status_hbox.addWidget(self.dtlabel)
        self.dtstatus = QtGui.QLabel("")
        self.dtstatus.setStyleSheet("border: 1px solid black;")
        status_hbox.addWidget(self.dtstatus)
        vbox.addLayout(status_hbox)
        viewer_hbox = QtGui.QHBoxLayout()
        # viewer_hbox.setContentsMargins(QtCore.QMargins(48, 10, 48, 2))
        viewer_hbox.setObjectName("viewer_hbox")
        self.viewerLspacer = QtGui.QSpacerItem(40, 100, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        viewer_hbox.addItem(self.viewerLspacer)
        w = fi.get_widget()
        w.setMinimumSize(QtCore.QSize(240, 240))
        viewer_hbox.addWidget(w, stretch=1)
        self.viewerRspacer = QtGui.QSpacerItem(40, 100, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        viewer_hbox.addItem(self.viewerRspacer)
        hw = QtGui.QWidget()
        hw.setLayout(viewer_hbox)
        vbox.addWidget(hw)
        readout_hbox = QtGui.QHBoxLayout()
        readout_hbox.setObjectName("readout_hbox")
        self.readout = QtGui.QLabel("X:                 Y:                    Value:")
        self.readout.setObjectName("readout")
        readout_hbox.addWidget(self.readout)
        self.wdesaturate = QtGui.QPushButton("Desaturate")
        self.wdesaturate.setObjectName("wdesaturate")
        self.wdesaturate.setMaximumSize(QtCore.QSize(100, 100))
        self.wdesaturate.clicked.connect(self.util.desaturate)
        readout_hbox.addWidget(self.wdesaturate)
        self.wcut = QtGui.QComboBox()
        for name in fi.get_autocut_methods():
            self.wcut.addItem(name)
        item = str(self.util.trickxsize.read())
        self.wcut.setCurrentText(item)
        self.wcut.currentIndexChanged.connect(self.cut_change)
        readout_hbox.addWidget(self.wcut)
        self.wcut.setVisible(False)
        self.wcut.setCurrentText('median')
        self.wcolor = QtGui.QComboBox()
        for name in fi.get_color_algorithms():
            self.wcolor.addItem(name)
        self.wcolor.currentIndexChanged.connect(self.color_change)
        readout_hbox.addWidget(self.wcolor)
        self.wcolor.setVisible(False)
        hw = QtGui.QWidget()
        hw.setLayout(readout_hbox)
        vbox.addWidget(hw)
        roi_hbox = QtGui.QHBoxLayout()
        roi_hbox.setObjectName("roi_hbox")
        roi_hbox.setContentsMargins(QtCore.QMargins(10, 0, 10, 0))
        self.roi_info = QtGui.QLabel("ROI: ")
        self.roi_info.setObjectName("roi_info")
        roi_hbox.addWidget(self.roi_info)
        filter_hbox = QtGui.QHBoxLayout()
        filter_hbox.setObjectName("roi_hbox")
        self.vid_filter = QtGui.QLabel("Filter: ")
        self.vid_filter.setObjectName("vid_filter")
        filter_hbox.addWidget(self.vid_filter)
        self.wchangefilter = QtGui.QPushButton(f"{self.util.targname}")
        self.wchangefilter.setObjectName("wchangefilter")
        self.wchangefilter.clicked.connect(self.filter_popup)
        filter_hbox.addWidget(self.wchangefilter)
        self.winitfilter = QtGui.QPushButton("Init Wheel")
        self.winitfilter.setObjectName("winitfilter")
        self.winitfilter.clicked.connect(self.util.init_filter)
        filter_hbox.addWidget(self.winitfilter)
        hw = QtGui.QWidget()
        hw.setLayout(filter_hbox)
        roi_hbox.addWidget(hw)
        self.box_readout = QtGui.QLabel("Amplitude:                  FWHM: ")
        self.box_readout.setMinimumSize(QtCore.QSize(200, 0))
        self.box_readout.setObjectName("box_readout")
        self.box_readout.setVisible(False)
        roi_hbox.addWidget(self.box_readout)
        hw = QtGui.QWidget()
        hw.setLayout(roi_hbox)
        vbox.addWidget(hw)
        buttons_hbox = QtGui.QHBoxLayout()
        buttons_hbox.setObjectName("buttons_hbox")
        buttons_vbox_left = QtGui.QVBoxLayout()
        buttons_vbox_left.setContentsMargins(QtCore.QMargins(0, 0, 10, 0))
        buttons_vbox_left.setObjectName("buttons_vbox_left")
        self.wtakeff = QtGui.QPushButton("Take Full Frame")
        self.wtakeff.setObjectName("wtakeff")
        self.wtakeff.clicked.connect(self.take_ff)
        self.wtakeff.setVisible(False)
        buttons_vbox_left.addWidget(self.wtakeff)
        self.wvideomode = QtGui.QPushButton("Video Mode")
        self.wvideomode.setObjectName("wvideomode")
        self.wvideomode.clicked.connect(self.video_mode)
        self.wvideomode.setVisible(False)
        buttons_vbox_left.addWidget(self.wvideomode)
        self.wtrickcontrol = QtGui.QPushButton("Control")
        self.wtrickcontrol.setObjectName("wtrickcontrol")
        self.wtrickcontrol.clicked.connect(self.control_popup)
        buttons_vbox_left.addWidget(self.wtrickcontrol)
        self.wfullframemode = QtGui.QPushButton("Full Frame Mode")
        self.wfullframemode.setObjectName("wfullframemode")
        self.wfullframemode.clicked.connect(self.full_frame_mode)
        buttons_vbox_left.addWidget(self.wfullframemode)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_vbox_left)
        buttons_hbox.addWidget(hw)
        buttons_vbox_mid = QtGui.QVBoxLayout()
        buttons_vbox_mid.setContentsMargins(QtCore.QMargins(10, 0, 0, 0))
        buttons_vbox_mid.setObjectName("buttons_vbox_mid")
        self.winittrick = QtGui.QPushButton("Init Trick")
        self.winittrick.setObjectName("winittrick")
        self.winittrick.clicked.connect(self.init_trick)
        buttons_vbox_mid.addWidget(self.winittrick)
        self.wreboottrick = QtGui.QPushButton("Reboot Trick")
        self.wreboottrick.setObjectName("wreboottrick")
        self.wreboottrick.clicked.connect(self.reboot_trick)
        self.wreboottrick.setEnabled(False)
        buttons_vbox_mid.addWidget(self.wreboottrick)
        self.wsetroi = QtGui.QPushButton("Set ROI")
        self.wsetroi.setObjectName("wsetroi")
        self.wsetroi.clicked.connect(self.set_roi)
        self.wsetroi.setVisible(False)
        buttons_vbox_mid.addWidget(self.wsetroi)
        self.wsky = QtGui.QPushButton("Load Sky")
        self.wsky.clicked.connect(self.load_sky)
        self.wsky.setVisible(False)
        self.wsky.setVisible(False)
        buttons_vbox_mid.addWidget(self.wsky)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_vbox_mid)
        buttons_hbox.addWidget(hw)
        image_info_vbox = QtGui.QVBoxLayout()
        image_info_vbox.setObjectName("image_info_vbox")
        self.sky_info = QtGui.QLabel("Sky: ")
        self.sky_info.setObjectName("sky_info")
        self.sky_info.setVisible(False)
        image_info_vbox.addWidget(self.sky_info)
        self.filt_info = QtGui.QLabel("Filter: ")
        self.filt_info.setObjectName("filt_info")
        self.filt_info.setVisible(False)
        image_info_vbox.addWidget(self.filt_info)
        self.image_info = QtGui.QLabel("Image: ")
        self.image_info.setMinimumSize(QtCore.QSize(175, 0))
        self.image_info.setObjectName("image_info")
        self.image_info.setVisible(False)
        image_info_vbox.addWidget(self.image_info)
        hw = QtGui.QWidget()
        hw.setLayout(image_info_vbox)
        buttons_hbox.addWidget(hw)
        buttons_vbox_right = QtGui.QVBoxLayout()
        buttons_vbox_right.setContentsMargins(QtCore.QMargins(0, 0, 0, 0))
        buttons_vbox_right.setObjectName("buttons_vbox_right")
        self.wopen = QtGui.QPushButton("Open File")
        self.wopen.setObjectName("wopen")
        self.wopen.clicked.connect(self.open_file)
        self.wopen.setVisible(False)
        buttons_vbox_right.addWidget(self.wopen)
        self.wrestartvideo = QtGui.QPushButton("Restart Video")
        self.wrestartvideo.setObjectName("wrestartvideo")
        self.wrestartvideo.clicked.connect(self.restart_video)
        buttons_vbox_right.addWidget(self.wrestartvideo)
        self.wquit = QtGui.QPushButton("Quit")
        self.wquit.setObjectName("wquit")
        self.wquit.clicked.connect(self.quit)
        buttons_vbox_right.addWidget(self.wquit)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_vbox_right)
        buttons_hbox.setSpacing(0)
        buttons_hbox.addWidget(hw)
        hw = QtGui.QWidget()
        hw.setLayout(buttons_hbox)
        vbox.addWidget(hw)

        vbox.setSpacing(0)
        vw = QtGui.QWidget()
        self.setCentralWidget(vw)
        vw.setLayout(vbox)

        fi.set_autocut_params('minmax')
        fi.set_color_algorithm('log')

        fi.set_callback('cursor-changed', self.motion_cb)
        fi.add_callback('cursor-down', self.btndown)
        roix = int(self.util.trickxsize.read())
        roiy = int(self.util.trickysize.read())
        self.fitsimage.set_limits(((1,-0.5),(roix-2, roiy-0.5)),coord='data')

        self.recdc, self.compdc, self.crossdc = self.add_canvas()
        self.boxtag = "roi-box"
        self.picktag = "pick-box"
        self.crosstag = "vid-crosshair"

        fi.get_canvas().add(self.crossdc(float(roix)/2-0.5, float(roiy)/2-0.5, color='skyblue', text=""), tag=self.crosstag)

        self.trick_filters = ['Ks', 'H', 'Home', 'Open', 'Block', 'DISMISS']

        self.c = None #ControlWindow)

    def add_canvas(self, tag=None):
        # add a canvas to the view
        my_canvas = self.fitsimage.get_canvas()
        RecCanvas = my_canvas.get_draw_class('rectangle')
        CompCanvas = my_canvas.get_draw_class('compass')
        CrossCanvas = my_canvas.get_draw_class('crosshair')
        return RecCanvas, CompCanvas, CrossCanvas


    def cut_change(self):
        self.fitsimage.set_autocut_params(self.wcut.currentText())

    def color_change(self):
        self.fitsimage.set_color_algorithm(self.wcolor.currentText())

    def motion_cb(self, viewer, button, data_x, data_y):

        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = int(viewer.get_data(int(data_x + 0.5), int(data_y + 0.5)))

        except Exception:
            value = None

        fits_x, fits_y = data_x, data_y

        # Calculate WCS RA
        try:
            # NOTE: image function operates on DATA space coords
            image = viewer.get_image()
            if image is None:
                # No image loaded
                return
            ra_txt, dec_txt = image.pixtoradec(fits_x, fits_y,
                                               format='str', coords='fits')
        except Exception as e:
            self.logger.warning("Bad coordinate conversion: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'

        if (fits_x+0.5 > 2048 or fits_x+0.5 <0) or (fits_y+0.5 > 2048 or fits_y+0.5 <0):
            text = "X: Y:  Value:"
            self.readout.setText(text)
        else:
            text = f"X: {int(fits_x+1.5) + int(self.util.trickxpos.read())} Y: {int(fits_y+1.5) + int(self.util.trickypos.read())}  Value: {value}"
            self.readout.setText(text)

    def quit(self, *args):
        self.logger.info("Attempting to shut down the application...")
        self.stop_video()
        self.stop_scan()
        time.sleep(1)
        self.threadpool = False
        QtGui.QApplication.instance().quit()

    ##TODO verify init trick, restart video, stop video, reboot trick, change_filter, desaturate


    def filter_popup(self):
        msg = QtGui.QMessageBox()
        msg.setWindowTitle("Change Trick Filter")
        msg.setText("Change filter to:")
        msg.setIcon(QtGui.QMessageBox.Information)
        for filter in self.trick_filters:
            msg.addButton(f'{filter}', QtGui.QMessageBox.YesRole)
        msg.buttonClicked.connect(self.util.change_filter)
        x = msg.exec_()

    def tiptilt_popup(self):
        msg = QtGui.QMessageBox()
        msg.setWindowTitle("Tip-tilt loop closed on TRICK")
        msg.setText("Open the loop or switch to STRAP before changing filter")
        msg.setIcon(QtGui.QMessageBox.Critical)
        y = msg.exec_()

    def init_trick(self):
        print("Initing TRICK")
        self.wquit.setEnabled(False)
        self.util.stopex.write(1)
        time.sleep(3)
        self.util.init.write(1)
        self.util.sampmode.write(5)
        self.util.cyclespr.write(50)
        self.util.cdsmode.write(1)
        self.util.readmode.write(3)
        self.util.go.write(1)
        time.sleep(3)
        if self.util.ops == "MGAO":
            self.util.trkenapx.write(0)
            self.util.trkfpspx.write('Passive')
        self.util.trkstop.write(1)
        time.sleep(1)
        if self.util.ops == "MGAO":
            self.util.trkfpspx.write('1 second')
            self.util.trkenapx.write(1)
        self.util.trkstsx.write(1)
        self.video = True
        video = Video(self.display_video)
        video.signals.load.connect(self.show_images)
        self.threadpool.start(video)
        self.wquit.setEnabled(True)
        self.mode = 'video'
        self.wrestartvideo.setEnabled(True)

    def restart_video(self):
        self.util.stopex.write(1)
        time.sleep(3)
        self.util.cdsmode.write(1)
        self.util.readmode.write(3)
        self.util.go.write(1)
        if self.util.ops == "MGAO":
            self.util.trkenapx.write(0)
            self.util.trkfpspx.write('Passive')
        self.util.trkstop.write(1)
        time.sleep(1)
        if self.util.ops == "MGAO":
            self.util.trkfpspx.write('1 second')
            self.util.trkenapx.write(1)
        self.util.trkstsx.write(1)
        self.start_video()


    ##TODO remove this when switching back to video mode is replaced with restart_video
    def start_video(self):
        self.video = True
        video = Video(self.display_video)
        video.signals.load.connect(self.show_images)
        left, right, up, down = self.getROI()
        self.threadpool.start(video)
        self.mode = 'video'

    #TODO make this actually stop video mode
    def stop_video(self):
        self.video = False

    def reboot_trick(self):
        print("Rebooting TRICK...")

    def display_video(self, file_callback):
        while self.video:
            file_callback.emit(self.util.trkro1px, self.util.trkro1ff, self.util.trkro1bg)
            left, right, up, down = self.getROI()
            time.sleep(1)

    def control_popup(self):
        self.c = ControlWindow()
        self.c.show()


    def show_images(self, pix, ff, bg):
        image = self.pixels_to_image(pix, ff, bg)
        roix = int(self.util.trickxsize.read())
        roiy = int(self.util.trickysize.read())
        self.img.load_data(image)
        self.fitsimage.set_image(self.img)
        self.fitsimage.set_limits(((1,-0.5),(roix-2, roiy-0.5)),coord='data')
        self.fitsimage.get_canvas().get_object_by_tag(self.crosstag)
        self.fitsimage.get_canvas().delete_object_by_tag(self.crosstag)
        self.fitsimage.get_canvas().add(self.crossdc(float(roix)/2-0.5, float(roiy)/2-0.5, color='blue', text=""), tag=self.crosstag)
        self.wchangefilter.setText(str(self.util.targname))
        self.dtstatus.setText(str(self.util.dtsensor))
        ts = str(self.util.trkstat)
        if ts == '2nd-channel video in progress':
            ts = '2C Video'
        self.trickstatus.setText(ts)
        m = self.util.readmode.read()
        if m == "3":
            self.rmstatus.setText("Video")
        elif m == "2":
            self.rmstatus.setText("Filmstrip")
        else:
            self.rmstatus.setText("FFrame")
        self.resize(240, 300)

    def pixels_to_image(self, pix, ff, bg):
        lst = str(pix).strip().replace(':', '').split()
        lst_ff = str(ff).strip().replace(':', '').split()
        lst_bg = str(bg).strip().replace(':', '').split()
        roix = int(self.util.trickxsize.read())
        roiy = int(self.util.trickysize.read())
        roi_size = roix * roiy
        pixelvalues = np.array(lst[1::2],dtype=float) # take every second value, since the first value is the pixel number
        pixelvalues= pixelvalues[:roi_size]
        pixelvalues_ff = np.array(lst_ff[1::2],dtype=float)
        pixelvalues_ff= pixelvalues_ff[:roi_size]
        pixelvalues_bg = np.array(lst_bg[1::2],dtype=float)
        pixelvalues_bg= pixelvalues_bg[:roi_size]
        dims = int(np.sqrt(pixelvalues.shape))
        image = np.reshape(pixelvalues,(dims,dims))
        ff = np.reshape(pixelvalues_ff,(dims,dims))
        bg = np.reshape(pixelvalues_bg,(dims,dims))
        new_image = image * ff - bg
        roi_image = new_image[0:roix, 0:roiy]
        return(roi_image)

    def full_frame_mode(self):
        if self.c != None:
            self.c = None
        self.wfullframemode.setVisible(False)
        self.wvideomode.setVisible(True)
        self.stop_video()
        self.fitsimage.clear()
        self.fitsimage.get_canvas().get_object_by_tag(self.crosstag)
        self.fitsimage.get_canvas().delete_object_by_tag(self.crosstag)
        self.resize(500, 700)
        self.viewerLspacer.changeSize(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.viewerRspacer.changeSize(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.fitsimage.get_widget().setMinimumSize(QtCore.QSize(512, 512))
        self.fitsimage.set_limits(None)
        self.readout.setMinimumSize(QtCore.QSize(240, 0))
        self.fitsimage.set_autocut_params('median')
        self.fitsimage.set_color_algorithm('linear')
        self.wdesaturate.setVisible(False)
        self.vid_filter.setVisible(False)
        self.wchangefilter.setVisible(False)
        self.winitfilter.setVisible(False)
        self.winittrick.setVisible(False)
        self.wrestartvideo.setVisible(False)
        self.wreboottrick.setVisible(False)
        self.wreboottrick.setVisible(False)
        self.wtrickcontrol.setVisible(False)
        self.rmstatus.setVisible(False)
        self.trickstatus.setVisible(False)
        self.dtstatus.setVisible(False)
        self.rmlabel.setVisible(False)
        self.tslabel.setVisible(False)
        self.dtlabel.setVisible(False)
        self.wcolor.setVisible(True)
        self.wcut.setVisible(True)
        self.wsky.setVisible(True)
        self.wtakeff.setVisible(True)
        self.wsetroi.setVisible(True)
        self.wsetroi.setEnabled(False)
        self.image_info.setVisible(True)
        self.sky_info.setVisible(True)
        self.filt_info.setVisible(True)
        self.wopen.setVisible(True)
        self.box_readout.setVisible(True)
        self.start_scan()
        self.mode = 'fullframe'

    def video_mode(self):
        self.wvideomode.setVisible(False)
        self.wfullframemode.setVisible(True)
        self.stop_scan()
        self.fitsimage.clear()
        self.fitsimage.rotate(0)
        self.fitsimage.get_widget().setMinimumSize(QtCore.QSize(240, 240))
        self.viewerLspacer.changeSize(40, 100, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.viewerRspacer.changeSize(40, 100, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        roix = int(self.util.trickxsize.read())
        roiy = int(self.util.trickysize.read())
        self.fitsimage.set_limits(((1,-0.5),(roix-2, roiy-0.5)),coord='data')
        self.fitsimage.get_canvas().add(self.crossdc(float(roix)/2-0.5, float(roiy)/2-0.5, color='blue', text=""), tag=self.crosstag)
        self.fitsimage.set_autocut_params('minmax')
        self.fitsimage.set_color_algorithm('log')
        self.wdesaturate.setVisible(True)
        self.vid_filter.setVisible(True)
        self.wchangefilter.setVisible(True)
        self.winitfilter.setVisible(True)
        self.winittrick.setVisible(True)
        self.wrestartvideo.setVisible(True)
        self.wreboottrick.setVisible(True)
        self.wquit.setVisible(True)
        self.wtrickcontrol.setVisible(True)
        self.rmstatus.setVisible(True)
        self.trickstatus.setVisible(True)
        self.dtstatus.setVisible(True)
        self.rmlabel.setVisible(True)
        self.tslabel.setVisible(True)
        self.dtlabel.setVisible(True)
        self.wcolor.setVisible(False)
        self.wcut.setVisible(False)
        self.wsky.setVisible(False)
        self.wtakeff.setVisible(False)
        self.wsetroi.setVisible(False)
        self.image_info.setVisible(False)
        self.sky_info.setVisible(False)
        self.filt_info.setVisible(False)
        self.box_readout.setVisible(False)
        self.wopen.setVisible(False)
        self.readout.setMinimumSize(QtCore.QSize(0, 0))
        self.resize(240, 300)
        #TODO replace this with restart video
        self.start_video()

    ##Full frame stuff
    def start_scan(self):
        self.scanning = True
        hdu = fits.PrimaryHDU()
        try:
            hdu.writeto('procImage.fits')
        except OSError:
            os.remove('procImage.fits')
            hdu.writeto('procImage.fits')
        self.cachedFiles = self.walkDirectory()
        scanner = Scanner(self.scan)
        scanner.signals.load.connect(self.processData)
        self.threadpool.start(scanner)

    def stop_scan(self):
        self.scanning = False

    def load_file(self, filepath):
        image = load_data(filepath, logger=self.logger)
        self.fitsimage.set_image(image)
        # self.setWindowTitle(filepath)
        left, right, up, down = self.getROI()
        try:
            self.fitsimage.get_canvas().get_object_by_tag(self.boxtag)
            self.fitsimage.get_canvas().delete_object_by_tag(self.boxtag)
            self.box = self.recdc(left, down, right, up, color='green')
            self.fitsimage.get_canvas().add(self.box, tag=self.boxtag, redraw=True)
        except KeyError:
            self.box = self.recdc(left, down, right, up, color='green')
            self.fitsimage.get_canvas().add(self.box, tag=self.boxtag, redraw=True)
        try:
            self.fitsimage.get_canvas().get_object_by_tag(self.picktag)
            self.fitsimage.get_canvas().delete_object_by_tag(self.picktag)
        except KeyError:
            pass
        width, height = image.get_size()
        data_x, data_y = width / 2.0, height / 2.0
        # x, y = self.fitsimage.get_canvas_xy(data_x, data_y)
        radius = float(max(width, height)) / 20
        self.fitsimage.get_canvas().add(self.compdc(data_x, data_y, radius, color='skyblue',
                                       fontsize=8))
        self.bd._orient(self.fitsimage, righthand=False, msg=True)

    def open_file(self):
        res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
                                                str(self.nightpath()))
        print(res)
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.processData(fileName)

    def load_sky(self):
        res = QtGui.QFileDialog.getOpenFileName(self, "Open Sky file",
                                                str(self.nightpath()))
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.subtract_sky(fileName)

    def subtract_sky(self, filename):
        skyname, skyheader, skyfitsData, skyfilter = self.addWcs(filename)
        name, header, fitsData, filter = self.addWcs(self.rawfile)
        with_sky = fitsData - skyfitsData
        mask = fits.getdata('/kroot/rel/ao/qfix/data/Trick/BadPix_1014Hz.fits', ext=0)
        text = f"Sky: {skyname}"
        self.sky_info.setText(text)
        self.load_file(self.writeFits(header, np.multiply(with_sky, mask)))

    def take_ff(self):
        self.wtakeff.setEnabled(False)
        self.util.stopex.write(1)
        time.sleep(3)
        self.util.timfile.write('/kroot/rel/default/data/nirtts_rel_v1.6.lod')
        self.util.cdsmode.write(1)
        self.util.readmode.write(1)
        self.util.itime.write(0)
        self.util.getkw.write(1)
        self.util.getdcskw.write(1)
        self.util.getaokw.write(1)
        self.util.go.write(1)

    def set_roi(self):
        xroi = self.xclick+1.0
        yroi = self.yclick+1.0
        self.util.trickxpos.write(xroi)
        self.util.trickypos.write(yroi)
        distcoeff = np.zeros(20)
        rows = csv.reader(open('/usr/local/qfix/data/Trick/setup_files/TRICK_DistCoeff.dat','r'))
        for idx,row in enumerate(rows):
            distcoeff[idx] = float(row[0][5:])
        self.util.trk_putxy_spoc(xroi, yroi, distcoeff, roisz=None)
        left, right, up, down = self.getROI()
        self.fitsimage.get_canvas().get_object_by_tag(self.boxtag)
        self.fitsimage.get_canvas().delete_object_by_tag(self.boxtag)
        self.box = self.recdc(left, down, right, up, color='green')
        self.fitsimage.get_canvas().add(self.box, tag=self.boxtag, redraw=True)
        self.wsetroi.setEnabled(False)

    ##Start of image find and processing code

    def scan(self, file_callback):
        while self.scanning:
            hasNewFiles, files, self.cachedFiles = self.updateFileCache(self.cachedFiles)
            if hasNewFiles and ('.fits' in files[0] or '.fits.gz' in files[0]):
                print("New Image Detected!")
                filen = files[0]
                self.waitForFileToBeUnlocked(filen, 1)
                file_callback.emit(filen)
            time.sleep(1)

    def walkDirectory(self):
        directory = self.nightpath()
        return [abspath(join(directory, f)) for f in listdir(directory) if isfile(join(directory, f))]

    def updateFileCache(self, cachedFiles):
        updatedFileList = self.walkDirectory()
        filtered = [i for i in updatedFileList if not i in cachedFiles]
        cachedFiles = updatedFileList
        return len(filtered) > 0, filtered, cachedFiles

    def fileIsCurrentlyLocked(self, filepath):
        locked = None
        hdulist = None
        file_object = None
        if os.path.exists(filepath):
            try:
                print("Trying to open %s." % filepath)
                #time.sleep(15) #place holder if OSError catch doesn't work
                hdulist = fits.open(filepath)

                file_object = np.sum([1 for hdu in hdulist if type(hdu) in
                        	[fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU]
                        	and hdu.data is not None])
                if file_object:
                    locked = False

            except TypeError:
                locked = True

            except OSError:
                locked = True

            finally:
                if file_object:
                    hdulist.close()

        else:
            print("%s not found." % filepath)

        return locked


    #    Checks if the files are ready.
    #    For a file to be ready it must exist and can be opened in append mode.

    def waitForFileToBeUnlocked(self, filename, wait_time):
        # if the file doesn't exist, wait wait_time seconds and try again until it's found
        while not os.path.exists(filename):
            print("%s hasn't arrived. Waiting %s seconds." % (filename, wait_time))
            time.sleep(wait_time)

        # if the file exists but locked, wait wait_time seconds and check
        # again until it's no longer locked by another process
        while self.fileIsCurrentlyLocked(filename):
            time.sleep(wait_time)

    def getROI(self):
        left = int(self.util.trickxpos.read()) - int(self.util.trickxsize.read())*3
        right = int(self.util.trickxpos.read()) + int(self.util.trickxsize.read())*3
        up = int(self.util.trickypos.read()) - int(self.util.trickysize.read())*3
        down = int(self.util.trickypos.read()) + int(self.util.trickysize.read())*3
        roix = int(self.util.trickxsize.read())
        roiy = int(self.util.trickysize.read())
        centerx = int(int(self.util.trickxpos.read()) + roix/2)
        centery = int(int(self.util.trickypos.read()) + roiy/2)
        self.roi_info.setText(f"ROI {centerx} {centery}")
        return left, right, up, down

    def nightpath(self):
        nightly = Path('/net/k1aoserver/k1aodata/nightly')
        date = datetime.datetime.utcnow()
        year, month, day = str(date.strftime("%y")), str(date.strftime("%m")), str(date.strftime("%d"))
        nightly = nightly / year / month / day / 'Trick'
        return nightly

    def processData(self, filename):
        self.rawfile = filename
        name, header, fitsData, filter = self.addWcs(filename)
        mask = fits.getdata('/kroot/rel/ao/qfix/data/Trick/BadPix_1014Hz.fits', ext=0)
        if filter == 'H':
            background = fits.getdata('/kroot/rel/ao/qfix/data/Trick/sky_H.fits')
            self.sky_info.setText('Sky: sky_H.fits')
        else:
            background = fits.getdata('/kroot/rel/ao/qfix/data/Trick/sky_Ks.fits')
            self.sky_info.setText('sky_Ks.fits')
        subtracted_data = fitsData-background
        self.load_file(self.writeFits(header, np.multiply(subtracted_data, mask)))
        text = f"Image: {name}"
        self.image_info.setText(text)
        text = f"Filter: {filter}"
        self.filt_info.setText(text)
        self.wsky.setEnabled(True)
        self.wtakeff.setEnabled(True)

    def addWcs(self, filen):
        w = wcs.WCS(naxis=2)
        fitsData = fits.getdata(filen, ext=0)
        header = fits.getheader(filen)
        ht, wd = fitsData.shape[:2]
        y = ht//2
        x = wd//2
        name = header['DATAFILE']
        ra = float(header['RA'])
        dec = float(header['DEC'])
        try:
            rot = float(header['ROTPOSN'])
        except VerifyError:
            self.rotator_invalid_popup()
            print("Invalid rotator angle in header, unable to rotate image")
            rot = 0.0
        filter = header['TRFWNAME']
        w.wcs.crpix = [y, x]
        w.wcs.cdelt = np.array([-0.05, 0.05])
        w.wcs.crota = np.array([0.05, rot])
        w.wcs.crval = [ra, dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        pixcrd = np.array([[0, 0], [24, 38], [45, 98]], dtype=np.float64)
        world = w.wcs_pix2world(pixcrd, 0)
        # Convert the same coordinates back to pixel coordinates.
        pixcrd2 = w.wcs_world2pix(world, 0)
        # These should be the same as the original pixel coordinates, modulo
        # some floating-point error.
        assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6
        # Now, write out the WCS object as a FITS header
        header = w.to_header()
        return name, header, fitsData, filter

    def rotator_invalid_popup(self):
        msg = QtGui.QMessageBox()
        msg.setWindowTitle("Invalid Header")
        msg.setText("Invalid rotator angle in header, unable to rotate image")
        msg.setIcon(QtGui.QMessageBox.Critical)
        y = msg.exec_()

    def writeFits(self, headerinfo, image_data):
        hdu = fits.PrimaryHDU(header=headerinfo, data=image_data)
        filename = 'procImage.fits'
        try:
            hdu.writeto(filename)
        except OSError:
            os.remove(filename)
            hdu.writeto(filename)
        return filename

    ##Find star stuff
    def cutdetail(self, image, shape_obj):
        view, mask = image.get_shape_view(shape_obj)

        data = image._slice(view)

        y1, y2 = view[0].start, view[0].stop
        x1, x2 = view[1].start, view[1].stop

        # mask non-containing members
        mdata = np.ma.array(data, mask=np.logical_not(mask))

        return x1, y1, x2, y2, mdata

    def findstar(self):
        image = self.fitsimage.get_image()
        obj = self.pickbox
        shape = obj
        x1, y1, x2, y2, data = self.cutdetail(image, shape)
        ht, wd = data.shape[:2]
        xc, yc = wd // 2, ht // 2
        radius = min(xc, yc)
        peaks = [(xc, yc)]
        peaks = self.iqcalc.find_bright_peaks(data,
                                              threshold=None,
                                              radius=radius)

        xc, yc = peaks[0]
        xc += 1
        yc += 1
        return int(xc), int(yc), data

    def fitstars(self, y_line):
        x = np.linspace(-30, 30, 60)
        model_gauss = models.Gaussian1D()
        fitter_gauss = fitting.LevMarLSQFitter()
        # gx = fitter_gauss(model_gauss, x, x_line)
        gy = fitter_gauss(model_gauss, x, y_line)
        # amplitude = (gx.amplitude.value+gy.amplitude.value)/2
        amplitude = gy.amplitude.value
        # fwhm = ((gx.stddev.value+gy.stddev.value)/2)*0.118 #2.355*PixelScale
        # fwhm = (gy.stddev.value)*0.118 #2.355*PixelScale
        fwhm = (gy.stddev.value)*2.355 #pixels instead of arcseconds
        return amplitude, fwhm

    def pickstar(self, viewer):
        warnings.warn("deprecated", DeprecationWarning)
        try:
            self.fitsimage.get_canvas().get_object_by_tag(self.picktag)
            self.fitsimage.get_canvas().delete_object_by_tag(self.picktag)
            self.pickbox = self.recdc(self.xclick-(int(self.util.trickxsize)*3), self.yclick-(int(self.util.trickysize)*3), self.xclick+(int(self.util.trickxsize)*3), self.yclick+(int(self.util.trickysize)*3), color='red')
            self.fitsimage.get_canvas().add(self.pickbox, tag=self.picktag, redraw=True)
        except KeyError:
            self.pickbox = self.recdc(self.xclick-30, self.yclick-30, self.xclick+30, self.yclick+30, color='red')
            self.fitsimage.get_canvas().add(self.pickbox, tag=self.picktag, redraw=True)
        image = self.fitsimage.get_image()
        try:
            xc, yc, data = self.findstar()
            # x_line = data[40-yc, 0:40] doesn't work well for some reason
            y_line = data[0:60, xc]
            # amplitude, fwhm = self.fitstars(x_line, y_line)
            amplitude, fwhm = self.fitstars(y_line)
            text = f"Amplitude: {amplitude:.2f} FWHM: {fwhm:.2f}"
            self.box_readout.setText(text)
        except IndexError:
            text = "Amplitude: N/A FWHM: N/A"
            self.box_readout.setText(text)


    def btndown(self, canvas, event, data_x, data_y):
        self.xclick = data_x
        self.yclick = data_y
        ##todo video mode adjusting ROI
        if self.mode == "video":
            xroi = int(1.5 + self.xclick) + int(self.util.trickxpos.read())
            yroi = int(1.5 + self.yclick) + int(self.util.trickypos.read())
            distcoeff = np.zeros(20)
            rows = csv.reader(open('/usr/local/qfix/data/Trick/setup_files/TRICK_DistCoeff.dat','r'))
            for idx,row in enumerate(rows):
                distcoeff[idx] = float(row[0][5:])
            self.util.trk_putxy_spoc(xroi, yroi, distcoeff, roisz=None)
        else:
            self.wsetroi.setEnabled(True)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.pickstar(self.fitsimage)



def main():
    ##Write dummy file so walkDirectory caches it in the beginning

    app = QtGui.QApplication([])

    # ginga needs a logger.
    # If you don't want to log anything you can create a null logger by
    # using null=True in this call instead of log_stderr=True
    logger = log.get_logger("TrickManager", log_stderr=True, level=40)

    w = FitsViewer(logger)
    w.show()
    app.setActiveWindow(w)
    w.raise_()
    w.activateWindow()
    w.start_video()
    sys.exit(app.exec_())

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)

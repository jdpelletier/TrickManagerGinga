import os, time, sys, threading, math
from os import listdir
from os.path import abspath, isfile, join
from pathlib import Path
import datetime

import numpy as np
from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.modeling import models, fitting
import PIL.Image as PILimage

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
    load = QtCore.Signal(object)

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


class FitsViewer(QtGui.QMainWindow):

    def __init__(self, logger, MainWindow):
        super(FitsViewer, self).__init__()
        self.logger = logger

        self.cachedFiles = None
        self.video = False
        #KTL stuff
        #Cache KTL keywords
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
        self.progress = ktl.cache('tds', 'progress')
        self.roipixels = ktl.cache('ao', 'TRKRO1PX')
        self.roipixels.monitor()
        self.cyclespr = ktl.cache('tds', 'cyclespr')
        #TODO add back in or delete
        # self.trkfpspx = ktl.cache('ao', 'trkfpspx')
        # self.trkenapx = ktl.cache('ao', 'trkenapx')
        self.trkstop = ktl.cache('trick', 'trkstop')
        self.trkstsx = ktl.cache('trick', 'trkstsx')

        self.rawfile = ''
        self.mode = ''

        self.img = AstroImage()

        self.threadpool = QtCore.QThreadPool()

        self.iqcalc = iqcalc.IQCalc(self.logger)

        MainWindow.setObjectName("MainWindow")
        self.MainWindow = MainWindow
        self.MainWindow.resize(330, 370)

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
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.hbox_layout = QtGui.QWidget(self.centralwidget)
        self.hbox_layout.setGeometry(QtCore.QRect(0, 0, 330, 240))
        self.hbox_layout.setObjectName("hbox_layout")
        self.display_hbox = QtGui.QHBoxLayout(self.hbox_layout)
        self.display_hbox.setContentsMargins(0, 0, 0, 0)
        self.display_hbox.setObjectName("display_hbox")
        w = fi.get_widget()
        w.setObjectName("w")
        self.display_hbox.addWidget(w, stretch=1)
        self.readout = QtGui.QLabel(self.centralwidget)
        self.readout.setGeometry(QtCore.QRect(30, 240, 251, 16))
        self.readout.setObjectName("readout")
        self.roi_info = QtGui.QLabel(self.centralwidget)
        self.roi_info.setGeometry(QtCore.QRect(20, 270, 171, 20))
        self.roi_info.setObjectName("roi_info")
        self.image_info = QtGui.QLabel(self.centralwidget)
        self.image_info.setGeometry(QtCore.QRect(480, 610, 111, 16))
        self.image_info.setObjectName("image_info")
        self.image_info.setVisible(False)
        self.sky_info = QtGui.QLabel(self.centralwidget)
        self.sky_info.setGeometry(QtCore.QRect(480, 570, 111, 16))
        self.sky_info.setObjectName("sky_info")
        self.sky_info.setVisible(False)
        self.filt_info = QtGui.QLabel(self.centralwidget)
        self.filt_info.setGeometry(QtCore.QRect(480, 590, 111, 16))
        self.filt_info.setObjectName("filt_info")
        self.filt_info.setVisible(False)
        self.ff_readout = QtGui.QLabel(self.centralwidget)
        self.ff_readout.setGeometry(QtCore.QRect(210, 520, 281, 16))
        self.ff_readout.setObjectName("readout")
        self.wcut = QtGui.QComboBox(self.centralwidget)
        self.wcut.setGeometry(QtCore.QRect(540, 520, 73, 22))
        self.wcut.setObjectName("wcut")
        for name in fi.get_autocut_methods():
            self.wcut.addItem(name)
        self.wcut.currentIndexChanged.connect(self.cut_change)
        self.wcolor = QtGui.QComboBox(self.centralwidget)
        self.wcolor.setGeometry(QtCore.QRect(620, 520, 73, 22))
        self.wcolor.setObjectName("wcolor")
        for name in fi.get_color_algorithms():
            self.wcolor.addItem(name)
        self.wcolor.currentIndexChanged.connect(self.color_change)
        self.wrestartvideo = QtGui.QPushButton(self.centralwidget)
        self.wrestartvideo.setGeometry(QtCore.QRect(230, 280, 93, 28))
        self.wrestartvideo.setObjectName("wrestartvideo")
        self.wfullframemode = QtGui.QPushButton(self.centralwidget)
        self.wfullframemode.setGeometry(QtCore.QRect(10, 310, 111, 28))
        self.wfullframemode.setObjectName("wfullframemode")
        self.wfullframemode.clicked.connect(self.full_frame_mode)
        self.winittrick = QtGui.QPushButton(self.centralwidget)
        self.winittrick.setGeometry(QtCore.QRect(130, 280, 93, 28))
        self.winittrick.setObjectName("winittrick")
        self.winittrick.clicked.connect(self.init_trick)
        self.wreboottrick = QtGui.QPushButton(self.centralwidget)
        self.wreboottrick.setGeometry(QtCore.QRect(130, 310, 93, 28))
        self.wreboottrick.setObjectName("wreboottrick")
        self.wreboottrick.clicked.connect(self.reboot_trick)
        self.wsky = QtGui.QPushButton(self.centralwidget)
        self.wsky.setGeometry(QtCore.QRect(330, 610, 93, 28))
        self.wsky.setObjectName("wsky")
        self.wsky.setVisible(False)
        self.wsky.clicked.connect(self.load_sky)
        self.wsetroi = QtGui.QPushButton(self.centralwidget)
        self.wsetroi.setGeometry(QtCore.QRect(330, 570, 93, 28))
        self.wsetroi.setObjectName("wsetroi")
        self.wsetroi.setVisible(False)
        self.wsetroi.clicked.connect(self.set_roi)
        self.wvideomode = QtGui.QPushButton(self.centralwidget)
        self.wvideomode.setGeometry(QtCore.QRect(10, 610, 91, 28))
        self.wvideomode.setObjectName("wvideomode")
        self.wvideomode.setVisible(False)
        self.wvideomode.clicked.connect(self.video_mode)
        self.wtakeff = QtGui.QPushButton(self.centralwidget)
        self.wtakeff.setGeometry(QtCore.QRect(10, 570, 111, 28))
        self.wtakeff.setObjectName("wtakeff")
        self.wtakeff.setVisible(False)
        self.wtakeff.clicked.connect(self.take_ff)
        self.box_readout = QtGui.QLabel(self.centralwidget)
        self.box_readout.setGeometry(QtCore.QRect(20, 550, 281, 16))
        self.box_readout.setObjectName("box_readout")
        self.wquit = QtGui.QPushButton(self.centralwidget)
        self.wquit.setGeometry(QtCore.QRect(230, 310, 93, 28))
        self.wquit.setObjectName("wquit")
        self.wquit.clicked.connect(self.quit)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 700, 26))
        self.menubar.setObjectName("menubar")
        self.menubar.setVisible(False)
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.wopenfile = QtGui.QAction(MainWindow)
        self.wopenfile.setObjectName("wopenfile")
        self.wopenfile.triggered.connect(self.open_file)
        self.wmquit = QtGui.QAction(MainWindow)
        self.wmquit.setObjectName("wmquit")
        self.wmquit.triggered.connect(self.quit)
        self.menuFile.addAction(self.wopenfile)
        self.menuFile.addAction(self.wmquit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        fi.set_callback('cursor-changed', self.motion_cb)
        fi.add_callback('cursor-down', self.btndown)

        self.recdc, self.compdc = self.add_canvas()
        self.boxtag = "roi-box"
        self.picktag = "pick-box"

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.readout.setText(_translate("MainWindow", "X:                 Y:                    Value:"))
        self.roi_info.setText(_translate("MainWindow", "ROI "))
        self.image_info.setText(_translate("MainWindow", "Image: "))
        self.sky_info.setText(_translate("MainWindow", "Sky:"))
        self.filt_info.setText(_translate("MainWindow", "Filter: "))
        self.wrestartvideo.setText(_translate("MainWindow", "Restart Video"))
        self.wfullframemode.setText(_translate("MainWindow", "Full Frame Mode"))
        self.winittrick.setText(_translate("MainWindow", "Init Trick"))
        self.wreboottrick.setText(_translate("MainWindow", "Reboot Trick"))
        self.wsky.setText(_translate("MainWindow", "Load Sky"))
        self.wsetroi.setText(_translate("MainWindow", "Set ROI"))
        self.wvideomode.setText(_translate("MainWindow", "Video Mode"))
        self.wtakeff.setText(_translate("MainWindow", "Take Full Frame"))
        self.box_readout.setText(_translate("MainWindow", "Amplitude:                  FWHM: "))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.wopenfile.setText(_translate("MainWindow", "Open"))
        self.wquit.setText(_translate("MainWindow", "Quit"))
        self.wmquit.setText(_translate("MainWindow", "Quit"))
        self.wvideomode.setVisible(False)
        self.wsky.setVisible(False)
        self.wtakeff.setVisible(False)
        self.wsetroi.setVisible(False)
        self.image_info.setVisible(False)
        self.sky_info.setVisible(False)
        self.filt_info.setVisible(False)
        self.wvideomode.setEnabled(False)
        self.menubar.setVisible(False)
        self.wcolor.setVisible(False)
        self.wcut.setVisible(False)

    def add_canvas(self, tag=None):
        # add a canvas to the view
        my_canvas = self.fitsimage.get_canvas()
        RecCanvas = my_canvas.get_draw_class('rectangle')
        CompCanvas = my_canvas.get_draw_class('compass')
        return RecCanvas, CompCanvas


    def cut_change(self):
        self.fitsimage.set_autocut_params(self.wcut.currentText())

    def color_change(self):
        self.fitsimage.set_color_algorithm(self.wcolor.currentText())

    def motion_cb(self, viewer, button, data_x, data_y):

        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = viewer.get_data(int(data_x + 0.5), int(data_y + 0.5))

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

        text = "X: %.2f  Y: %.2f  Value: %s" % (fits_x, fits_y, value)
        self.readout.setText(text)

    def quit(self, *args):
        self.logger.info("Attempting to shut down the application...")
        self.stop_video()
        self.stop_scan()
        time.sleep(2)
        self.threadpool = False
        QtGui.QApplication.instance().quit()

    ##TODO verify init trick, restart video, stop video, reboot trick

    def init_trick(self):
        print("Initing TRICK")
        self.wquit.setEnabled(False)
        self.stopex.write(1)
        time.sleep(3)
        self.init.write(1)
        self.sampmode.write(5)
        self.cyclespr.write(50)
        self.cdsmode.write(1)
        self.readmode.write(3)
        self.go.write(1)
        time.sleep(3)
        # self.trkenapx.write(0)
        # self.trkfpspx.write('Passive')
        self.trkstop.write(1)
        time.sleep(1)
        # self.trkfpspx.write('1 second')
        # self.trkenapx.write(1)
        self.trkstsx.write(1)
        self.video = True
        video = Video(self.display_video)
        video.signals.load.connect(self.show_images)
        self.threadpool.start(video)
        self.wquit.setEnabled(True)
        self.mode = 'video'
        self.wrestartvideo.setEnabled(True)

    def restart_video(self):
        self.stopex.write(1)
        time.sleep(3)
        self.cdsmode.write(1)
        self.readmode.write(3)
        self.go.write(1)
        # self.trkenapx.write(0)
        # self.trkfpspx.write('Passive')
        self.trkstop.write(1)
        time.sleep(1)
        # self.trkfpspx.write('1 second')
        # self.trkenapx.write(1)
        self.trkstsx.write(1)
        self.start_video()


    ##TODO remove this when switching back to video mode is replaced with restart_video
    def start_video(self):
        self.video = True
        print("video started")
        video = Video(self.display_video)
        video.signals.load.connect(self.show_images)
        left, right, up, down = self.getROI()
        self.roi_info.setText(f"ROI: {left} {right} {up} {down}")
        self.threadpool.start(video)
        self.mode = 'video'

    #TODO make this actually stop video mode
    def stop_video(self):
        print("video stopped")
        self.video = False

    def reboot_trick(self):
        print("Rebooting TRICK...")

    def display_video(self, file_callback):
        while self.video:
            file_callback.emit(self.roipixels)
            time.sleep(1)

    def show_images(self, pix):
        image = self.pixels_to_image(pix)
        self.img.load_data(image)
        self.fitsimage.set_image(self.img)
        self.MainWindow.resize(330, 370)

    def pixels_to_image(self, pix):
        lst = str(pix).strip().replace(':', '').split()
        dct = {int(lst[i]): float(lst[i + 1]) for i in range(0, len(lst), 2)}
        size = int(math.sqrt(len(dct)))

        image = []

        for i in range(size):
            row = []
            for k in range(size):
                row.append(dct[(i*16)+k])
            image.append(row)

        image = np.array(image)

        return(image)

    def full_frame_mode(self):
        self.wfullframemode.setVisible(False)
        self.wvideomode.setVisible(True)
        self.stop_video()
        self.fitsimage.clear()
        self.MainWindow.resize(700, 700)
        self.winittrick.setVisible(False)
        self.wrestartvideo.setVisible(False)
        self.wreboottrick.setVisible(False)
        self.wreboottrick.setVisible(False)
        self.display_hbox.setGeometry(QtCore.QRect(0, 0, 700, 512))
        self.fitsimage.get_widget().resize(512,512)
        self.roi_info.setGeometry(QtCore.QRect(20, 520, 171, 20))
        self.readout.setGeometry(QtCore.QRect(210, 520, 281, 16))
        self.wcolor.setVisible(True)
        self.wcut.setVisible(True)
        self.wsky.setVisible(True)
        self.wtakeff.setVisible(True)
        self.wsetroi.setVisible(True)
        self.image_info.setVisible(True)
        self.sky_info.setVisible(True)
        self.filt_info.setVisible(True)
        self.wfullframemode.setEnabled(False)
        self.wquit.setVisible(False)
        self.wvideomode.setEnabled(True)
        self.menubar.setVisible(True)
        self.start_scan()
        self.mode = 'fullframe'

    def video_mode(self):
        self.wvideomode.setVisible(False)
        self.wfullframemode.setVisible(True)
        self.stop_scan()
        self.fitsimage.clear()
        self.fitsimage.rotate(0)
        self.hbox_layout.setGeometry(QtCore.QRect(0, 0, 310, 240))
        self.fitsimage.get_widget().resize(240,240)
        self.readout.setGeometry(QtCore.QRect(30, 240, 251, 16))
        self.winittrick.setVisible(True)
        self.wrestartvideo.setVisible(True)
        self.wreboottrick.setVisible(True)
        self.wfullframemode.setEnabled(True)
        self.roi_info.setGeometry(QtCore.QRect(30, 20, 171, 20))
        self.wquit.setVisible(True)
        self.wcolor.setVisible(False)
        self.wcut.setVisible(False)
        self.wsky.setVisible(False)
        self.wtakeff.setVisible(False)
        self.wsetroi.setVisible(False)
        self.image_info.setVisible(False)
        self.sky_info.setVisible(False)
        self.filt_info.setVisible(False)
        self.wvideomode.setEnabled(False)
        self.menubar.setVisible(False)
        self.restart_video()

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
        print("scan started...")
        scanner = Scanner(self.scan)
        scanner.signals.load.connect(self.processData)
        self.threadpool.start(scanner)

    def stop_scan(self):
        self.scanning = False
        print('Scanning stopped.')

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
        self.stopex.write(1)
        time.sleep(3)
        self.timfile.write('/kroot/rel/default/data/nirtts_rel_v1.6.lod')
        # self.init.write(1)
        # self.sampmode.write(5)
        self.cdsmode.write(1)
        self.readmode.write(1)
        self.itime.write(0)
        self.getkw.write(1)
        self.getdcskw.write(1)
        self.getaokw.write(1)
        self.go.write(1)

    def set_roi(self):
        self.trickxpos.write(self.xclick-8)
        self.trickypos.write(self.yclick-8)
        print("TRICK ROI set")
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
            print(self.progress.read())
            time.sleep(wait_time)

    def getROI(self):
        left = int(self.trickxpos.read())  + 8 - int(self.trickxsize.read())*3
        right = int(self.trickxpos.read()) + 8 + int(self.trickxsize.read())*3
        up = int(self.trickypos.read()) + 8 - int(self.trickysize.read())*3
        down = int(self.trickypos.read()) + 8 + int(self.trickysize.read())*3
        print("ROI box: %d %d %d %d" %(left, right, up, down))
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
        rot = float(header['ROTPOSN'])
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
        try:
            self.fitsimage.get_canvas().get_object_by_tag(self.picktag)
            self.fitsimage.get_canvas().delete_object_by_tag(self.picktag)
            self.pickbox = self.recdc(self.xclick-30, self.yclick-30, self.xclick+30, self.yclick+30, color='red')
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
        # self.fitsimage.set_pan(data_x, data_y)
        self.xclick = data_x
        self.yclick = data_y
        ##todo video mode adjusting ROI
        if self.mode == "video":
            x = float(self.trickxpos.read()) - (float(self.trickxsize.read())/2.0 - self.xclick)
            y = float(self.trickypos.read()) - (float(self.trickysize.read())/2.0 - float(self.yclick))
            self.trickxpos.write(x)
            self.trickypos.write(y)
            # self.trkenapx.write(0)
            # self.trkfpspx.write('Passive')
            self.trkstop.write(1)
            time.sleep(1)
            # self.trkfpspx.write('1 second')
            # self.trkenapx.write(1)
            self.trkstsx.write(1)
        else:
            self.wsetroi.setEnabled(True)
            self.pickstar(self.fitsimage)
        # self.wsetroi.setEnabled(True)
        # self.pickstar(self.fitsimage)



def main():
    ##Write dummy file so walkDirectory caches it in the beginning

    app = QtGui.QApplication([])

    # ginga needs a logger.
    # If you don't want to log anything you can create a null logger by
    # using null=True in this call instead of log_stderr=True
    logger = log.get_logger("TrickManager", log_stderr=True, level=40)

    MainWindow = QtGui.QMainWindow()

    w = FitsViewer(logger, MainWindow)
    MainWindow.show()
    w.start_video()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()

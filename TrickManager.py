import time, threading
import math

import numpy as np

from ginga import Bindings
from ginga.misc import log
from ginga.qtw.QtHelp import QtGui, QtCore
from ginga.qtw.ImageViewQt import CanvasView
from ginga.util.loader import load_data
from ginga.AstroImage import AstroImage

import ktl

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

    def __init__(self, logger):
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

        self.img = AstroImage()

        self.threadpool = QtCore.QThreadPool()

        # create the ginga viewer and configure it
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

        w = fi.get_widget()
        w.resize(512, 512)

        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(QtCore.QMargins(2, 2, 2, 2))
        vbox.setSpacing(1)
        vbox.addWidget(w, stretch=1)

        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))

        self.roi_info = QtGui.QLabel("")

        hbox.addStretch(1)
        hbox.addWidget(self.roi_info, stretch = 0)

        text = f"ROI {self.trickxpos.read()} {self.trickypos.read()}"
        self.roi_info.setText(text)

        hw = QtGui.QWidget()
        hw.setLayout(hbox)
        vbox.addWidget(hw, stretch=0)

        hbox2 = QtGui.QHBoxLayout()
        hbox2.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))

        self.readout = QtGui.QLabel("")

        hbox2.addStretch(1)
        hbox2.addWidget(self.readout, stretch = 0)

        self.wcut = QtGui.QComboBox()
        for name in fi.get_autocut_methods():
            self.wcut.addItem(name)
        self.wcut.currentIndexChanged.connect(self.cut_change)
        self.wcolor = QtGui.QComboBox()
        for name in fi.get_color_algorithms():
            self.wcolor.addItem(name)
        self.wcolor.currentIndexChanged.connect(self.color_change)
        hbox2.addStretch(1)

        for w in (self.wcut, self.wcolor):
            hbox2.addWidget(w, stretch=0)

        hw2 = QtGui.QWidget()
        hw2.setLayout(hbox2)
        vbox.addWidget(hw2, stretch=0)

        hbox3 = QtGui.QHBoxLayout()
        hbox3.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))
        self.winittrick = QtGui.QPushButton("Init Trick")
        self.winittrick.clicked.connect(self.init_trick)
        self.wrestartvideo = QtGui.QPushButton("Restart Video")
        self.wrestartvideo.clicked.connect(self.restart_video)
        self.wreboottrick = QtGui.QPushButton("Reboot Trick")
        self.wreboottrick.clicked.connect(self.reboot_trick)
        fi.set_callback('cursor-changed', self.motion_cb)
        hbox3.addStretch(1)
        for w in (self.winittrick, self.wrestartvideo, self.wreboottrick):
            hbox3.addWidget(w, stretch=0)

        hw3 = QtGui.QWidget()
        hw3.setLayout(hbox3)
        vbox.addWidget(hw3, stretch=0)

        hbox4 = QtGui.QHBoxLayout()
        hbox4.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))
        self.wfullframemode = QtGui.QPushButton("Full Frame Mode")
        self.wfullframemode.clicked.connect(self.full_frame_mode)
        self.wfullframemode = QtGui.QPushButton("Video Mode")
        self.wfullframemode.clicked.connect(self.video_mode)
        self.wfullframemode.setEnabled(False)
        self.wstopvideo = QtGui.QPushButton("Stop Video")
        self.wstopvideo.clicked.connect(self.stop_video)
        self.wstopvideo.setEnabled(False)
        self.wquit = QtGui.QPushButton("Quit")
        self.wquit.clicked.connect(self.quit)
        hbox4.addStretch(1)
        for w in (self.wfullframemod, self.wstopvideo, self.wstopvideo, self.wquit):
            hbox4.addWidget(w, stretch=0)

        hw4 = QtGui.QWidget()
        hw4.setLayout(hbox4)
        vbox.addWidget(hw4, stretch=0)

        vw = QtGui.QWidget()
        self.setCentralWidget(vw)
        vw.setLayout(vbox)
        self.recdc, self.compdc = self.add_canvas()
        self.boxtag = "roi-box"
        self.picktag = "pick-box"


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
        time.sleep(0.5)
        self.threadpool = False
        self.deleteLater()

    ##TODO verify init trick, restart video, stop video, reboot trick

    def init_trick(self):
        print("Initing TRICK")
        self.wstopvideo.setEnabled(True)
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
        self.trkenapx.write(0)
        self.trkfpspx.write('Passive')
        self.trkstop.write(1)
        time.sleep(1)
        self.trkfpspx.write('1 second')
        self.trkenapx.write(1)
        self.trkstsx.write(1)
        self.video = True
        video = Video(self.display_video)
        video.signals.load.connect(self.show_images)
        self.threadpool.start(video)

    def restart_video(self):
        self.stopex.write(1)
        time.sleep(3)
        self.cdsmode.write(1)
        self.readmode.write(3)
        self.go.write(1)
        time.sleep(3)
        self.trkenapx.write(0)
        self.trkfpspx.write('Passive')
        self.trkstop.write(1)
        time.sleep(1)
        self.trkfpspx.write('1 second')
        self.trkenapx.write(1)
        self.trkstsx.write(1)

    def stop_video(self):
        print("video stopped")
        self.wstopvideo.setEnabled(False)
        self.wquit.setEnabled(True)
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
        fi = CanvasView(self.logger, render='widget')
        w = fi.get_widget()
        w.resize(700, 700)
        self.wstopvideo.setEnabled(False)
        self.winittrick.setEnabled(False)
        self.wrestartvideo.setEnabled(False)
        self.wreboottrick.setEnabled(False)
        self.wreboottrick.setEnabled(False)
        self.wvideomode.setEnabled(True)

    def video_mode(self):
        fi = CanvasView(self.logger, render='widget')
        w = fi.get_widget()
        w.resize(500, 500)
        self.wstopvideo.setEnabled(True)
        self.winittrick.setEnabled(True)
        self.wrestartvideo.setEnabled(True)
        self.wreboottrick.setEnabled(True)
        self.wreboottrick.setEnabled(True)
        self.wvideomode.setEnabled(False)



def main():
    ##Write dummy file so walkDirectory caches it in the beginning

    app = QtGui.QApplication([])

    # ginga needs a logger.
    # If you don't want to log anything you can create a null logger by
    # using null=True in this call instead of log_stderr=True
    logger = log.get_logger("example1", log_stderr=True, level=40)



    w = FitsViewer(logger)
    w.resize(500, 500)
    w.show()
    app.setActiveWindow(w)
    w.raise_()
    w.activateWindow()
    app.exec_()

if __name__ == "__main__":
    main()

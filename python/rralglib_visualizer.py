import sys, os, time
from PyQt6.QtCore import Qt, QRect, QTimer, pyqtSignal
from PyQt6.QtGui import QAction, QKeyEvent, QDoubleValidator, QIntValidator, QRegularExpressionValidator
from PyQt6.QtWidgets import QApplication, QHBoxLayout, QLabel, QSlider, QVBoxLayout, QWidget, QLineEdit, QPushButton, QComboBox, QCheckBox, QMessageBox, QToolBar, QMenuBar, QFileDialog, QDialog, QTextEdit, QMenu, QGridLayout, QFrame
import pyqtgraph as pg
import numpy as np
import pandas as pd
from scipy.signal import butter, ellip, cheby1, bessel
from rralglib import srmac, cwt_peaks, terma, sos_filt, sqi_full, z_norm, find_peaks, cwt_peaks_oa

### Most important default data parameters
sample_rate = 64
### Other parameters can be changed as needed within MainWindow.__init__

### Helper function to read a .csv file
def read_file_csv_like(filepath, key_data=None, key_time=None, delimiter=None):
    """
    Generic function to read a .csv file with predefined keys
    """
    data_file = pd.read_csv(filepath, delimiter=delimiter, encoding="utf-8")

    if key_data is None:
        key_data = 1

    if key_time is None:
        key_time = 0

    key_time = data_file.keys()[key_time]
    key_data = data_file.keys()[key_data]

    time_raw = data_file[key_time]
    data_raw = data_file[key_data]

    data_raw = np.array(data_raw.tolist())
    time_raw = np.array(time_raw.tolist())

    return data_raw, time_raw
    
### Helper function to read a .csv file without using pandas
def read_file_csv_like_native(filepath, key_data=None, key_time=None, delimiter=None, start_index=0):
    """
    Generic function to read a .csv file with predefined keys implemented using native python file processing
    """

    if delimiter is None:
        delimiter = ","

    if key_data is None:
        key_data = 1

    if key_time is None:
        key_time = 0

    time_raw = []
    data_raw = []

    with open(filepath, "r", encoding="utf-8") as file:
        for i, line in enumerate(file):
            if i < start_index:
                continue
            try:
                line = line.split(delimiter)
                try:
                    time_raw.append(float(line[key_time]))
                except: 
                    ### Neulog style time column
                    time = line[key_time].split("'")[1].split(":")
                    time = float(time[0]) * 3600.0 + float(time[1]) * 60.0 + float(time[2])
                    time_raw.append(time)
                try:
                    data_raw.append(float(line[key_data]))
                except:
                    data_raw.append(float(line[key_data].split("\n")[0]))
            except:
                print("Couldn't parse line "+str(i)+"... "+str(line))

    data_raw = np.array(data_raw)
    time_raw = np.array(time_raw)

    return data_raw, time_raw

### GUI ###################################################

### Error dialog window
class ErrorDialog(QDialog):
    def __init__(self, parent, message):
        super(ErrorDialog, self).__init__(parent=parent)

        self.setWindowTitle("ERROR")
        self.main_layout = QVBoxLayout(self)

        self.setMinimumWidth(200)
        self.setMinimumHeight(100)

        self.msg = QLabel(str(message))
        self.main_layout.addWidget(self.msg)

        ### Confirm button
        self.button_confirm = QPushButton("OK")
        self.button_confirm.pressed.connect(self.ok_button_pressed)
        self.main_layout.addWidget(self.button_confirm)

    def ok_button_pressed(self):
        self.accept()

### About window
class AboutDialog(QDialog):
    def __init__(self, parent):
        super(AboutDialog, self).__init__(parent=parent)

        self.setWindowTitle("About")
        self.main_layout = QVBoxLayout(self)

        self.setFixedWidth(300)
        self.setFixedHeight(200)

        self.header = QLabel("About")
        self.header.setStyleSheet("font-weight:bold;")
        self.main_layout.addWidget(self.header)

        self.paragraph = QLabel("The rralglib visualizer is a PyQt6 based Python tool for signal playback and peak detection using some of the algorithms implemented in the Python module of the rralglib library. This tool was created as part of the Qualification Work for the University of Latvia.")
        self.paragraph.setWordWrap(True)
        self.main_layout.addWidget(self.paragraph)
        
        self.end = QLabel("Author: Toms Racinskis")
        self.end.setStyleSheet("font-weight:italic; color: #999999;")
        self.main_layout.addWidget(self.end)

    def ok_button_pressed(self):
        self.accept()

### Read file dialog window
class ReadFileDialog(QDialog):
    def __init__(self, parent, selected_file):
        super(ReadFileDialog, self).__init__(parent=parent)

        self.preview_length = 32
        self.file_preview = ""

        self.selected_file = selected_file
        self.xdata = None
        self.ydata = None
        self.sample_rate = None
        self.filename = None
        self.delimiter = None
        self.key_time = None
        self.key_data = None

        if self.selected_file is None:
            ErrorDialog(self, "Error: No valid file selected").exec()
            self.reject()
            self.close()

        self.setWindowTitle("Open file")
        self.main_layout = QVBoxLayout(self)

        self.setMinimumWidth(500)
        self.setMinimumHeight(450)

        ### Preview layout
        self.preview_layout = QVBoxLayout()
        self.main_layout.addLayout(self.preview_layout)

        ### Selected file name label
        self.selection_label = QLabel("Previewing selected file: "+str(self.selected_file))
        self.selection_label.setWordWrap(True)
        self.selection_label.setStyleSheet("font: bold 14px;")
        self.main_layout.addWidget(self.selection_label)

        ### Selected file preview
        self.selection_preview = QTextEdit()
        self.selection_preview.setDisabled(True)
        self.main_layout.addWidget(self.selection_preview)

        ### Input frame
        self.input_frame = QFrame()
        self.input_frame.setFrameShape(QFrame.Shape.StyledPanel)
        self.input_frame.setLineWidth(3)
        self.main_layout.addWidget(self.input_frame)

        ### Input layout
        self.input_layout = QGridLayout(self.input_frame)

        ### Sample rate input
        self.input_samplerate_label_1 = QLabel("Sample rate: ")
        self.input_samplerate_label_1.setFixedWidth(75)
        self.input_samplerate = QLineEdit()
        self.input_samplerate.setMaximumWidth(50)
        self.input_samplerate.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.input_samplerate.setValidator(QIntValidator(1,10000))
        self.input_samplerate.textChanged.connect(self.update_sample_rate)
        self.input_samplerate_label_2 = QLabel("(Leave blank to calculate sample rate automatically)")
        self.input_samplerate_label_2.setStyleSheet("color: grey; font: italic")
        self.input_layout.addWidget(self.input_samplerate_label_1, 0, 1)
        self.input_layout.addWidget(self.input_samplerate, 0, 2)
        self.input_layout.addWidget(self.input_samplerate_label_2, 0, 3)

        ### Delimiter input
        self.input_delimiter_label_1 = QLabel("Delimiter: ")
        self.input_delimiter_label_1.setFixedWidth(75)
        self.input_delimiter = QLineEdit()
        self.input_delimiter.setMaximumWidth(50)
        self.input_delimiter.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.input_delimiter.setMaxLength(1)
        self.input_delimiter.setText(",")
        self.input_delimiter.textChanged.connect(self.update_delimiter)
        self.input_delimiter_label_2 = QLabel("(Default: ,)")
        self.input_delimiter_label_2.setStyleSheet("color: grey; font: italic")
        self.input_layout.addWidget(self.input_delimiter_label_1, 1, 1)
        self.input_layout.addWidget(self.input_delimiter, 1, 2)
        self.input_layout.addWidget(self.input_delimiter_label_2, 1, 3)

        ### Time column input
        self.input_timecolumn_label_1 = QLabel("Time column: ")
        self.input_timecolumn_label_1.setFixedWidth(75)
        self.input_timecolumn = QLineEdit()
        self.input_timecolumn.setMaximumWidth(50)
        self.input_timecolumn.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.input_timecolumn.setMaxLength(2)
        self.input_timecolumn.setText("0")
        self.input_timecolumn.setValidator(QIntValidator(1,10000))
        self.input_timecolumn.textChanged.connect(self.update_time_column)
        self.input_timecolumn_label_2 = QLabel("(Default: 0)")
        self.input_timecolumn_label_2.setStyleSheet("color: grey; font: italic")
        self.input_layout.addWidget(self.input_timecolumn_label_1, 2, 1)
        self.input_layout.addWidget(self.input_timecolumn, 2, 2)
        self.input_layout.addWidget(self.input_timecolumn_label_2, 2, 3)

        ### Data column input
        self.input_datacolumn_label_1 = QLabel("Data column: ")
        self.input_datacolumn_label_1.setFixedWidth(75)
        self.input_datacolumn = QLineEdit()
        self.input_datacolumn.setMaximumWidth(50)
        self.input_datacolumn.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.input_datacolumn.setMaxLength(2)
        self.input_datacolumn.setText("1")
        self.input_datacolumn.setValidator(QIntValidator(1,10000))
        self.input_datacolumn.textChanged.connect(self.update_data_column)
        self.input_datacolumn_label_2 = QLabel("(Default: 1)")
        self.input_datacolumn_label_2.setStyleSheet("color: grey; font: italic")
        self.input_layout.addWidget(self.input_datacolumn_label_1, 3, 1)
        self.input_layout.addWidget(self.input_datacolumn, 3, 2)
        self.input_layout.addWidget(self.input_datacolumn_label_2, 3, 3)

        ## Bottom layout
        self.bottom_layout = QHBoxLayout()
        self.main_layout.addLayout(self.bottom_layout)

        ### Cancel button
        self.button_cancel = QPushButton("Cancel")
        self.button_cancel.pressed.connect(self.cancel_button_pressed)
        self.bottom_layout.addWidget(self.button_cancel)

        ### Confirm button
        self.button_confirm = QPushButton("Confirm")
        self.button_confirm.pressed.connect(self.confirm_button_pressed)
        self.bottom_layout.addWidget(self.button_confirm)

        ### Set file preview
        self.set_file_preview()

    ### Try to update the sample rate
    def update_sample_rate(self):
        try:
            self.sample_rate = int(self.input_samplerate.text())
        except:
            self.sample_rate = None
            ErrorDialog(self, "Invalid sample rate.").exec()

    ### Try to update the delimiter
    def update_delimiter(self):
        try:
            self.delimiter = str(self.input_delimiter.text())
            if len(self.delimiter) != 1:
                self.delimiter = None
                ErrorDialog(self, "Delimiter must be exactly one character.").exec()
        except:
            self.delimiter = None
            ErrorDialog(self, "Invalid delimiter.").exec()

    ### Try to update the time column key
    def update_time_column(self):
        try:
            self.key_time = int(self.input_timecolumn.text())
        except:
            self.key_time = None
            ErrorDialog(self, "Invalid time column.").exec()

    ### Try to update the data column key
    def update_data_column(self):
        try:
            self.key_data = int(self.input_datacolumn.text())
        except:
            self.key_data = None
            ErrorDialog(self, "Invalid data column.").exec()

    ### Read the first few lines of the file to show a file preview
    def set_file_preview(self):
        try:
            with open(self.selected_file, encoding="utf-8") as file:
                for i, line in enumerate(file):
                    if i > self.preview_length:
                        self.file_preview += "..."
                        break
                    self.file_preview += line
            self.selection_preview.setText(self.file_preview)
        except Exception as e:
            self.selection_preview.setText("Could not open the selected file. ["+str(e)+"]")

    ### Cancel this dialog
    def cancel_button_pressed(self):
        self.reject()
        self.close()

    ### Read the full file and confirm this dialog
    def confirm_button_pressed(self):
        key_data = None
        key_time = None
        delimiter = None

        if self.key_data != None:
            key_data = self.key_data
        else:
            key_data = 0

        if self.key_time != None:
            key_time = self.key_time
        else:
            key_data = 1

        if self.delimiter != None:
            delimiter = self.delimiter    
        else:
            delimiter = ","

        try:
            self.ydata, self.xdata = read_file_csv_like_native(self.selected_file, key_data=key_data, key_time=key_time, delimiter=delimiter)
            self.filename = self.selected_file

            if len(self.input_samplerate.text()) > 0:
                self.sample_rate = int(self.input_samplerate.text())

            if self.sample_rate is None:
                if len(self.xdata) > 1:
                    self.sample_rate = int(1.0 / self.xdata[1]-self.xdata[0])

            self.accept()
            self.close()
        except Exception as e:
            ErrorDialog(self, "Could not open the selected file. ["+str(e)+"]").exec()
            self.reject()
            self.close()

### The main program window
class MainWindow(QWidget):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent=parent)
        # self.setWindowFlag(Qt.WindowType.FramelessWindowHint, True) ### Remove window frame
        self.main_layout = QVBoxLayout(self)
        self.setWindowTitle("rralglib signal visualizer")

        ### Main loop stuff
        self.file_loaded = False
        self.running = False
        self.analyze = False
        self.tick = 0

        ### Framerate stuff
        self.framerate = 240  ### Approximate framerate at which to run the main application loop ### Set to zero to update as fast as possible
        self.seconds_per_second = 5 ### The seconds of data to process per second of playback
        self.samples_per_second = 512 ### The amount of samples per second to increment in the live plot ### Make this a slider
        self.algorithm_delta_ticks = 1

        ### Raw data storage
        self.data_raw = []
        self.time_raw = []
        self.sample_rate = sample_rate

        ### Data window variables
        self.data_cutoff = 0
        self.data_window_size = 40
        self.data_index = 0
        self.max_data_index = 0

        ### Selected algorithm
        self.algorithms = ["srmac","cwt","terma","find_peaks"]
        self.algorithm = self.algorithms[0]
        self.coef_resolution = 10000

        ### Default filter cutoffs
        self.lpf_cutoff = 0.8
        self.hpf_cutoff = 0.05

        ### GUI
        #region

        ### Stylesheets
        self.stylesheet_lineedit_valid = "background-color: #FFFFFF;"
        self.stylesheet_lineedit_invalid = "background-color: #FF0000;"

        # Top bar layout
        self.top_layout = QHBoxLayout()
        self.main_layout.addLayout(self.top_layout)

        # Mid view layout
        self.mid_layout = QVBoxLayout()
        self.main_layout.addLayout(self.mid_layout)

        # Bottom bar layout
        self.bot_layout = QVBoxLayout()
        self.main_layout.addLayout(self.bot_layout)

        ### Menu bar (file, edit etc)
        self.menu = QMenuBar()
        self.menu_file = self.menu.addMenu("&File")
        self.menu_edit = self.menu.addMenu("&Edit")
        self.menu_view = self.menu.addMenu("&View")
        self.menu_help = self.menu.addMenu("&Help")
        self.top_layout.addWidget(self.menu)

        ### File menu actions
        self.action_open_file = QAction("Open file (O)")
        self.action_open_file.setStatusTip("Open file (O)")
        self.action_open_file.triggered.connect(self.button_open_new_file)
        self.menu_file.addAction(self.action_open_file)

        ### Edit menu actions
        self.action_start = QAction("Start/stop playback (Q)")
        self.action_start.setStatusTip("Start/stop playback (Q)")
        self.action_start.triggered.connect(self.play_button_pressed)
        self.menu_edit.addAction(self.action_start)

        self.action_reset = QAction("Reset plot (R)")
        self.action_reset.setStatusTip("Reset plot (R)")
        self.action_reset.triggered.connect(self.reset_button_pressed)
        self.menu_edit.addAction(self.action_reset)

        ### View menu actions
        self.action_show_algs = QAction("Show algorithm bar")
        self.action_show_algs.setStatusTip("Show algorithm bar")
        self.action_show_algs.setCheckable(True)
        self.action_show_algs.setChecked(False)
        self.action_show_algs.triggered.connect(self.toggle_algorithm_bar)
        self.menu_view.addAction(self.action_show_algs)

        self.action_show_midline = QAction("Show midline")
        self.action_show_midline.setStatusTip("Show midline")
        self.action_show_midline.setCheckable(True)
        self.action_show_midline.setChecked(False)
        self.action_show_midline.triggered.connect(self.toggle_midline)
        self.menu_view.addAction(self.action_show_midline)

        self.action_show_peaks = QAction("Show detected peaks")
        self.action_show_peaks.setStatusTip("Show detected peaks")
        self.action_show_peaks.setCheckable(True)
        self.action_show_peaks.setChecked(True)
        self.action_show_peaks.triggered.connect(self.toggle_peaks)
        self.menu_view.addAction(self.action_show_peaks)

        self.action_show_signal = QAction("Show signal")
        self.action_show_signal.setStatusTip("Show signal")
        self.action_show_signal.setCheckable(True)
        self.action_show_signal.setChecked(True)
        self.action_show_signal.triggered.connect(self.toggle_signal)
        self.menu_view.addAction(self.action_show_signal)

        ### Help menu actions
        self.action_about = QAction("About")
        self.action_about.setStatusTip("About")
        self.action_about.triggered.connect(self.show_about)
        self.menu_help.addAction(self.action_about)

        # Main plot widget
        self.plot_main = pg.PlotWidget()
        self.plot_main.setMinimumHeight(200)
        self.plot_main_legend = self.plot_main.plotItem.addLegend()
        self.mid_layout.addWidget(self.plot_main)
        self.plot_main.plotItem.getAxis(name='left').setLabel(text='Signal amplitude')
        self.plot_main.plotItem.getAxis(name='bottom').setLabel(text='Time (s)')

        # Data plot
        self.plot_signal = self.plot_main.plot(pen=pg.mkPen(width=1, color="#922"), name="Signal data")

        # Peak scatter
        self.peak_scatter = pg.ScatterPlotItem(pen=pg.mkPen(width=3, color="#1F1"), symbol='o', size=2, name="Detected peaks")
        self.plot_main.addItem(self.peak_scatter)

        # Plot midline
        self.midline = pg.InfiniteLine(angle=90, movable=False)
        self.midline.setVisible(False)
        self.plot_main.addItem(self.midline)

        # Progress bar
        self.progress_bar = QSlider()
        self.progress_bar.setOrientation(Qt.Orientation.Horizontal)
        progress_bar_stylesheet = "QSlider::handle{border-radius: 5px; border: 1px solid #999999}"
        progress_bar_stylesheet += "QSlider::handle:hover{border: 2px solid #999999}"
        self.progress_bar.setStyleSheet(progress_bar_stylesheet)
        # self.progress_bar.mousePressEvent.connect(self.start_dragging_progress_bar)
        self.progress_bar.valueChanged.connect(self.dragging_progress_bar)
        # self.progress_bar..connect(self.stop_dragging_progress_bar)
        self.mid_layout.addWidget(self.progress_bar)

        # First row
        self.row_1 = QHBoxLayout()
        self.bot_layout.addLayout(self.row_1)

        # Time indicator label
        self.time_indicator_samples = QLabel()
        self.time_indicator_samples.setText("0 / 0")
        self.row_1.addWidget(self.time_indicator_samples)

        # Start button
        self.button_play = QPushButton("Start (Q)")
        self.button_play.pressed.connect(self.play_button_pressed)
        self.row_1.addWidget(self.button_play)

        # Reset button
        self.button_reset = QPushButton("Reset (R)")
        self.button_reset.pressed.connect(self.reset_button_pressed)
        self.row_1.addWidget(self.button_reset)

        # Algorithm layout
        self.algorithm_layout = QHBoxLayout()

        # Algorithm container
        self.algorithm_container = QWidget()
        self.algorithm_container.setLayout(self.algorithm_layout)
        self.algorithm_container.hide()
        self.bot_layout.addWidget(self.algorithm_container)

        # Algorithm selector
        self.algselector = QComboBox()
        self.algselector.setMaximumWidth(200)
        self.algselector.setMinimumWidth(200)
        self.algselector.addItems(self.algorithms)
        self.algselector.setCurrentIndex(0)
        self.algselector.currentIndexChanged.connect(self.change_algorithm)
        self.algorithm_layout.addWidget(self.algselector)

        # Analyze checkbox
        self.analyze_check = QCheckBox(text="Analyze")
        self.analyze_check.setChecked(False)
        self.analyze_check.checkStateChanged.connect(self.analyze_check_toggled)
        self.algorithm_layout.addWidget(self.analyze_check)

        # Set to defaults
        self.button_default = QPushButton(text="Default")
        self.button_default.pressed.connect(self.default_button_pressed)
        self.algorithm_layout.addWidget(self.button_default)

        # SRMAC parameter container
        self.srmac_parameter_layout = QGridLayout()
        self.srmac_parameter_container = QWidget()
        self.srmac_parameter_container.setLayout(self.srmac_parameter_layout)
        self.algorithm_layout.addWidget(self.srmac_parameter_container)
        self.srmac_parameter_container.hide()

        self.params_srmac_coef_fast = None
        self.params_srmac_coef_slow = None
        self.params_srmac_coef_cross = None
        self.params_srmac_th = None
        self.params_srmac_width = None

        # coef_fast, coef_slow, coef_cross, width, threshold
        self.srmac_coef_fast_line_edit = QLineEdit()
        self.srmac_coef_fast_line_edit_label = QLabel("coef_fast")
        self.srmac_coef_fast_line_edit.setText("0.0")
        self.srmac_coef_fast_line_edit.setValidator(QDoubleValidator(0.0, 1.0, 5, notation=QDoubleValidator.Notation.StandardNotation))
        self.srmac_coef_fast_line_edit.textChanged.connect(lambda: self.update_srmac_parameters("params_srmac_coef_fast"))
        self.srmac_parameter_layout.addWidget(self.srmac_coef_fast_line_edit, 0, 0)
        self.srmac_parameter_layout.addWidget(self.srmac_coef_fast_line_edit_label, 0, 1)
        
        self.srmac_coef_slow_line_edit = QLineEdit()
        self.srmac_coef_slow_line_edit_label = QLabel("coef_slow")
        self.srmac_coef_slow_line_edit.setText("0.0")
        self.srmac_coef_slow_line_edit.setValidator(QDoubleValidator(0.0, 1.0, 5, notation=QDoubleValidator.Notation.StandardNotation))
        self.srmac_coef_slow_line_edit.textChanged.connect(lambda: self.update_srmac_parameters("params_srmac_coef_slow"))
        self.srmac_parameter_layout.addWidget(self.srmac_coef_slow_line_edit, 1, 0)
        self.srmac_parameter_layout.addWidget(self.srmac_coef_slow_line_edit_label, 1, 1)

        self.srmac_coef_cross_line_edit = QLineEdit()
        self.srmac_coef_cross_line_edit_label = QLabel("coef_cross")
        self.srmac_coef_cross_line_edit.setText("0.0")
        self.srmac_coef_cross_line_edit.setValidator(QDoubleValidator(0.0, 1.0, 5, notation=QDoubleValidator.Notation.StandardNotation))
        self.srmac_coef_cross_line_edit.textChanged.connect(lambda: self.update_srmac_parameters("params_srmac_coef_cross"))
        self.srmac_parameter_layout.addWidget(self.srmac_coef_cross_line_edit, 2, 0)
        self.srmac_parameter_layout.addWidget(self.srmac_coef_cross_line_edit_label, 2, 1)

        self.srmac_th_line_edit = QLineEdit()
        self.srmac_th_line_edit_label = QLabel("threshold")
        self.srmac_th_line_edit.setText("0.0")
        self.srmac_th_line_edit.setValidator(QDoubleValidator(0.0, 100000000000.0, 20, notation=QDoubleValidator.Notation.StandardNotation))
        self.srmac_th_line_edit.textChanged.connect(lambda: self.update_srmac_parameters("params_srmac_coef_th"))
        self.srmac_parameter_layout.addWidget(self.srmac_th_line_edit, 3, 0)
        self.srmac_parameter_layout.addWidget(self.srmac_th_line_edit_label, 3, 1) 

        self.srmac_width_line_edit = QLineEdit()
        self.srmac_width_line_edit_label = QLabel("width")
        self.srmac_width_line_edit.setText("0.5")
        self.srmac_width_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.srmac_width_line_edit.textChanged.connect(lambda: self.update_srmac_parameters("params_srmac_coef_width"))
        self.srmac_parameter_layout.addWidget(self.srmac_width_line_edit, 4, 0)
        self.srmac_parameter_layout.addWidget(self.srmac_width_line_edit_label, 4, 1)

        # TERMA parameter container
        self.terma_parameter_layout = QGridLayout()
        self.terma_parameter_container = QWidget()
        self.terma_parameter_container.setLayout(self.terma_parameter_layout)
        self.algorithm_layout.addWidget(self.terma_parameter_container)
        self.terma_parameter_container.hide()

        self.params_terma_w1 = None
        self.params_terma_w2 = None
        self.params_terma_width = None
        self.params_terma_b = None

        # w_1, w_2, b, width
        self.terma_w1_line_edit = QLineEdit()
        self.terma_w1_line_edit_label = QLabel("window_event")
        self.terma_w1_line_edit.setText("0.0")
        self.terma_w1_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.terma_w1_line_edit.textChanged.connect(lambda: self.update_terma_parameters("params_terma_w1"))
        self.terma_parameter_layout.addWidget(self.terma_w1_line_edit, 0, 0)
        self.terma_parameter_layout.addWidget(self.terma_w1_line_edit_label, 0, 1)

        self.terma_w2_line_edit = QLineEdit()
        self.terma_w2_line_edit_label = QLabel("window_cycle")
        self.terma_w2_line_edit.setText("0.0")
        self.terma_w2_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.terma_w2_line_edit.textChanged.connect(lambda: self.update_terma_parameters("params_terma_w2"))
        self.terma_parameter_layout.addWidget(self.terma_w2_line_edit, 1, 0)
        self.terma_parameter_layout.addWidget(self.terma_w2_line_edit_label, 1, 1)

        self.terma_b_line_edit = QLineEdit()
        self.terma_b_line_edit_label = QLabel("b_coef")
        self.terma_b_line_edit.setText("0.0")
        self.terma_b_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.terma_b_line_edit.textChanged.connect(lambda: self.update_terma_parameters("params_terma_b"))
        self.terma_parameter_layout.addWidget(self.terma_b_line_edit, 2, 0)
        self.terma_parameter_layout.addWidget(self.terma_b_line_edit_label, 2, 1)

        self.terma_width_line_edit = QLineEdit()
        self.terma_width_line_edit_label = QLabel("width")
        self.terma_width_line_edit.setText("0.0")
        self.terma_width_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.terma_width_line_edit.textChanged.connect(lambda: self.update_terma_parameters("params_terma_width"))
        self.terma_parameter_layout.addWidget(self.terma_width_line_edit, 3, 0)
        self.terma_parameter_layout.addWidget(self.terma_width_line_edit_label, 3, 1)

        # find_peaks parameter container
        self.findpeaks_parameter_layout = QGridLayout()
        self.findpeaks_parameter_container = QWidget()
        self.findpeaks_parameter_container.setLayout(self.findpeaks_parameter_layout)
        self.algorithm_layout.addWidget(self.findpeaks_parameter_container)
        self.findpeaks_parameter_container.hide()

        self.params_findpeaks_prominence = None
        self.params_findpeaks_heval = None
        self.params_findpeaks_width = None
        self.params_findpeaks_proximity = None

        # prominence, heval_ratio, width, proximity
        self.findpeaks_prominence_line_edit = QLineEdit()
        self.findpeaks_prominence_line_edit_label = QLabel("prominence")
        self.findpeaks_prominence_line_edit.setText("0.0")
        self.findpeaks_prominence_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.findpeaks_prominence_line_edit.textChanged.connect(lambda: self.update_findpeaks_parameters("params_findpeaks_prominence"))
        self.findpeaks_parameter_layout.addWidget(self.findpeaks_prominence_line_edit, 0, 0)
        self.findpeaks_parameter_layout.addWidget(self.findpeaks_prominence_line_edit_label, 0, 1)

        self.findpeaks_heval_line_edit = QLineEdit()
        self.findpeaks_heval_line_edit_label = QLabel("heval_ratio")
        self.findpeaks_heval_line_edit.setText("0.0")
        self.findpeaks_heval_line_edit.setValidator(QDoubleValidator(0.0, 1.0, 5, notation=QDoubleValidator.Notation.StandardNotation))
        self.findpeaks_heval_line_edit.textChanged.connect(lambda: self.update_findpeaks_parameters("params_findpeaks_heval"))
        self.findpeaks_parameter_layout.addWidget(self.findpeaks_heval_line_edit, 1, 0)
        self.findpeaks_parameter_layout.addWidget(self.findpeaks_heval_line_edit_label, 1, 1)

        self.findpeaks_width_line_edit = QLineEdit()
        self.findpeaks_width_line_edit_label = QLabel("width")
        self.findpeaks_width_line_edit.setText("0.0")
        self.findpeaks_width_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.findpeaks_width_line_edit.textChanged.connect(lambda: self.update_findpeaks_parameters("params_findpeaks_width"))
        self.findpeaks_parameter_layout.addWidget(self.findpeaks_width_line_edit, 2, 0)
        self.findpeaks_parameter_layout.addWidget(self.findpeaks_width_line_edit_label, 2, 1)

        self.findpeaks_proxim_line_edit = QLineEdit()
        self.findpeaks_proxim_line_edit_label = QLabel("proximity")
        self.findpeaks_proxim_line_edit.setText("0.0")
        self.findpeaks_proxim_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.findpeaks_proxim_line_edit.textChanged.connect(lambda: self.update_findpeaks_parameters("params_findpeaks_proximity"))
        self.findpeaks_parameter_layout.addWidget(self.findpeaks_proxim_line_edit, 3, 0)
        self.findpeaks_parameter_layout.addWidget(self.findpeaks_proxim_line_edit_label, 3, 1)

        # CWToa parameter container
        self.cwtoa_parameter_layout = QGridLayout()
        self.cwtoa_parameter_container = QWidget()
        self.cwtoa_parameter_container.setLayout(self.cwtoa_parameter_layout)
        self.algorithm_layout.addWidget(self.cwtoa_parameter_container)
        self.cwtoa_parameter_container.hide()

        self.params_cwtoa_f_min = None
        self.params_cwtoa_f_max = None
        self.params_cwtoa_scales = None
        self.params_cwtoa_width = None
        self.params_cwtoa_th = None
        self.params_cwtoa_margin = None

        # f_min, f_max, scales, width, threshold
        self.cwtoa_fmin_line_edit = QLineEdit()
        self.cwtoa_fmin_line_edit_label = QLabel("min_freq")
        self.cwtoa_fmin_line_edit.setText("0.0")
        self.cwtoa_fmin_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.cwtoa_fmin_line_edit.textChanged.connect(lambda: self.update_cwtoa_parameters("params_cwt_fmin"))
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_fmin_line_edit, 0, 0)
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_fmin_line_edit_label, 0, 1)

        self.cwtoa_fmax_line_edit = QLineEdit()
        self.cwtoa_fmax_line_edit_label = QLabel("max_freq")
        self.cwtoa_fmax_line_edit.setText("0.0")
        self.cwtoa_fmax_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.cwtoa_fmax_line_edit.textChanged.connect(lambda: self.update_cwtoa_parameters("params_cwt_fmax"))
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_fmax_line_edit, 1, 0)
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_fmax_line_edit_label, 1, 1)

        self.cwtoa_scales_line_edit = QLineEdit()
        self.cwtoa_scales_line_edit_label = QLabel("resolution")
        self.cwtoa_scales_line_edit.setText("0")
        self.cwtoa_scales_line_edit.setValidator(QIntValidator(1, 200))
        self.cwtoa_scales_line_edit.textChanged.connect(lambda: self.update_cwtoa_parameters("params_cwt_scales"))
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_scales_line_edit, 2, 0)
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_scales_line_edit_label, 2, 1)

        self.cwtoa_width_line_edit = QLineEdit()
        self.cwtoa_width_line_edit_label = QLabel("width")
        self.cwtoa_width_line_edit.setText("0.0")
        self.cwtoa_width_line_edit.setValidator(QDoubleValidator(0.0, 10000.0, 10, notation=QDoubleValidator.Notation.StandardNotation))
        self.cwtoa_width_line_edit.textChanged.connect(lambda: self.update_cwtoa_parameters("params_cwt_width"))
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_width_line_edit, 3, 0)
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_width_line_edit_label, 3, 1)

        self.cwtoa_th_line_edit = QLineEdit()
        self.cwtoa_th_line_edit_label = QLabel("threshold")
        self.cwtoa_th_line_edit.setText("0.0")
        self.cwtoa_th_line_edit.setValidator(QDoubleValidator(0.0, 100000000000.0, 20, notation=QDoubleValidator.Notation.StandardNotation))
        self.cwtoa_th_line_edit.textChanged.connect(lambda: self.update_cwtoa_parameters("params_cwt_th"))
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_th_line_edit, 4, 0)
        self.cwtoa_parameter_layout.addWidget(self.cwtoa_th_line_edit_label, 4, 1)

        #endregion

        ### Main loop timer
        self.loop_timer = QTimer()
        self.loop_timer.timeout.connect(self.main_loop)
        if self.framerate == 0:
            self.loop_timer.start(0)
        else:
            self.loop_timer.start(int(1000/self.framerate))

        ### Sliding data timer
        self.data_timer = QTimer()
        self.data_timer.timeout.connect(self.increment_index)
        if self.samples_per_second == 0:
            self.data_timer.start(0)
        else:
            self.data_timer.start(int(1000/self.samples_per_second))

        ### Disable all widgets which require data to be loaded
        self.button_play.setEnabled(False)
        self.button_reset.setEnabled(False)
        self.menu_edit.setEnabled(False)
        self.progress_bar.setEnabled(False)
        self.algorithm_container.setEnabled(False)

        ### Set all algorithm parameters to default
        self.algorithm = self.algorithms[0]
        self.default_button_pressed()
        self.algorithm = self.algorithms[1]
        self.default_button_pressed()
        self.algorithm = self.algorithms[2]
        self.default_button_pressed()
        self.algorithm = self.algorithms[3]
        self.default_button_pressed()
        self.algorithm = self.algorithms[0]

    ### Reset data from the raw arrays and update anything that changes based on data parameters
    def reset_data(self):
        ### Make copies of the raw arrays
        self.time = np.array(self.time_raw)
        self.data = np.array(self.data_raw)

        ### Apply a bandpass filter to the data
        lpf = butter(3, self.lpf_cutoff / (self.sample_rate * 0.5), analog=False, btype="low", output="sos")
        hpf = butter(3, self.hpf_cutoff / (self.sample_rate * 0.5), analog=False, btype="high", output="sos")
        # lpf = ellip(3, 0.1, 20, self.lpf_cutoff / (self.sample_rate * 0.5), analog=False, btype="low", output="sos")
        # hpf = ellip(3, 0.1, 20, self.hpf_cutoff / (self.sample_rate * 0.5), analog=False, btype="high", output="sos")
        # lpf = cheby1(3, 0.1, self.lpf_cutoff / (self.sample_rate * 0.5), analog=False, btype="low", output="sos")
        # hpf = cheby1(3, 0.1, self.hpf_cutoff / (self.sample_rate * 0.5), analog=False, btype="high", output="sos")
        self.data -= np.mean(self.data) # try to mitigate the initial filter transient
        self.data = sos_filt(self.data, hpf)
        self.data = sos_filt(self.data, lpf)

        ### Apply cutoff
        self.data = self.data[self.data_cutoff:]
        self.time = self.time[self.data_cutoff:]

        ### Update the progress bar range
        self.max_data_index = len(self.data)-self.data_window_size*self.sample_rate
        self.progress_bar.setRange(0, self.max_data_index)

        ### Update the data timer interval
        self.data_timer.stop()
        self.samples_per_second = int(self.sample_rate * self.seconds_per_second)
        self.algorithm_delta_ticks = self.sample_rate
        if self.samples_per_second == 0:
            self.data_timer.start(0)
        else:
            self.data_timer.start(int(1000/self.samples_per_second))

        ### Update the plot once
        self.update_plot()

    def update_srmac_parameters(self, param):
        if param == "params_srmac_coef_fast":
            ### EWMA fast coefficient
            try:
                if self.srmac_coef_fast_line_edit.text() is not None:
                    self.params_srmac_coef_fast = float(self.srmac_coef_fast_line_edit.text())
                    self.srmac_coef_fast_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: SRMAC fast coefficient could not be converted to float")
                self.params_srmac_coef_fast = 0.0508
                self.srmac_coef_fast_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_srmac_coef_slow":
            ### EWMA slow coefficient
            try:
                if self.srmac_coef_slow_line_edit.text() is not None:
                    self.params_srmac_coef_slow = float(self.srmac_coef_slow_line_edit.text())
                    self.srmac_coef_slow_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: SRMAC slow coefficient could not be converted to float")
                self.params_srmac_coef_slow = 0.0103
                self.srmac_coef_slow_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_srmac_coef_cross":
            ### EWMA cross coefficient
            try:
                if self.srmac_coef_cross_line_edit.text() is not None:
                    self.params_srmac_coef_cross = float(self.srmac_coef_cross_line_edit.text())
                    self.srmac_coef_cross_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: SRMAC cross coefficient could not be converted to float")
                self.params_srmac_coef_cross = 0.91
                self.srmac_coef_cross_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_srmac_coef_width":
            ### SRMAC width
            try:
                if self.srmac_width_line_edit.text() is not None:
                    self.params_srmac_width = float(self.srmac_width_line_edit.text())
                    self.srmac_width_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: SRMAC width could not be converted to Float")
                self.params_srmac_width = 0.5
                self.srmac_width_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_srmac_coef_th":
            ### SRMAC threshold
            try:
                if self.srmac_th_line_edit.text() is not None:
                    self.params_srmac_th = float(self.srmac_th_line_edit.text())
                    self.srmac_th_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: SRMAC threshold could not be converted to Float")
                self.params_srmac_th = 0
                self.srmac_th_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)

    def update_terma_parameters(self, param):
        if param == "params_terma_w1":
            ### Event window/W1
            try:
                if self.terma_w1_line_edit.text() is not None:
                    self.params_terma_w1 = float(self.terma_w1_line_edit.text())
                    self.terma_w1_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: TERMA event window value could not be converted to float")
                self.params_terma_w1 = 1
                self.terma_w1_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_terma_w2":
            ### Event window/W1
            try:
                if self.terma_w2_line_edit.text() is not None:
                    self.params_terma_w2 = float(self.terma_w2_line_edit.text())
                    self.terma_w2_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: TERMA cycle window value could not be converted to float")
                self.params_terma_w2 = 4
                self.terma_w2_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_terma_width":
            ### EWMA cross coefficient
            try:
                if self.terma_width_line_edit.text() is not None:
                    self.params_terma_width = float(self.terma_width_line_edit.text())
                    self.terma_width_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: TERMA width could not be converted to float")
                self.params_terma_width = 1
                self.terma_width_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_terma_b":
            ### SRMAC width
            try:
                if self.terma_b_line_edit.text() is not None:
                    self.params_terma_b = float(self.terma_b_line_edit.text())
                    self.terma_b_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: TERMA b coefficient could not be converted to Float")
                self.params_terma_b = 0
                self.terma_b_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)

    def update_findpeaks_parameters(self, param):
        if param == "params_findpeaks_prominence":
            ### Prominence
            try:
                if self.findpeaks_prominence_line_edit.text() is not None:
                    self.params_findpeaks_prominence = float(self.findpeaks_prominence_line_edit.text())
                    self.findpeaks_prominence_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: Find_peaks prominence value could not be converted to float")
                self.params_findpeaks_prominence = 0.1
                self.findpeaks_prominence_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_findpeaks_heval":
            ### Heval
            try:
                if self.findpeaks_heval_line_edit.text() is not None:
                    self.params_findpeaks_heval = float(self.findpeaks_heval_line_edit.text())
                    self.findpeaks_heval_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: Find_peaks heval value could not be converted to float")
                self.params_findpeaks_heval = 0.8
                self.findpeaks_heval_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_findpeaks_width":
            ### Width
            try:
                if self.findpeaks_width_line_edit.text() is not None:
                    self.params_findpeaks_width = float(self.findpeaks_width_line_edit.text())
                    self.findpeaks_width_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: Find_peaks width value could not be converted to float")
                self.params_findpeaks_width = 0.5
                self.findpeaks_width_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_findpeaks_proximity":
            ### Proximity
            try:
                if self.findpeaks_proxim_line_edit.text() is not None:
                    self.params_findpeaks_proximity = float(self.findpeaks_proxim_line_edit.text())
                    self.findpeaks_proxim_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: Find_peaks proximity value could not be converted to float")
                self.params_findpeaks_proximity = 0.5
                self.findpeaks_proxim_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)

    def update_cwtoa_parameters(self, param):
        if param == "params_cwt_fmin":
            ### Fmin
            try:
                if self.cwtoa_fmin_line_edit.text() is not None:
                    self.params_cwtoa_f_min = float(self.cwtoa_fmin_line_edit.text())
                    self.cwtoa_fmin_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: CWT minimum frequency value could not be converted to float")
                self.params_cwtoa_f_min = 0.1
                self.cwtoa_fmin_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_cwt_fmax":
            ### Fmax
            try:
                if self.cwtoa_fmax_line_edit.text() is not None:
                    self.params_cwtoa_f_max = float(self.cwtoa_fmax_line_edit.text())
                    self.cwtoa_fmax_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: CWT maximum frequency value could not be converted to float")
                self.params_cwtoa_f_max = 1
                self.cwtoa_fmax_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_cwt_scales":
            ### Scales
            try:
                if self.cwtoa_scales_line_edit.text() is not None:
                    self.params_cwtoa_scales = int(self.cwtoa_scales_line_edit.text())
                    self.cwtoa_scales_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: CWT minimum frequency value could not be converted to int")
                self.params_cwtoa_scales = 3
                self.cwtoa_scales_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_cwt_width":
            ### Width
            try:
                if self.cwtoa_width_line_edit.text() is not None:
                    self.params_cwtoa_width = float(self.cwtoa_width_line_edit.text())
                    self.cwtoa_width_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: CWT minimum frequency value could not be converted to float")
                self.params_cwtoa_width = 1
                self.cwtoa_width_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)
        elif param == "params_cwt_th":
            ### Threshold
            try:
                if self.cwtoa_th_line_edit.text() is not None:
                    self.params_cwtoa_th = float(self.cwtoa_th_line_edit.text())
                    self.cwtoa_th_line_edit.setStyleSheet(self.stylesheet_lineedit_valid)
            except:
                print("Error: CWT minimum frequency value could not be converted to float")
                self.params_cwtoa_th = 0
                self.cwtoa_th_line_edit.setStyleSheet(self.stylesheet_lineedit_invalid)

    ### Wrapper function for calling respiratory rate algorithms
    def run_algorithm(self, data):
        if self.algorithm == "srmac":
            rr, peaks = srmac(data, self.sample_rate, len(data), coef_fast=self.params_srmac_coef_fast, coef_slow=self.params_srmac_coef_slow, coef_cross=self.params_srmac_coef_cross, width=self.params_srmac_width, threshold=self.params_srmac_th)
        elif self.algorithm == "terma":
            rr, peaks = terma(data, self.sample_rate, len(data), window_event=self.params_terma_w1, window_cycle=self.params_terma_w2, b_coef=self.params_terma_b, width=self.params_terma_width)
        elif self.algorithm == "find_peaks":
            rr, peaks = find_peaks(data, self.sample_rate, len(data), prominence=self.params_findpeaks_prominence, heval_ratio=self.params_findpeaks_heval, width=self.params_findpeaks_width, proximity=self.params_findpeaks_proximity)
        elif self.algorithm == "cwt":
            rr, peaks = cwt_peaks(data, self.sample_rate, len(data), min_freq=self.params_cwtoa_f_min, max_freq=self.params_cwtoa_f_max, resolution=self.params_cwtoa_scales, width=self.params_cwtoa_width, threshold=self.params_cwtoa_th)
        else:
            rr, peaks = 0, []

        return rr, peaks

    def algorithm_wrapper(self):
        ### If analyze is true, run the selected respiratory rate algorithm on this data window
        if self.analyze:
            time_window = self.time[self.data_index:self.data_index+self.data_window_size*self.sample_rate]
            data_window = self.data[self.data_index:self.data_index+self.data_window_size*self.sample_rate]

            rr, peaks = self.run_algorithm(z_norm(data_window))

            ### Show peaks
            if self.action_show_peaks.isChecked():
                self.peak_scatter.setData(time_window[peaks], data_window[peaks])

            ### Show respiratory rate
            #...
        else:
            ### Clear peak plot
            self.peak_scatter.setData([],[])

    ### Update the plot and related items
    def update_plot(self):
        ### Plot data
        time_window = self.time[self.data_index:self.data_index+self.data_window_size*self.sample_rate]
        data_window = self.data[self.data_index:self.data_index+self.data_window_size*self.sample_rate]
        self.plot_signal.setData(time_window, data_window)

        ### Set plot range
        self.plot_main.setYRange(min(data_window)-abs(min(data_window))/2, max(data_window)+abs(max(data_window))/2)
        self.plot_main.setXRange(min(time_window), max(time_window))

        ### Update the midline position
        self.midline.setX(time_window[len(time_window)//2])

        ### Update the progress bar value
        self.progress_bar.setValue(self.data_index)

        ### Update the sample time indicator
        # self.time_indicator_samples.setText(str(round(self.data_index/self.sample_rate, 2))+" / "+str(round(self.max_data_index/self.sample_rate, 2)))
        self.time_indicator_samples.setText(str(self.data_index)+" / "+str(self.max_data_index))
    
    ### When data is loaded for the first time, enable all features that start disabled
    def enable_widgets(self):
        self.button_play.setEnabled(True)
        self.button_reset.setEnabled(True)
        self.menu_edit.setEnabled(True)
        self.progress_bar.setEnabled(True)
        self.algorithm_container.setEnabled(True)

    ### Toggle data playback
    def play_button_pressed(self):
        if self.file_loaded:
            if self.running == True:
                self.running = False
            else:
                self.running = True

    ### Reset data index and stop plotting
    def reset_button_pressed(self):
        self.running = False
        self.data_index = 0
        self.update_plot()

    def default_button_pressed(self):
        if self.algorithm == "srmac":
            self.srmac_coef_fast_line_edit.setText("0.0508")
            self.srmac_coef_slow_line_edit.setText("0.0103")
            self.srmac_coef_cross_line_edit.setText("0.91")
            self.srmac_width_line_edit.setText("0.5")
            self.srmac_th_line_edit.setText("0")
        elif self.algorithm == "terma":
            self.terma_w1_line_edit.setText("1")
            self.terma_w2_line_edit.setText("4")
            self.terma_width_line_edit.setText("1")
            self.terma_b_line_edit.setText("0")
        elif self.algorithm == "find_peaks":
            self.findpeaks_prominence_line_edit.setText("0.1")
            self.findpeaks_heval_line_edit.setText("0.8")
            self.findpeaks_width_line_edit.setText("0.5")
            self.findpeaks_proxim_line_edit.setText("0.5")
        elif self.algorithm == "cwt":
            self.cwtoa_fmin_line_edit.setText("0.1")
            self.cwtoa_fmax_line_edit.setText("1")
            self.cwtoa_scales_line_edit.setText("3")
            self.cwtoa_width_line_edit.setText("1")
            self.cwtoa_th_line_edit.setText("0")
        else:
            ErrorDialog(self, "Error: invalid algorithm selected").exec()

    ### Process analyze checkbox toggle
    def analyze_check_toggled(self):
        if self.analyze_check.isChecked():
            self.analyze = True
        else:
            self.analyze = False
        self.update_plot()

    ### Process progress bar drag event
    def dragging_progress_bar(self):
        ### Update the data index
        self.data_index = self.progress_bar.value()

        ### Update the plot once, but only if playback is stopped
        if not self.running:
            self.update_plot()
            self.algorithm_wrapper()

    ### Change the selected respiratory rate algorithm
    def change_algorithm(self):
        self.algorithm = self.algorithms[self.algselector.currentIndex()]

        ### Hide all algorithm parameter containers
        self.srmac_parameter_container.hide()
        self.terma_parameter_container.hide()
        self.findpeaks_parameter_container.hide()
        self.cwtoa_parameter_container.hide()

        ### Show only the one that has been selected
        if self.algorithm == "srmac":
            self.srmac_parameter_container.show()
        elif self.algorithm == "terma":
            self.terma_parameter_container.show()
        elif self.algorithm == "find_peaks":
            self.findpeaks_parameter_container.show()
        elif self.algorithm == "cwt":
            self.cwtoa_parameter_container.show()

    ### Toggle algorithm bar visibility
    def toggle_algorithm_bar(self):
        if self.action_show_algs.isChecked():
            self.algorithm_container.show()
        else:
            self.algorithm_container.hide()
        self.change_algorithm()

    ### Toggle midline visibility
    def toggle_midline(self):
        if self.action_show_midline.isChecked():
            self.midline.show()
        else:
            self.midline.hide()

    ### Toggle whether or not peaks are displayed
    def toggle_peaks(self):
        if self.action_show_peaks.isChecked():
            self.peak_scatter.show()
        else:
            self.peak_scatter.hide()

    ### Toggle whether or not the signal is displayed
    def toggle_signal(self):
        if self.action_show_signal.isChecked():
            self.plot_signal.show()
        else:
            self.plot_signal.hide()

    ### Show general info about the application
    def show_about(self):
        if self.running:
            self.play_button_pressed()
        AboutDialog(self).exec()

    ### Open the dialog for opening a new file
    def button_open_new_file(self):
        ### Open the native file selection dialog window
        dialog = QFileDialog(self)
        selected_file = dialog.getOpenFileName(filter="*.csv")
        dialog.close()

        ### If valid pass the selected file to the custom file dialog window for reading data
        selected_file = selected_file[0]
        if len(selected_file) == 0:
            ErrorDialog(self, "Error: No file selected.").exec()
            return
        if selected_file.split(".")[-1] != "csv":
            ErrorDialog(self, "Error: Incorrect file format.").exec()
            return
        else:
            dialog_open_file = ReadFileDialog(self, selected_file)
            if dialog_open_file.exec():
                ### Retrieve all necessary data from the dialog before it is closed
                self.data_raw = dialog_open_file.ydata
                self.time_raw = dialog_open_file.xdata
                self.sample_rate = dialog_open_file.sample_rate

                ### Validate the data
                if len(self.data_raw) != len(self.time_raw):
                    ErrorDialog(self, "Error: data and time arrays are of different sizes.").exec()
                    return
                if len(self.data_raw) < self.data_window_size*self.sample_rate:
                    ErrorDialog(self, "Error: data array too short.").exec()
                    return

                ### If this is the first file to be loaded, finish GUI initialization
                if not self.file_loaded:
                    self.file_loaded = True
                    self.enable_widgets()

                ### Update data
                self.reset_data()
                return
        return

    ### Increment the data_index
    def increment_index(self):
        if self.running and self.file_loaded:
            try:
                if self.data_index < self.max_data_index:
                    self.data_index += 1
                else:
                    self.running = False
                    self.update_plot()
            except:
                self.running = False
                self.update_plot()

    ### Process keypresses
    def keyPressEvent(self, event):
        if type(event) == QKeyEvent:
            if event.key() == Qt.Key.Key_Q:
                if self.file_loaded:
                    self.play_button_pressed()
            elif event.key() == Qt.Key.Key_R:
                if self.file_loaded:
                    self.reset_button_pressed()
            elif event.key() == Qt.Key.Key_O:
                self.button_open_new_file()

    ### Main application loop
    def main_loop(self):
        ### Increment the tick counter
        self.tick += 1

        ### Update plots, process data etc.
        if self.running and self.file_loaded:
            try:
                ### Update plot
                self.update_plot()

                ### Analyze data to find peaks
                if self.tick % self.algorithm_delta_ticks == 0:
                    self.algorithm_wrapper()

                ### Update play button text
                self.button_play.setText("Pause (Q)")
            except Exception as e:
                ErrorDialog(self, "Error: "+str(e)).exec()
                self.running = False
        else:
            ### Update play button text
            self.button_play.setText("Start (Q)")

        ### Do these need to be here?
        QApplication.processEvents()
        self.show()

    ### Override the closeEvent method to add a confirmation dialog
    def closeEvent(self, event):
        self.running = False
        if QMessageBox.question(self, "", "Are you sure you want to quit?", QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No) == QMessageBox.StandardButton.Yes:
            event.ignore()
            self.close()
            sys.exit()
        else:
            event.ignore()

### Launch application
app = QApplication(sys.argv)
window = MainWindow()
window.showNormal()
sys.exit(app.exec())

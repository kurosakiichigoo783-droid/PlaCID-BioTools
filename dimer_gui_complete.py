import sys
import csv
import json # Import json for parameter saving/loading
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
import os
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTextEdit, QPushButton, QTableWidget, QTableWidgetItem, QFormLayout,
    QLineEdit, QLabel, QGroupBox, QHeaderView, QAbstractItemView,
    QMessageBox, QScrollArea, QSplitter, QFrame, QTabWidget, QFileDialog, QCheckBox,
    QDialog, QInputDialog, QDialogButtonBox
)
#from email_dialog import EmailDialog
from PyQt5.QtGui import QFont, QColor, QPalette
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from Bio import Entrez, SeqIO
import io
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from typing import Dict, List # Import Dict and List

from docking_backend import PrimerDesignerBackend

# Default parameters, reflecting a professional standard
DEFAULT_PARAMS = {
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 58.0,
    'PRIMER_MAX_TM': 62.0,
    'PRIMER_MAX_DIFF_TM': 2.0,
    'PRIMER_OPT_GC_PERCENT': 50.0,
    'PRIMER_MIN_GC': 40.0,
    'PRIMER_MAX_GC': 60.0,
    'PRIMER_PRODUCT_SIZE_RANGE': '100-300',
    'PRIMER_NUM_RETURN': 5,
    'NA_CONC': 50.0, # mM
    'MG_CONC': 1.5, # mM
    'DNTP_CONC': 0.2, # mM
    'PRIMER_CONC': 250.0, # nM
    'MASK_SNPS': False,
    'MASK_REPEATS': False,
    'TARGET_START': '',
    'TARGET_END': '',
}

# --- Worker Threads ---

class NCBIFetchWorker(QThread):
    """
    Worker thread to fetch sequences from NCBI without freezing the GUI.
    """
    finished = pyqtSignal(object)
    error = pyqtSignal(str)

    def __init__(self, gene_id, email):
        super().__init__()
        self.gene_id = gene_id
        self.email = email

    def run(self):
        """
        Execute the NCBI fetch process.
        """
        try:
            Entrez.email = self.email
            handle = Entrez.efetch(db="nucleotide", id=self.gene_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            self.finished.emit(str(record.seq))
        except Exception as e:
            self.error.emit(str(e))

class PrimerDesignWorker(QThread):
    """
    Worker thread to run primer design without freezing the GUI.
    """
    finished = pyqtSignal(object)
    error = pyqtSignal(str)

    def __init__(self, sequence, params):
        super().__init__()
        self.sequence = sequence
        self.params = params
        self.backend = PrimerDesignerBackend(params)

    def run(self):
        """
        Execute the primer design process using the backend.
        """
        try:
            results = self.backend.design_primers(self.sequence)
            self.finished.emit(results)
        except Exception as e:
            self.error.emit(str(e))

# --- Dialogs ---

class EmailDialog(QDialog):
    """
    Dialog for sending email with results.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Email Results")
        layout = QFormLayout(self)

        self.recipient_input = QLineEdit()
        self.subject_input = QLineEdit()
        self.body_input = QTextEdit()
        self.smtp_host_input = QLineEdit()
        self.smtp_port_input = QLineEdit()
        self.smtp_user_input = QLineEdit()
        self.smtp_pass_input = QLineEdit()
        self.smtp_pass_input.setEchoMode(QLineEdit.Password)

        layout.addRow(QLabel("Recipient:"), self.recipient_input)
        layout.addRow(QLabel("Subject:"), self.subject_input)
        layout.addRow(QLabel("Body:"), self.body_input)
        layout.addRow(QLabel("SMTP Host:"), self.smtp_host_input)
        layout.addRow(QLabel("SMTP Port:"), self.smtp_port_input)
        layout.addRow(QLabel("SMTP User:"), self.smtp_user_input)
        layout.addRow(QLabel("SMTP Password:"), self.smtp_pass_input)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def get_details(self):
        return {
            "recipient": self.recipient_input.text(),
            "subject": self.subject_input.text(),
            "body": self.body_input.toPlainText(),
            "smtp_host": self.smtp_host_input.text(),
            "smtp_port": int(self.smtp_port_input.text()),
            "smtp_user": self.smtp_user_input.text(),
            "smtp_pass": self.smtp_pass_input.text(),
        }

# --- Frontend GUI ---

class PrimerDesignApp(QMainWindow):
    """
    Main application window for the Primer Designer.
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Automated Primer Designer")
        self.setGeometry(100, 100, 1600, 900) # Increased size for new panels
        self.results_data = None
        self.setup_ui()

    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        left_panel = QFrame()
        left_panel.setFrameShape(QFrame.StyledPanel)
        left_layout = QVBoxLayout(left_panel)
        input_tabs = QTabWidget()
        self.create_paste_tab(input_tabs)
        self.create_file_tab(input_tabs)
        self.create_ncbi_tab(input_tabs)
        left_layout.addWidget(input_tabs)
        
        # Parameters
        params_box = QGroupBox("Design Parameters")
        params_layout = QFormLayout()
        self.param_inputs = {}
        for key, value in DEFAULT_PARAMS.items():
            if isinstance(value, bool):
                checkbox = QCheckBox(key.replace('_', ' ').title())
                checkbox.setChecked(value)
                self.param_inputs[key] = checkbox
                params_layout.addRow(checkbox)
            else:
                line_edit = QLineEdit(str(value))
                self.param_inputs[key] = line_edit
                params_layout.addRow(QLabel(key.replace('_', ' ').title() + ":"), line_edit)
        
        params_box.setLayout(params_layout)
        scroll_area = QScrollArea()
        scroll_area.setWidget(params_box)
        scroll_area.setWidgetResizable(True)
        left_layout.addWidget(scroll_area)

        # Save/Load buttons for parameters
        param_buttons_layout = QHBoxLayout()
        save_params_button = QPushButton("Save Parameters")
        save_params_button.clicked.connect(self.save_parameters)
        param_buttons_layout.addWidget(save_params_button)
        load_params_button = QPushButton("Load Parameters")
        load_params_button.clicked.connect(self.load_parameters)
        param_buttons_layout.addWidget(load_params_button)
        left_layout.addLayout(param_buttons_layout)

        self.run_button = QPushButton("Design Primers")
        self.run_button.clicked.connect(self.run_primer_design)
        left_layout.addWidget(self.run_button)
        splitter.addWidget(left_panel)

        # Right Panel - Tab Widget for Results, Viz, Stats
        right_tabs = QTabWidget()
        self.create_results_tab(right_tabs)
        self.create_visualization_tab(right_tabs)
        self.create_analysis_tab(right_tabs)
        splitter.addWidget(right_tabs)
        splitter.setSizes([400, 1200]) # Initial size distribution

    def create_paste_tab(self, tabs):
        paste_tab = QWidget()
        layout = QVBoxLayout(paste_tab)
        self.sequence_input = QTextEdit()
        self.sequence_input.setPlaceholderText("Paste your DNA sequence here (FASTA format or raw sequence)...")
        layout.addWidget(self.sequence_input)
        tabs.addTab(paste_tab, "Paste Sequence")

    def create_file_tab(self, tabs):
        file_tab = QWidget()
        layout = QVBoxLayout(file_tab)
        self.file_path_label = QLabel("No file selected.")
        layout.addWidget(self.file_path_label)
        browse_button = QPushButton("Browse...")
        browse_button.clicked.connect(self.browse_file)
        layout.addWidget(browse_button)
        tabs.addTab(file_tab, "File Upload")

    def create_ncbi_tab(self, tabs):
        ncbi_tab = QWidget()
        layout = QFormLayout(ncbi_tab)
        self.email_input = QLineEdit()
        self.email_input.setPlaceholderText("your.email@example.com")
        layout.addRow(QLabel("Your Email:"), self.email_input)
        self.gene_id_input = QLineEdit()
        self.gene_id_input.setPlaceholderText("e.g., NM_000546.6 for TP53")
        layout.addRow(QLabel("Gene ID/Accession:"), self.gene_id_input)
        self.fetch_button = QPushButton("Fetch Sequence")
        self.fetch_button.clicked.connect(self.fetch_sequence_from_ncbi)
        layout.addWidget(self.fetch_button)
        tabs.addTab(ncbi_tab, "NCBI Fetch")

    def create_results_tab(self, tabs):
        results_tab = QWidget()
        results_layout = QVBoxLayout(results_tab)
        
        results_splitter = QSplitter(Qt.Vertical)
        
        results_box = QGroupBox("Primer Design Results")
        table_layout = QVBoxLayout()
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(15) 
        self.results_table.setHorizontalHeaderLabels([
            "Rank", "Quality", "Custom Score", "Fwd Sequence", "Fwd Tm", "Fwd GC%", 
            "Fwd Hairpin dG", "Fwd Self-Dimer dG", 
            "Rev Sequence", "Rev Tm", "Rev GC%", 
            "Rev Hairpin dG", "Rev Self-Dimer dG", 
            "Cross-Dimer dG", "Product Size"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.results_table.selectionModel().selectionChanged.connect(self.update_visualization)
        table_layout.addWidget(self.results_table)
        
        export_layout = QHBoxLayout()
        self.export_csv_button = QPushButton("Export to CSV")
        self.export_csv_button.clicked.connect(self.export_to_csv)
        self.export_fasta_button = QPushButton("Export to FASTA")
        self.export_fasta_button.clicked.connect(self.export_to_fasta)
        # self.email_results_button = QPushButton("Email Results")
        # self.email_results_button.clicked.connect(self.email_results)
        export_layout.addWidget(self.export_csv_button)
        export_layout.addWidget(self.export_fasta_button)
        # export_layout.addWidget(self.email_results_button)
        table_layout.addLayout(export_layout)

        results_box.setLayout(table_layout)
        results_splitter.addWidget(results_box)
        results_layout.addWidget(results_splitter)
        
        tabs.addTab(results_tab, "Results Table")

    def create_visualization_tab(self, tabs):
        vis_tab = QWidget()
        vis_layout = QVBoxLayout(vis_tab)
        vis_box = QGroupBox("Primer Visualization")
        vis_content_layout = QVBoxLayout()
        self.vis_panel = QTextEdit()
        self.vis_panel.setReadOnly(True)
        self.vis_panel.setFont(QFont("Courier", 10))
        vis_content_layout.addWidget(self.vis_panel)
        vis_box.setLayout(vis_content_layout)
        vis_layout.addWidget(vis_box)
        tabs.addTab(vis_tab, "Primer Visualization")

    def create_analysis_tab(self, tabs):
        analysis_tab = QWidget()
        analysis_layout = QVBoxLayout(analysis_tab)
        
        # Summary Statistics Table
        stats_box = QGroupBox("Summary Statistics (All Candidates)")
        stats_layout = QVBoxLayout()
        self.stats_table = QTableWidget()
        self.stats_table.setColumnCount(6)
        self.stats_table.setHorizontalHeaderLabels(["Metric", "Mean", "Median", "Std Dev", "Min", "Max"])
        self.stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        stats_layout.addWidget(self.stats_table)
        stats_box.setLayout(stats_layout)
        analysis_layout.addWidget(stats_box)

        # Plotting Area
        plots_box = QGroupBox("Distribution Plots (All Candidates)")
        plots_layout = QVBoxLayout()
        self.figure, self.ax = plt.subplots(figsize=(10, 6))
        self.canvas = FigureCanvas(self.figure)
        plots_layout.addWidget(self.canvas)
        plots_box.setLayout(plots_layout)
        analysis_layout.addWidget(plots_box)
        
        tabs.addTab(analysis_tab, "Analysis & Statistics")


    def browse_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Open Sequence File", "", "All Files (*);;FASTA Files (*.fa *.fasta);;GenBank Files (*.gb *.gbk)")
        if file_path:
            self.file_path_label.setText(file_path)
            try:
                with open(file_path, 'r') as f:
                    content = f.read()
                    try:
                        record = SeqIO.read(io.StringIO(content), "genbank")
                        sequence = str(record.seq)
                    except ValueError:
                        try:
                            record = SeqIO.read(io.StringIO(content), "fasta")
                            sequence = str(record.seq)
                        except ValueError:
                            sequence = content.split('\n', 1)[-1].replace('\n', '') if '>' in content else content
                self.sequence_input.setText(sequence)
                self.show_info("File loaded successfully.")
            except Exception as e:
                self.show_error(f"Error reading file: {e}")

    def fetch_sequence_from_ncbi(self):
        email = self.email_input.text().strip()
        if not email or "@" not in email:
            self.show_error("Please enter a valid email address for NCBI Entrez.")
            return
        
        gene_id = self.gene_id_input.text().strip()
        if not gene_id:
            self.show_error("Please enter a Gene ID or Accession number.")
            return

        self.fetch_button.setText("Fetching...")
        self.fetch_button.setEnabled(False)
        self.fetch_worker = NCBIFetchWorker(gene_id, email)
        self.fetch_worker.finished.connect(self.populate_sequence_from_fetch)
        self.fetch_worker.error.connect(self.fetch_error)
        self.fetch_worker.start()

    def populate_sequence_from_fetch(self, sequence):
        self.sequence_input.setText(sequence)
        self.show_info("Sequence fetched successfully.")
        self.fetch_button.setText("Fetch Sequence")
        self.fetch_button.setEnabled(True)

    def fetch_error(self, message):
        self.show_error(f"Failed to fetch sequence from NCBI: {message}")
        self.fetch_button.setText("Fetch Sequence")
        self.fetch_button.setEnabled(True)

    def get_current_params(self):
        """Collects current parameter values from UI widgets."""
        current_params = {}
        for key, value_widget in self.param_inputs.items():
            if isinstance(value_widget, QCheckBox):
                current_params[key] = value_widget.isChecked()
            else:
                value_text = value_widget.text()
                if key in ['TARGET_START', 'TARGET_END'] and value_text == '':
                    current_params[key] = None
                elif '.' in value_text:
                    try:
                        current_params[key] = float(value_text)
                    except ValueError:
                        current_params[key] = value_text # Keep as string if not a valid float
                elif value_text.isdigit():
                    try:
                        current_params[key] = int(value_text)
                    except ValueError:
                        current_params[key] = value_text # Keep as string if not a valid int
                else:
                    current_params[key] = value_text
        return current_params

    def set_params_from_dict(self, params_dict: Dict):
        """Sets UI widgets with values from a dictionary."""
        for key, value in params_dict.items():
            if key in self.param_inputs:
                widget = self.param_inputs[key]
                if isinstance(widget, QCheckBox):
                    if isinstance(value, bool):
                        widget.setChecked(value)
                elif isinstance(widget, QLineEdit):
                    widget.setText(str(value))
            
    def save_parameters(self):
        params_to_save = self.get_current_params()
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Parameters", "primer_params.json", "JSON Files (*.json)")
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    json.dump(params_to_save, f, indent=4)
                self.show_info(f"Parameters saved to {file_path}")
            except Exception as e:
                self.show_error(f"Error saving parameters: {e}")

    def load_parameters(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Parameters", "", "JSON Files (*.json)")
        if file_path:
            try:
                with open(file_path, 'r') as f:
                    loaded_params = json.load(f)
                self.set_params_from_dict(loaded_params)
                self.show_info(f"Parameters loaded from {file_path}")
            except Exception as e:
                self.show_error(f"Error loading parameters: {e}")

    def run_primer_design(self):
        sequence = self.sequence_input.toPlainText().strip()
        if sequence.startswith('>'):
            parts = sequence.split('\n', 1)
            sequence = parts[1].replace('\n', '') if len(parts) > 1 else ''

        params = self.get_current_params() # Get parameters directly from UI

        self.run_button.setText("Designing...")
        self.run_button.setEnabled(False)
        self.worker = PrimerDesignWorker(sequence, params)
        self.worker.finished.connect(self.display_results)
        self.worker.error.connect(self.show_error)
        self.worker.start()

    def display_results(self, results):
        self.results_data = results
        self.results_table.setRowCount(0)
        num_primers = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
        if num_primers == 0:
            self.show_error("No primer pairs found. Try relaxing parameters.")
            self.run_button.setText("Design Primers")
            self.run_button.setEnabled(True)
            self.display_statistics(results.get('summary_statistics', {}))
            self.plot_distributions(results.get('all_candidate_metrics', []))
            return

        self.results_table.setRowCount(num_primers)
        self.results_table.setColumnCount(15) # Ensure column count matches headers
        self.results_table.setHorizontalHeaderLabels([
            "Rank", "Quality", "Custom Score", "Fwd Sequence", "Fwd Tm", "Fwd GC%", 
            "Fwd Hairpin dG", "Fwd Self-Dimer dG", 
            "Rev Sequence", "Rev Tm", "Rev GC%", 
            "Rev Hairpin dG", "Rev Self-Dimer dG", 
            "Cross-Dimer dG", "Product Size"
        ])
        
        for i in range(num_primers):
            quality_score = results.get(f'PRIMER_PAIR_{i}_CUSTOM_SCORE', 0)
            quality_badge = results.get(f'PRIMER_PAIR_{i}_QUALITY_BADGE', 'N/A')
            
            self.results_table.setItem(i, 0, QTableWidgetItem(str(i + 1)))
            self.results_table.setItem(i, 1, QTableWidgetItem(quality_badge))
            self.results_table.setItem(i, 2, QTableWidgetItem(f"{quality_score:.2f}"))
            self.results_table.setItem(i, 3, QTableWidgetItem(results[f'PRIMER_LEFT_{i}_SEQUENCE']))
            self.results_table.setItem(i, 4, QTableWidgetItem(f"{results[f'PRIMER_LEFT_{i}_TM']:.2f}"))
            self.results_table.setItem(i, 5, QTableWidgetItem(f"{results[f'PRIMER_LEFT_{i}_GC_PERCENT']:.2f}"))
            self.results_table.setItem(i, 6, QTableWidgetItem(f"{results.get(f'PRIMER_LEFT_{i}_HAIRPIN_DG', 0):.2f}"))
            self.results_table.setItem(i, 7, QTableWidgetItem(f"{results.get(f'PRIMER_LEFT_{i}_SELF_DIMER_DG', 0):.2f}"))
            self.results_table.setItem(i, 8, QTableWidgetItem(results[f'PRIMER_RIGHT_{i}_SEQUENCE']))
            self.results_table.setItem(i, 9, QTableWidgetItem(f"{results[f'PRIMER_RIGHT_{i}_TM']:.2f}"))
            self.results_table.setItem(i, 10, QTableWidgetItem(f"{results[f'PRIMER_RIGHT_{i}_GC_PERCENT']:.2f}"))
            self.results_table.setItem(i, 11, QTableWidgetItem(f"{results.get(f'PRIMER_RIGHT_{i}_HAIRPIN_DG', 0):.2f}"))
            self.results_table.setItem(i, 12, QTableWidgetItem(f"{results.get(f'PRIMER_RIGHT_{i}_SELF_DIMER_DG', 0):.2f}"))
            self.results_table.setItem(i, 13, QTableWidgetItem(f"{results.get(f'PRIMER_PAIR_{i}_CROSS_DIMER_DG', 0):.2f}"))
            self.results_table.setItem(i, 14, QTableWidgetItem(str(results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'])))

        self.run_button.setText("Design Primers")
        self.run_button.setEnabled(True)
        if num_primers > 0:
            self.results_table.selectRow(0)
        
        # Update statistics and plots
        self.display_statistics(results.get('summary_statistics', {}))
        self.plot_distributions(results.get('all_candidate_metrics', []))


    def update_visualization(self, selected, deselected):
        if not self.results_data or not selected.indexes():
            self.vis_panel.setHtml("")
            return
        
        row = selected.indexes()[0].row()
        template = self.sequence_input.toPlainText().strip()
        if template.startswith('>'):
            template = template.split('\n', 1)[1].replace('\n', '')

        # primer3 returns 0-based start, and length (e.g., [0, 20]) for forward
        fwd_start_in_template = self.results_data[f'PRIMER_LEFT_{row}'][0]
        fwd_len = self.results_data[f'PRIMER_LEFT_{row}'][1]
        
        # primer3 returns 0-based start of the reverse primer on the template, and length (e.g., [200, 20])
        # The 'start' is the 0-based index of the rightmost base of the oligo on the template.
        rev_end_in_template = self.results_data[f'PRIMER_RIGHT_{row}'][0]
        rev_len = self.results_data[f'PRIMER_RIGHT_{row}'][1]
        rev_start_in_template = rev_end_in_template - rev_len + 1 # Calculate true start on template

        html = f"""<pre><font face=\"Courier\" size=\"2\">
               <font color=\"gray\">{template[:fwd_start_in_template]}</font>\
               <font color=\"green\" style=\"background-color:#2E8B57;\">{template[fwd_start_in_template : fwd_start_in_template + fwd_len]}</font>\
               <font color=\"gray\">{template[fwd_start_in_template + fwd_len : rev_start_in_template]}</font>\
               <font color=\"red\" style=\"background-color:#A52A2A;\">{template[rev_start_in_template : rev_start_in_template + rev_len]}</font>\
               <font color=\"gray\">{template[rev_start_in_template + rev_len:]}</font>\
               </font></pre>"""
        
        # Add target region highlighting
        target_start_param = self.param_inputs['TARGET_START'].text()
        target_end_param = self.param_inputs['TARGET_END'].text()

        if target_start_param and target_end_param:
            try:
                t_start = int(target_start_param)
                t_end = int(target_end_param)
                if 0 <= t_start < t_end <= len(template):
                    # Create a list of (text, color) tuples for flexible highlighting
                    highlight_segments = []
                    current_pos = 0

                    def add_segment(text, color):
                        if text:
                            highlight_segments.append(f'<font color="{color}">{text}</font>')

                    # Before target region
                    add_segment(template[current_pos:t_start], "gray")
                    current_pos = t_start

                    # Within target region
                    # Before forward primer in target
                    add_segment(template[current_pos:min(fwd_start_in_template, t_end)], "cyan")
                    current_pos = max(current_pos, fwd_start_in_template)

                    # Forward primer
                    if current_pos < t_end:
                        add_segment(template[current_pos:min(current_pos + fwd_len, t_end)], "green")
                    current_pos += fwd_len
                    
                    # Between forward and reverse primer in target
                    if current_pos < t_end:
                        add_segment(template[current_pos:min(rev_start_in_template, t_end)], "cyan")
                    current_pos = max(current_pos, rev_start_in_template)

                    # Reverse primer
                    if current_pos < t_end:
                        add_segment(template[current_pos:min(current_pos + rev_len, t_end)], "red")
                    current_pos += rev_len

                    # After reverse primer in target
                    if current_pos < t_end:
                        add_segment(template[current_pos:t_end], "cyan")
                    current_pos = max(current_pos, t_end)

                    # After target region
                    add_segment(template[current_pos:], "gray")

                    html = f"""<pre><font face=\"Courier\" size=\"2\">{''.join(highlight_segments)}</font></pre>"""


            except ValueError:
                pass # Invalid target region input, fall back to just primer highlighting
        
        self.vis_panel.setHtml(html)

    def display_statistics(self, stats: Dict):
        self.stats_table.setRowCount(0)
        if not stats:
            return

        row = 0
        for metric, data in stats.items():
            self.stats_table.setRowCount(row + 1)
            self.stats_table.setItem(row, 0, QTableWidgetItem(metric.replace('_', ' ').title()))
            self.stats_table.setItem(row, 1, QTableWidgetItem(f"{data['mean']:.2f}"))
            self.stats_table.setItem(row, 2, QTableWidgetItem(f"{data['median']:.2f}"))
            self.stats_table.setItem(row, 3, QTableWidgetItem(f"{data['std_dev']:.2f}"))
            self.stats_table.setItem(row, 4, QTableWidgetItem(f"{data['min']:.2f}"))
            self.stats_table.setItem(row, 5, QTableWidgetItem(f"{data['max']:.2f}"))
            row += 1

    def plot_distributions(self, all_candidates: List[Dict]):
        self.figure.clear()
        if not all_candidates:
            self.canvas.draw()
            return
        
        # Create a subplot grid for multiple plots
        gs = self.figure.add_gridspec(2, 2, hspace=0.4, wspace=0.3)
        
        # Plot 1: Histogram of Custom Scores
        scores = [c['CUSTOM_SCORE'] for c in all_candidates]
        if scores:
            ax1 = self.figure.add_subplot(gs[0, 0])
            ax1.hist(scores, bins=20, color='skyblue', edgecolor='black')
            ax1.set_title('Custom Scores', fontsize=10)
            ax1.set_xlabel('Score', fontsize=8)
            ax1.set_ylabel('Frequency', fontsize=8)
            ax1.tick_params(axis='both', which='major', labelsize=7)

        # Plot 2: Histogram of Forward Primer Tm
        fwd_tms = [c['PRIMER_LEFT_TM'] for c in all_candidates]
        if fwd_tms:
            ax2 = self.figure.add_subplot(gs[0, 1])
            ax2.hist(fwd_tms, bins=20, color='lightcoral', edgecolor='black')
            ax2.set_title('Fwd Primer Tm', fontsize=10)
            ax2.set_xlabel('Tm (°C)', fontsize=8)
            ax2.set_ylabel('Frequency', fontsize=8)
            ax2.tick_params(axis='both', which='major', labelsize=7)

        # Plot 3: Histogram of Forward Primer GC%
        fwd_gcs = [c['PRIMER_LEFT_GC_PERCENT'] for c in all_candidates]
        if fwd_gcs:
            ax3 = self.figure.add_subplot(gs[1, 0])
            ax3.hist(fwd_gcs, bins=20, color='lightgreen', edgecolor='black')
            ax3.set_title('Fwd Primer GC%', fontsize=10)
            ax3.set_xlabel('GC%', fontsize=8)
            ax3.set_ylabel('Frequency', fontsize=8)
            ax3.tick_params(axis='both', which='major', labelsize=7)

        # Plot 4: Histogram of Cross-Dimer dG
        cross_dgs = [c['PRIMER_PAIR_CROSS_DIMER_DG'] for c in all_candidates]
        if cross_dgs:
            ax4 = self.figure.add_subplot(gs[1, 1])
            ax4.hist(cross_dgs, bins=20, color='gold', edgecolor='black')
            ax4.set_title('Cross-Dimer dG', fontsize=10)
            ax4.set_xlabel('dG (kcal/mol)', fontsize=8)
            ax4.set_ylabel('Frequency', fontsize=8)
            ax4.tick_params(axis='both', which='major', labelsize=7)
        
        self.figure.tight_layout()
        self.canvas.draw()

    def export_to_csv(self):
        if not self.results_data:
            self.show_error("No results to export.")
            return
        
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "primer_results.csv", "CSV Files (*.csv)")
        if path:
            try:
                with open(path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    headers = [self.results_table.horizontalHeaderItem(i).text() for i in range(self.results_table.columnCount())]
                    writer.writerow(headers)
                    for row in range(self.results_table.rowCount()):
                        writer.writerow([self.results_table.item(row, col).text() for col in range(self.results_table.columnCount())])
                self.show_info("Results exported to CSV successfully.")
            except Exception as e:
                self.show_error(f"Error exporting to CSV: {e}")

    def export_to_fasta(self):
        if not self.results_data:
            self.show_error("No results to export.")
            return

        path, _ = QFileDialog.getSaveFileName(self, "Save FASTA", "primer_results.fasta", "FASTA Files (*.fasta)")
        if path:
            try:
                with open(path, 'w') as f:
                    num_primers = self.results_data.get('PRIMER_PAIR_NUM_RETURNED', 0)
                    for i in range(num_primers):
                        f.write(f">Pair_{i+1}_Forward_Score_{self.results_data.get(f'PRIMER_PAIR_{i}_CUSTOM_SCORE', 0):.2f}\n")
                        f.write(f"{self.results_data[f'PRIMER_LEFT_{i}_SEQUENCE']}\n")
                        f.write(f">Pair_{i+1}_Reverse_Score_{self.results_data.get(f'PRIMER_PAIR_{i}_CUSTOM_SCORE', 0):.2f}\n")
                        f.write(f"{self.results_data[f'PRIMER_RIGHT_{i}_SEQUENCE']}\n")
                self.show_info("Results exported to FASTA successfully.")
            except Exception as e:
                self.show_error(f"Error exporting to FASTA: {e}")

    def show_error(self, message):
        QMessageBox.critical(self, "Error", message)
        self.run_button.setText("Design Primers")
        self.run_button.setEnabled(True)

    def show_info(self, message):
        QMessageBox.information(self, "Info", message)


if __name__ == '__main__':
    # Add a global exception hook to catch and display uncaught exceptions
    def except_hook(cls, exception, traceback):
        sys.__excepthook__(cls, exception, traceback)
    sys.excepthook = except_hook

    app = QApplication(sys.argv)
    
    app.setStyle('Fusion')
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.WindowText, Qt.white)
    palette.setColor(QPalette.Base, QColor(25, 25, 25))
    palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    palette.setColor(QPalette.ToolTipBase, Qt.white)
    palette.setColor(QPalette.ToolTipText, Qt.white)
    palette.setColor(QPalette.Text, Qt.white)
    palette.setColor(QPalette.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ButtonText, Qt.white)
    palette.setColor(QPalette.BrightText, Qt.red)
    palette.setColor(QPalette.Link, QColor(42, 130, 218))
    palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(palette)
    
    main_win = PrimerDesignApp()
    main_win.show()
    sys.exit(app.exec_())

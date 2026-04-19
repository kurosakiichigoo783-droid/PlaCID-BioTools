
import sys
import pytest
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import Qt
from docking_gui_complete import PrimerDesignApp

def test_primer_designer_e2e(qtbot):
    """
    End-to-end test for the primer designer application.
    """
    app = PrimerDesignApp()
    qtbot.addWidget(app)

    # 1. Enter sequence directly
    sequence_input = "ATTACAATTTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT"
    app.sequence_input.setText(sequence_input)

    # 2. Test NCBI fetch
    app.email_input.setText("gemini@example.com")
    app.gene_id_input.setText("NM_000546.6")

    # Click fetch button and wait
    qtbot.mouseClick(app.fetch_button, Qt.LeftButton)
    qtbot.wait(10000) # Wait for fetch to complete

    # Assert that the sequence is fetched
    assert len(app.sequence_input.toPlainText()) > 0

    # 3. After fetching, design primers
    qtbot.mouseClick(app.run_button, Qt.LeftButton)

    def worker_created():
        assert hasattr(app, 'worker')
    qtbot.waitUntil(worker_created)

    with qtbot.waitSignal(app.worker.finished, timeout=30000):
        pass # worker is already running


    # 4. Check for results
    assert app.results_table.rowCount() > 0

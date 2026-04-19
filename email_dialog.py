
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

# Medical Records Digitizer

A desktop application designed to scan, digitize, and review laboratory follow-up forms into structured, per-patient records. Built using Python with a Tkinter graphical interface, it leverages Advanced Image Processing (OpenCV) and PaddleOCR to streamline data entry for healthcare metrics.

---

## Features

* **Advanced Document Preprocessing**: 
  * Perspective-correction ("flat scan" effect) using automatic contour boundary detection.
  * Adaptive uneven lighting and shadow removal.
  * Contrast optimization via CLAHE (Contrast Limited Adaptive Histogram Equalization).
  * Automated orientation check (rotates landscape photos back to portrait format).
* **Robust Text Extraction & Alignment**:
  * Powered by PaddleOCR for high-precision text detection.
  * Intelligently groups fragmented word blocks horizontally to handle dot leaders (e.g., `HEMATIES ..... 4.90`) and extracts accurate values, labels, and units.
* **Smart Vocabulary & Synonym Learning**:
  * Map regional lab terminology or custom acronyms to standard diagnostic categories.
  * Automatically learns new lab shorthand variants directly from review corrections.
* **Interactive Data Grid**: 
  * Real-time range checking and outlier identification for specialized metrics (HbA1c, Glycemia, BMI, Blood Pressure, etc.).
  * Direct manual field creation alongside existing data loading and comparison capabilities.
* **Flexible Exports**: Writes organized structured profiles downstream into standalone, per-patient JSON profiles as well as an inclusive master `records.csv` compilation.

---

## Getting Started

### Screenshot

![screenshot 1](screenshots/screenshot_1.png)

![screenshot 1](screenshots/screenshot_2.png)

### Installation

1. Clone the repository and navigate into it:
   ```bash
   git clone https://github.com/abdoulahmr/Medical-Report-Scanner
   cd Medical-Report-Scanner

2. Create virtual environment and activate it:
   ```bash
   python -m venv venv
   source venv/bin/activate (macOS / Linux)
   .\venv\Scripts\activate (windows)

3. install required libraries and run:
   ```bash
   pip install -r requirements.txt
   python main.py
import os
import re
import csv
import json
import glob
import queue
import string
import hashlib
import difflib
import threading
import unicodedata
from datetime import datetime

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import cv2
from PIL import Image, ImageTk  

# 1. PaddleOCR (Completely separate, fails gracefully)
try:
    from paddleocr import PaddleOCR
except ImportError:
    PaddleOCR = None

# 2. Matplotlib (Separate so it doesn't break if sklearn is missing)
try:
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# 3. Scikit-Learn
try:
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer, KNNImputer
    
    # MICE is experimental in sklearn, so it requires this enable flag first
    from sklearn.experimental import enable_iterative_imputer  
    from sklearn.impute import IterativeImputer
    
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

# 4. UMAP
try:
    from umap import UMAP
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

IMAGE_EXTENSIONS = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff", ".webp"}

KNOWN_FIELDS = {
    "glycemie a jeun": {"unit": "g/L", "min": 0.3, "max": 6.0},
    "glycemie": {"unit": "g/L", "min": 0.3, "max": 6.0},
    "glycemie post prandiale": {"unit": "g/L", "min": 0.3, "max": 6.0},
    "glycemie postprandiale": {"unit": "g/L", "min": 0.3, "max": 6.0},
    "hba1c": {"unit": "%", "min": 3.0, "max": 18.0},
    "hemoglobine glyquee": {"unit": "%", "min": 3.0, "max": 18.0},
    "poids": {"unit": "kg", "min": 25.0, "max": 220.0},
    "taille": {"unit": "cm", "min": 100.0, "max": 210.0},
    "imc": {"unit": "kg/m2", "min": 12.0, "max": 60.0},
    "tension arterielle systolique": {"unit": "mmHg", "min": 70.0, "max": 260.0},
    "tension arterielle diastolique": {"unit": "mmHg", "min": 40.0, "max": 150.0},
    "creatinine": {"unit": "mg/L", "min": 2.0, "max": 100.0},
    "cholesterol total": {"unit": "g/L", "min": 0.5, "max": 5.0},
    "triglycerides": {"unit": "g/L", "min": 0.2, "max": 10.0},
}

LABEL_SUGGESTIONS = [
    # Hémogramme (NFS) - Lignée Rouge
    "Hématies", 
    "Hémoglobine", 
    "Hématocrite", 
    "V.G.M.", 
    "T.C.M.H.", 
    "C.C.M.H.", 
    "I.D.R. (R.D.W.)", 

    # Hémogramme (NFS) - Lignée Blanche
    "Leucocytes", 
    "Polynucléaires Neutrophiles %", 
    "Polynucléaires Neutrophiles (mm3)", 
    "Lymphocytes %", 
    "Lymphocytes (mm3)", 
    "Monocytes %", 
    "Monocytes (mm3)", 
    "Eosinophiles %", 
    "Eosinophiles (mm3)", 
    "Polynucléaires Basophiles %", 
    "Polynucléaires Basophiles (mm3)", 
    "Cellules de taille moyenne (MID) %",
    "Cellules de taille moyenne (MID) (mm3)",
    "Granulocytes %",
    "Granulocytes (mm3)",

    # Hémogramme (NFS) - Lignée Plaquettaire
    "Plaquettes", 
    "V.P.M. (MPV)", 
    "P.D.W.", 
    "P-LCR", 

    # Métabolisme des Glucides
    "Glycémie à jeun (g/l)", 
    "Glycémie à jeun (mmol/l)", 
    "Glycémie post-prandiale",
    "Hémoglobine glyquée (HbA1c)", 

    # Bilan Rénal & Urinaire
    "Urée sanguine", 
    "Créatinine sanguine", 
    "Volume urinaire des 24 heures", 
    "Microalbuminurie (mg/l)", 
    "Microalbuminurie (mg/24h)", 
    "Protéinurie des 24 heures",

    # Bilan Lipidique
    "Cholestérol total (g/l)", 
    "Cholestérol total (mmol/l)", 
    "Triglycérides (g/l)", 
    "Triglycérides (mmol/l)", 
    "Cholestérol HDL", 
    "Cholestérol LDL", 
    "Rapport Cholestérol Total / HDL",

    # Bilan Hépatique
    "Transaminases SGOT (ASAT)", 
    "Transaminases SGPT (ALAT)", 

    # Endocrinologie & Métabolisme
    "T.S.H. ultra-sensible", 
    "Calcémie",
    "Fer sérique",
    "Ferritine",

    # Inflammation
    "Protéine C-réactive (CRP)",

    # Pourcentage de cellules de taille moyenne
    "Mid %",
    "Mid (mm3)",
]

BUILTIN_ALIASES = {
    "HEMATIES": ["GR", "Hematies", "RBC"],
    "Hémoglobine": ["Hb", "HGB", "Hb g/dl"],
    "Hematocrite": ["Hte", "HCT", "Ht"],
    "V.G.M": ["VGM", "MCV"],
    "T.C.M.H": ["TCMH", "MCH"],
    "C.C.M.H.": ["CCMH", "MCHC"],
    "I.D.R(R.D.W)": ["IDR", "RDW"],
    "LEUCOCYTES": ["GB", "WBC", "Leuco"],
    "Polynucléaires Neutrophiles %": ["PNN %", "Neutro %", "Neut %"],
    "Polynucléaires Neutrophiles mm3": ["PNN mm3", "Neutro mm3", "Neut mm3"],
    "Lymphocytes %": ["Lympho %", "Lym %"],
    "Lymphocytes mm3": ["Lympho mm3", "Lym mm3"],
    "Monocytes %": ["Mono %"],
    "Monocytes mm3": ["Mono mm3"],
    "Eosinophiles %": ["Eosino %", "Eo %"],
    "Eosinophiles mm3": ["Eosino mm3", "Eo mm3"],
    "Polynucléaires Basophiles %": ["PNB %", "Baso %"],
    "Polynucléaires Basophiles mm3": ["PNB mm3", "Baso mm3"],
    "Glycime a jeun g/l": ["Gly", "Glyc", "GAJ", "Glycemie"],
    "Glycime a jeun mmol/l": ["Gly mmol", "GAJ mmol"],
    "Uree sanguine": ["Uree", "BUN"],
    "Creatinine sanguine": ["Creat", "Crea", "Creatinine"],
    "Transaminases SGOT": ["SGOT", "ASAT", "TGO"],
    "Transaminases SGPT": ["SGPT", "ALAT", "TGP"],
    "Hemoglobine glyquee HbA1c": ["HbA1c", "Hb glyquee", "Hb A1c"],
    "PLAQUETTES": ["Plaq", "PLT"],
    "CHOLESTEROLTOTAL g/l": ["Chol T", "CT", "Cholesterol total"],
    "TRIGLYCERIDES g/l": ["TG", "Triglycerides"],
    "CHOLESTEROL HDL": ["HDL", "HDL-C", "HDLc"],
    "CHOLESTEROL传统 LDL": ["LDL", "LDL-C", "LDLc"],
    "MICROALBUMINURIE mg/l": ["MAU", "Microalb"],
    "T.S.H. ultra-sensible": ["TSH", "TSHus"],
}

DATE_PATTERNS = [
    re.compile(r"\b(\d{1,2})[/\-.](\d{1,2})[/\-.](\d{2,4})\b"),
]


def normalize_label(label):
    label = label.strip().lower()
    label = unicodedata.normalize("NFKD", label)
    label = "".join(c for c in label if not unicodedata.combining(c))
    label = re.sub(r"\s+", " ", label)
    return label


def try_parse_date(text):
    for pat in DATE_PATTERNS:
        m = pat.search(text)
        if m:
            d, mth, y = m.groups()
            try:
                d, mth, y = int(d), int(mth), int(y)
                if y < 100:
                    y += 2000 if y < 50 else 1900
                dt = datetime(y, mth, d)
                return dt.strftime("%Y-%m-%d")
            except ValueError:
                continue
    return None


DOT_LEADER_RE = re.compile(r"(?:\.\s*){2,}")
NUMBER_RE = re.compile(r"-?\d+(?:[.,]\d+)?")

_LABEL_LOOKUP_NORM = {}
for _suggestion in LABEL_SUGGESTIONS:
    _LABEL_LOOKUP_NORM[normalize_label(_suggestion)] = _suggestion

_BUILTIN_ALIAS_LOOKUP_NORM = {}
for _canonical, _synonyms in BUILTIN_ALIASES.items():
    for _syn in _synonyms:
        _BUILTIN_ALIAS_LOOKUP_NORM[normalize_label(_syn)] = _canonical


def group_ocr_lines(raw_boxes):
    if not raw_boxes:
        return []

    items = []
    for box in raw_boxes:
        xs = [p[0] for p in box["bbox"]]
        ys = [p[1] for p in box["bbox"]]
        items.append({
            "text": box["text"],
            "confidence": box["confidence"],
            "x0": min(xs), "x1": max(xs),
            "y0": min(ys), "y1": max(ys),
            "cy": (min(ys) + max(ys)) / 2.0,
            "h": max(1.0, max(ys) - min(ys)),
        })
    items.sort(key=lambda it: it["cy"])

    rows = []  
    for it in items:
        placed = False
        for row in rows:
            row_cy = row["cy_sum"] / row["count"]
            row_h = row["h_sum"] / row["count"]
            tolerance = max(row_h, it["h"]) * 0.7
            if abs(it["cy"] - row_cy) <= tolerance:
                row["items"].append(it)
                row["cy_sum"] += it["cy"]
                row["h_sum"] += it["h"]
                row["count"] += 1
                placed = True
                break
        if not placed:
            rows.append({"items": [it], "cy_sum": it["cy"], "h_sum": it["h"], "count": 1})

    grouped = []
    for row in rows:
        row_items = sorted(row["items"], key=lambda it: it["x0"])
        text = " ".join(it["text"] for it in row_items)
        text = DOT_LEADER_RE.sub(" ", text)
        text = re.sub(r"\s+", " ", text).strip()
        if not text:
            continue
        x0 = min(it["x0"] for it in row_items)
        x1 = max(it["x1"] for it in row_items)
        y0 = min(it["y0"] for it in row_items)
        y1 = max(it["y1"] for it in row_items)
        conf = sum(it["confidence"] for it in row_items) / len(row_items)
        grouped.append({
            "text": text,
            "confidence": conf,
            "bbox": [[x0, y0], [x1, y0], [x1, y1], [x0, y1]],
        })

    grouped.sort(key=lambda g: g["bbox"][0][1])
    return grouped


def parse_label_value_unit(text):
    cleaned = DOT_LEADER_RE.sub(" ", text)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    m = NUMBER_RE.search(cleaned)
    if not m:
        return cleaned, "", ""
    label = cleaned[:m.start()].strip(" .:-\u2013\u2014")
    value = m.group(0)
    unit = cleaned[m.end():].strip(" .:-\u2013\u2014")
    return label, value, unit


def match_known_label(raw_label, extra_aliases=None):
    if not raw_label:
        return raw_label
    norm = normalize_label(raw_label)
    if not norm:
        return raw_label

    if norm in _LABEL_LOOKUP_NORM:
        return _LABEL_LOOKUP_NORM[norm]

    if norm in _BUILTIN_ALIAS_LOOKUP_NORM:
        return _BUILTIN_ALIAS_LOOKUP_NORM[norm]

    if extra_aliases and norm in extra_aliases:
        return extra_aliases[norm]["canonical"]

    best = difflib.get_close_matches(norm, _LABEL_LOOKUP_NORM.keys(), n=1, cutoff=0.7)
    if best:
        return _LABEL_LOOKUP_NORM[best[0]]

    return raw_label


def validate_value(label, value_text):
    key = normalize_label(label)
    spec = KNOWN_FIELDS.get(key)
    if not spec:
        return None
    num_match = re.search(r"-?\d+[.,]?\d*", value_text)
    if not num_match:
        return "expected a numeric value"
    try:
        val = float(num_match.group(0).replace(",", "."))
    except ValueError:
        return None
    if val < spec["min"] or val > spec["max"]:
        return f"unusual value for {label} (typical range {spec['min']}-{spec['max']} {spec['unit']})"
    return None


def file_hash(path):
    st = os.stat(path)
    key = f"{os.path.abspath(path)}|{st.st_size}|{int(st.st_mtime)}"
    return hashlib.sha1(key.encode("utf-8")).hexdigest()


def order_points(pts):
    rect = np.zeros((4, 2), dtype="float32")
    s = pts.sum(axis=1)
    rect[0] = pts[np.argmin(s)]
    rect[2] = pts[np.argmax(s)]
    diff = np.diff(pts, axis=1)
    rect[1] = pts[np.argmin(diff)]
    rect[3] = pts[np.argmax(diff)]
    return rect


def find_document_contour(img):
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    edged = cv2.Canny(blurred, 50, 150)
    edged = cv2.dilate(edged, np.ones((5, 5), np.uint8), iterations=2)
    edged = cv2.erode(edged, np.ones((5, 5), np.uint8), iterations=1)

    contours, _ = cv2.findContours(edged, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None
    contours = sorted(contours, key=cv2.contourArea, reverse=True)[:10]

    img_area = img.shape[0] * img.shape[1]
    for c in contours:
        peri = cv2.arcLength(c, True)
        approx = cv2.approxPolyDP(c, 0.02 * peri, True)
        if len(approx) == 4 and cv2.contourArea(approx) > 0.2 * img_area:
            return approx.reshape(4, 2).astype("float32")
    return None


def four_point_transform(img, pts):
    rect = order_points(pts)
    (tl, tr, br, bl) = rect

    width_a = np.linalg.norm(br - bl)
    width_b = np.linalg.norm(tr - tl)
    max_width = max(int(width_a), int(width_b), 1)

    height_a = np.linalg.norm(tr - br)
    height_b = np.linalg.norm(tl - bl)
    max_height = max(int(height_a), int(height_b), 1)

    dst = np.array([
        [0, 0],
        [max_width - 1, 0],
        [max_width - 1, max_height - 1],
        [0, max_height - 1],
    ], dtype="float32")

    m = cv2.getPerspectiveTransform(rect, dst)
    return cv2.warpPerspective(img, m, (max_width, max_height))


def compute_skew_angle(img, max_correction=15.0):
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    gray = cv2.bitwise_not(gray)
    thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)[1]

    coords = np.column_stack(np.where(thresh > 0))
    if coords.shape[0] < 20:
        return 0.0

    angle = cv2.minAreaRect(coords)[-1]
    if angle < -45:
        angle = -(90 + angle)
    else:
        angle = -angle

    if abs(angle) > max_correction:
        return 0.0
    return angle


def remove_shadows(img):
    planes = cv2.split(img)
    result_planes = []
    for plane in planes:
        dilated = cv2.dilate(plane, np.ones((7, 7), np.uint8))
        bg = cv2.medianBlur(dilated, 21)
        diff = 255 - cv2.absdiff(plane, bg)
        norm = cv2.normalize(diff, None, alpha=0, beta=255,
                              norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_8UC1)
        result_planes.append(norm)
    return cv2.merge(result_planes)


class ImagePipeline:

    def __init__(self, path):
        self.path = path
        self.original = cv2.imread(path, cv2.IMREAD_COLOR)
        if self.original is None:
            raise ValueError(f"Could not read image: {path}")
        self.rotation_90 = 0
        self.fine_angle = 0.0
        self.crop_rect = None
        self.denoise = False
        self.contrast = False

        self.scanned_base = None
        self.shadow_removal = False

    def reset(self):
        self.rotation_90 = 0
        self.fine_angle = 0.0
        self.crop_rect = None
        self.denoise = False
        self.contrast = False
        self.scanned_base = None
        self.shadow_removal = False

    def auto_scan(self):
        self.crop_rect = None
        self.rotation_90 = 0
        self.fine_angle = 0.0

        working = self.original.copy()

        quad = find_document_contour(working)
        found = quad is not None
        scanned = four_point_transform(working, quad) if found else working

        h, w = scanned.shape[:2]
        if w > h * 1.15:
            scanned = cv2.rotate(scanned, cv2.ROTATE_90_CLOCKWISE)

        skew = compute_skew_angle(scanned)
        if abs(skew) > 0.3:
            hh, ww = scanned.shape[:2]
            center = (ww / 2, hh / 2)
            m = cv2.getRotationMatrix2D(center, skew, 1.0)
            scanned = cv2.warpAffine(scanned, m, (ww, hh),
                                      flags=cv2.INTER_LINEAR,
                                      borderMode=cv2.BORDER_REPLICATE)

        self.scanned_base = scanned
        self.shadow_removal = True
        self.contrast = True
        return found

    def render(self):
        img = (self.scanned_base if self.scanned_base is not None else self.original).copy()

        if self.crop_rect:
            x1, y1, x2, y2 = self.crop_rect
            x1, x2 = sorted((max(0, x1), max(0, x2)))
            y1, y2 = sorted((max(0, y1), max(0, y2)))
            h, w = img.shape[:2]
            x2, y2 = min(w, x2), min(h, y2)
            if x2 - x1 > 5 and y2 - y1 > 5:
                img = img[y1:y2, x1:x2]

        if self.rotation_90 == 90:
            img = cv2.rotate(img, cv2.ROTATE_90_CLOCKWISE)
        elif self.rotation_90 == 180:
            img = cv2.rotate(img, cv2.ROTATE_180)
        elif self.rotation_90 == 270:
            img = cv2.rotate(img, cv2.ROTATE_90_COUNTERCLOCKWISE)

        if abs(self.fine_angle) > 0.01:
            h, w = img.shape[:2]
            center = (w / 2, h / 2)
            m = cv2.getRotationMatrix2D(center, self.fine_angle, 1.0)
            img = cv2.warpAffine(img, m, (w, h),
                                  flags=cv2.INTER_LINEAR,
                                  borderMode=cv2.BORDER_REPLICATE)

        if self.shadow_removal:
            if len(img.shape) == 2:
                img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
            elif img.shape[2] == 4:
                img = cv2.cvtColor(img, cv2.COLOR_BGRA2BGR)
            img = remove_shadows(img)

        if self.denoise:
            if len(img.shape) == 2:
                img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
            elif img.shape[2] == 4:
                img = cv2.cvtColor(img, cv2.COLOR_BGRA2BGR)
            img = cv2.fastNlMeansDenoisingColored(
                src=img,
                dst=None,
                h=7,
                hColor=7,
                templateWindowSize=7,
                searchWindowSize=21
            )

        if self.contrast:
            lab = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)
            l, a, b = cv2.split(lab)
            clahe = cv2.createCLAHE(clipLimit=2.5, tileGridSize=(8, 8))
            l = clahe.apply(l)
            lab = cv2.merge((l, a, b))
            img = cv2.cvtColor(lab, cv2.COLOR_LAB2BGR)

        return img


class OCREngine:

    def __init__(self):
        self._ocr = None
        self._lock = threading.Lock()

    @staticmethod
    def _paddleocr_cache_dir():
        return os.path.join(os.path.expanduser("~"), ".paddleocr")

    def _init_paddleocr(self):
        return PaddleOCR(use_angle_cls=True, lang="fr", show_log=False)

    def ensure_loaded(self):
        with self._lock:
            if self._ocr is None:
                if PaddleOCR is None:
                    raise RuntimeError(
                        "PaddleOCR is not installed. Run:\n"
                        "  pip install paddleocr paddlepaddle"
                    )
                try:
                    self._ocr = self._init_paddleocr()
                except Exception as e:
                    msg = str(e).lower()
                    corruption_hints = ("unexpected end of data", "tarfile",
                                         "eof", "not a gzip file", "truncat")
                    cache_dir = self._paddleocr_cache_dir()
                    if any(h in msg for h in corruption_hints) and os.path.isdir(cache_dir):
                        import shutil
                        shutil.rmtree(cache_dir, ignore_errors=True)
                        try:
                            self._ocr = self._init_paddleocr()
                        except Exception as e2:
                            raise RuntimeError(
                                "PaddleOCR model files looked corrupted, so the "
                                f"cache at {cache_dir} was cleared and a fresh "
                                "download was attempted, but it still failed:\n"
                                f"{e2}\n\n"
                                "Check your internet connection and try Run OCR again."
                            ) from e2
                    else:
                        raise
        return self._ocr

    def run(self, bgr_image):
        ocr = self.ensure_loaded()
        result = ocr.ocr(bgr_image, cls=True)
        lines = []
        if not result:
            return lines
        page = result[0] if isinstance(result[0], list) or result[0] is None else result
        page = page or []
        for entry in page:
            try:
                bbox, (text, conf) = entry
            except (ValueError, TypeError):
                continue
            lines.append({
                "bbox": [[float(x), float(y)] for x, y in bbox],
                "text": text,
                "confidence": float(conf),
            })
        return lines


class DataStore:
    def __init__(self, folder):
        self.folder = folder
        self.records_path = os.path.join(folder, "records.json")
        self.manifest_path = os.path.join(folder, "manifest.json")
        self.csv_path = os.path.join(folder, "records.csv")
        self.patients_dir = os.path.join(folder, "patients")
        self.aliases_path = os.path.join(folder, "label_aliases.json")
        self.data = {"patients": {}}
        self.manifest = {"processed": {}}
        self.aliases = {}   
        self._load()

    ###################################
    # Developer / maintenance methods #
    ###################################
    def get_all_variable_names(self):
        """Return a sorted list of tuples: [(label, count), ...]"""
        from collections import Counter
        counts = Counter()
        for p in self.data.get("patients", {}).values():
            for visit in p.get("visits", []):
                for f in visit.get("fields", []):
                    lbl = f.get("label", "").strip()
                    if lbl:
                        counts[lbl] += 1
        # Sort by count (descending), then alphabetically
        return sorted(counts.items(), key=lambda x: (-x[1], x[0]))
    ###################################

    def _load(self):
        if os.path.exists(self.records_path):
            with open(self.records_path, "r", encoding="utf-8") as f:
                self.data = json.load(f)
        if os.path.exists(self.manifest_path):
            with open(self.manifest_path, "r", encoding="utf-8") as f:
                self.manifest = json.load(f)
        if os.path.exists(self.aliases_path):
            with open(self.aliases_path, "r", encoding="utf-8") as f:
                self.aliases = json.load(f)
        self._migrate()

    def _migrate(self):
        for p in self.data.get("patients", {}).values():
            for visit in p.get("visits", []):
                if "images" not in visit:
                    old = visit.pop("source_image", None)
                    visit["images"] = [old] if old else []

    def learn_alias(self, raw_text, canonical):
        norm = normalize_label(raw_text)
        if not norm or normalize_label(canonical) == norm:
            return False
        existing = self.aliases.get(norm)
        if existing and existing.get("canonical") == canonical:
            return False
        self.aliases[norm] = {"canonical": canonical, "example": raw_text}
        return True

    def remove_alias(self, norm_key):
        self.aliases.pop(norm_key, None)

    def save_aliases(self):
        os.makedirs(self.folder, exist_ok=True)
        with open(self.aliases_path, "w", encoding="utf-8") as f:
            json.dump(self.aliases, f, ensure_ascii=False, indent=2)

    def standardize_labels(self, mapping):
        """
        mapping = {
            "CCMH": "MCHC",
            "C.C.M.H.": "MCHC",
            "TCMH": "MCH",
            ...
        }

        Returns number of modified fields.
        """

        normalized = {
            normalize_label(k): v
            for k, v in mapping.items()
        }

        count = 0

        for patient in self.data["patients"].values():
            for visit in patient["visits"]:
                for field in visit["fields"]:

                    current = field.get("label", "")
                    key = normalize_label(current)

                    if key in normalized:
                        new_name = normalized[key]

                        if current != new_name:
                            field["label"] = new_name
                            count += 1

        self.save()

        return count
        
    def save(self):
        os.makedirs(self.folder, exist_ok=True)
        os.makedirs(self.patients_dir, exist_ok=True)
        with open(self.records_path, "w", encoding="utf-8") as f:
            json.dump(self.data, f, ensure_ascii=False, indent=2)
        with open(self.manifest_path, "w", encoding="utf-8") as f:
            json.dump(self.manifest, f, ensure_ascii=False, indent=2)
        self.save_aliases()
        self._export_csv()
        self._export_patient_json()

    def is_processed(self, path):
        return file_hash(path) in self.manifest["processed"]

    def find_by_image(self, path):
        entry = self.manifest["processed"].get(file_hash(path))
        if not entry:
            return None, None
        return entry["patient_id"], entry["visit_index"]

    def patient_ids(self):
        return sorted(self.data["patients"].keys())

    def patient_name(self, pid):
        p = self.data["patients"].get(pid)
        return p["name"] if p else ""

    def ensure_patient(self, pid, name):
        if pid not in self.data["patients"]:
            self.data["patients"][pid] = {"name": name, "visits": []}
        elif name:
            self.data["patients"][pid]["name"] = name

    def visits(self, pid):
        p = self.data["patients"].get(pid)
        return p["visits"] if p else []

    def record_options(self, pid):
        options = []
        for idx, visit in enumerate(self.visits(pid)):
            date = visit.get("visit_date") or "no date"
            n_images = len(visit.get("images", []))
            n_fields = len(visit.get("fields", []))
            options.append((idx, f"{date}  -  {n_images} image(s), {n_fields} field(s)"))
        options.sort(key=lambda t: self.visits(pid)[t[0]].get("visit_date") or "")
        return options

    def get_visit(self, pid, visit_index):
        visits = self.visits(pid)
        if visit_index is None or not (0 <= visit_index < len(visits)):
            return None
        return visits[visit_index]

    def save_or_update_visit(self, pid, name, visit_index, image_path, visit_date, fields):
        self.ensure_patient(pid, name)
        visits = self.data["patients"][pid]["visits"]
        abs_path = os.path.abspath(image_path)

        if visit_index is None:
            visit = {
                "images": [abs_path],
                "visit_date": visit_date,
                "processed_at": datetime.now().isoformat(timespec="seconds"),
                "fields": fields,
            }
            visits.append(visit)
            visit_index = len(visits) - 1
        else:
            visit = visits[visit_index]
            if abs_path not in visit["images"]:
                visit["images"].append(abs_path)
            if visit_date:
                visit["visit_date"] = visit_date
            visit["fields"] = fields
            visit["processed_at"] = datetime.now().isoformat(timespec="seconds")

        h = file_hash(image_path)
        self.manifest["processed"][h] = {
            "path": abs_path,
            "patient_id": pid,
            "visit_index": visit_index,
        }
        return visit_index

    def _export_csv(self):
        with open(self.csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["patient_id", "patient_name", "visit_date",
                        "source_images", "label", "value", "confidence",
                        "unit", "flag"])
            for pid, p in sorted(self.data["patients"].items()):
                for visit in p["visits"]:
                    images_joined = "; ".join(os.path.basename(i) for i in visit.get("images", []))
                    for field in visit["fields"]:
                        w.writerow([
                            pid, p["name"], visit["visit_date"],
                            images_joined,
                            field.get("label", ""),
                            field.get("value", ""),
                            field.get("confidence", ""),
                            field.get("unit", ""),
                            field.get("flag", ""),
                        ])

    def _export_patient_json(self):
        for pid, p in self.data["patients"].items():
            path = os.path.join(self.patients_dir, f"{pid}.json")
            visits = sorted(p["visits"], key=lambda v: v.get("visit_date") or "")
            payload = {
                "patient_id": pid,
                "name": p.get("name", ""),
                "visits": [
                    {
                        "visit_date": v.get("visit_date"),
                        "images": [os.path.basename(i) for i in v.get("images", [])],
                        "processed_at": v.get("processed_at"),
                        "fields": v.get("fields", []),
                    }
                    for v in visits
                ],
            }
            with open(path, "w", encoding="utf-8") as f:
                json.dump(payload, f, ensure_ascii=False, indent=2)

class MultivariateSetupDialog(tk.Toplevel):
    """Dialog to select variables and imputation method before running multivariate plots."""

    def __init__(self, parent, variables_with_counts, visit_coverage):
        super().__init__(parent)
        self.title("Multivariate Analysis Setup")
        self.geometry("600x680")
        self.resizable(True, True)
        self.transient(parent)
        self.grab_set()
        
        self.result = None
        self._visit_coverage = visit_coverage  # Dict mapping (pid, v_idx) -> set(labels)

        self._build_ui(variables_with_counts)
        self.protocol("WM_DELETE_WINDOW", self.on_cancel)

        # Center the dialog over the parent window
        self.update_idletasks()
        x = parent.winfo_x() + (parent.winfo_width() // 2) - (self.winfo_width() // 2)
        y = parent.winfo_y() + (parent.winfo_height() // 2) - (self.winfo_height() // 2)
        self.geometry(f"+{x}+{y}")
        
        # Trigger initial calculation
        self._update_status()

    def _auto_detect(self):
        """Automatically find the best combination of variables and imputation method."""
        # Get all variables and sort by occurrence count (ascending) 
        # so we can drop the rarest variables first if needed.
        all_items = []
        for iid in self.tree.get_children():
            var_name, count_str = self.tree.item(iid, "values")
            all_items.append((var_name, int(count_str), iid))
            
        all_items.sort(key=lambda x: x[1]) # Sort rarest to most common

        min_vars = 2
        min_samples = 3
        best_combo_iids = []
        best_method = "Drop missing visits (Complete Cases)"

        # Strategy 1: Try to use Complete Cases (Drop missing) by progressively 
        # removing the rarest variables until we have enough valid visits.
        for i in range(len(all_items), min_vars - 1, -1):
            # Take the 'i' most common variables
            current_vars = [item[0] for item in all_items[-i:]]
            current_iids = [item[2] for item in all_items[-i:]]
            
            valid_samples = 0
            for key, vars_in_visit in self._visit_coverage.items():
                if set(current_vars).issubset(vars_in_visit):
                    valid_samples += 1
            
            if valid_samples >= min_samples:
                best_combo_iids = current_iids
                best_method = "Drop missing visits (Complete Cases)"
                break  # Found the maximum possible variables for Complete Cases

        # Strategy 2: If even the top 2 most common variables don't share 3 visits,
        # fallback to MICE imputation using ALL variables.
        if not best_combo_iids:
            best_combo_iids = [item[2] for item in all_items]
            best_method = "Multiple Imputation by Chained Equations (MICE)"
            
            # Verify if we even have enough visits with AT LEAST ONE variable
            valid_samples = sum(1 for vars_in_visit in self._visit_coverage.values() if vars_in_visit)
            if valid_samples < min_samples:
                messagebox.showwarning(
                    "Auto Detect Failed", 
                    "There is not enough overlapping data in the entire dataset to run multivariate plots.", 
                    parent=self
                )
                return

        # Apply the found settings to the UI
        self.tree.selection_remove(self.tree.selection())
        self.tree.selection_add(best_combo_iids)
        self.impute_var.set(best_method)
        
        # Force status update to show the new green checkmarks
        self._update_status()

    def _build_ui(self, variables_with_counts):
        main_frame = ttk.Frame(self, padding=15)
        main_frame.pack(fill=tk.BOTH, expand=True)

        ttk.Label(main_frame, text="Select variables to include in the analysis:", font=("", 10, "bold")).pack(anchor="w", pady=(0, 5))

        # --- Occurrence Table ---
        tree_frame = ttk.Frame(main_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        cols = ("variable", "count")
        self.tree = ttk.Treeview(tree_frame, columns=cols, show="headings", selectmode="extended", height=9)
        self.tree.heading("variable", text="Variable Name")
        self.tree.heading("count", text="Occurrences")
        self.tree.column("variable", width=350, anchor="w")
        self.tree.column("count", width=100, anchor="center")

        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

        for var, count in variables_with_counts:
            iid = self.tree.insert("", tk.END, values=(var, count))
            self.tree.selection_add(iid)

        # Bind tree selection to update the status dynamically
        self.tree.bind("<<TreeviewSelect>>", lambda e: self._update_status())

        # --- Imputation Options ---
        opt_frame = ttk.LabelFrame(main_frame, text="Missing Value Handling")
        opt_frame.pack(fill=tk.X, pady=(0, 10))

        self.impute_var = tk.StringVar(value="Drop missing visits (Complete Cases)")
        methods = [
            "Drop missing visits (Complete Cases)",
            "Median Replacement",
            "Mean Replacement",
            "K-Nearest Neighbors (KNN)",
            "Multiple Imputation by Chained Equations (MICE)"
        ]
        for m in methods:
            ttk.Radiobutton(opt_frame, text=m, variable=self.impute_var, value=m).pack(anchor="w", padx=10, pady=1)
            
        # Bind radio buttons to update the status dynamically
        self.impute_var.trace_add("write", lambda *args: self._update_status())

        # --- Live Requirement Check Widget ---
        status_frame = ttk.LabelFrame(main_frame, text="Algorithm Data Requirements")
        status_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.status_text = tk.Text(status_frame, height=6, state="disabled", wrap="word", 
                                   font=("Consolas", 9), bg="#f9f9f9", relief="flat", bd=0)
        self.status_text.pack(fill=tk.X, padx=5, pady=5)

        # --- Buttons ---
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill=tk.X)
        ttk.Button(btn_frame, text="Cancel", command=self.on_cancel, width=15).pack(side=tk.RIGHT, padx=(5, 0))
        ttk.Button(btn_frame, text="Run Analysis", command=self.on_ok, width=15).pack(side=tk.RIGHT)
        ttk.Button(btn_frame, text="Auto Detect Best Combo", command=self._auto_detect).pack(side=tk.LEFT)

    def _update_status(self, *args):
        """Dynamically calculate if selected variables are enough for each plot type."""
        selected_vars = set([self.tree.item(iid, "values")[0] for iid in self.tree.selection()])
        impute_method = self.impute_var.get()
        n_vars = len(selected_vars)
        
        valid_samples = 0
        for key, vars_in_visit in self._visit_coverage.items():
            if impute_method == "Drop missing visits (Complete Cases)":
                # A visit is only valid if it has ALL selected variables
                if selected_vars.issubset(vars_in_visit):
                    valid_samples += 1
            else:
                # If imputing, a visit is valid if it has AT LEAST ONE selected variable
                if selected_vars.intersection(vars_in_visit):
                    valid_samples += 1

        # Define requirements: (min_vars, min_samples)
        requirements = {
            "Correlation Matrix": (2, 3),
            "PCA (2D)": (2, 3),
            "t-SNE (2D)": (2, 3),
            "UMAP (2D)": (2, 3)
        }
        
        self.status_text.config(state="normal")
        self.status_text.delete("1.0", tk.END)
        
        self.status_text.insert(tk.END, f"Variables selected: {n_vars}  |  Usable visits: {valid_samples}\n")
        self.status_text.insert(tk.END, "-" * 45 + "\n")
        
        for name, (req_vars, req_samples) in requirements.items():
            has_vars = n_vars >= req_vars
            has_samples = valid_samples >= req_samples
            
            if has_vars and has_samples:
                self.status_text.insert(tk.END, f"✔ {name}: OK\n")
            else:
                missing = []
                if not has_vars: missing.append(f"need {req_vars}+ vars")
                if not has_samples: missing.append(f"need {req_samples}+ visits")
                self.status_text.insert(tk.END, f"✖ {name}: NOT ENOUGH DATA ({', '.join(missing)})\n")
            
        self.status_text.config(state="disabled")

    def on_ok(self):
        selected_iids = self.tree.selection()
        if not selected_iids:
            messagebox.showwarning("Selection Required", "Please select at least one variable.", parent=self)
            return

        selected_vars = [self.tree.item(iid, "values")[0] for iid in selected_iids]
        impute_method = self.impute_var.get()
        
        self.result = (selected_vars, impute_method)
        self.destroy()

    def on_cancel(self):
        self.result = None
        self.destroy()

    def on_ok(self):
        selected_iids = self.tree.selection()
        if not selected_iids:
            messagebox.showwarning("Selection Required", "Please select at least one variable.", parent=self)
            return

        selected_vars = [self.tree.item(iid, "values")[0] for iid in selected_iids]
        impute_method = self.impute_var.get()
        
        self.result = (selected_vars, impute_method)
        self.destroy()

    def on_cancel(self):
        self.result = None
        self.destroy()

class EDAWindow(tk.Toplevel):
    """Exploratory Data Analysis window with Univariate Carousel and Multivariate Plots."""

    PLOT_TYPES = [
        ("Time Series (line)", "timeseries"),
        ("Histogram", "histogram"),
        ("Box Plot (per patient)", "boxplot"),
        ("Bar Chart (mean per patient)", "bar"),
        ("--- Multivariate ---", None),
        ("Correlation Matrix", "correlation"),
        ("PCA (2D)", "pca"),
        ("t-SNE (2D)", "tsne"),
        ("UMAP (2D)", "umap"),
        ("--- Other ---", None),
        ("Summary Table Only", "table"),
    ]

    def __init__(self, parent, store: DataStore):
        super().__init__(parent)
        self.title("Exploratory Data Analysis")
        self.geometry("1100x800")
        self.store = store

        self._filtered_rows = []
        self._filtered_labels = []
        self._current_idx = 0
        self._is_multivariate = False
        self._mv_settings_applied = False      
        self._mv_selected_vars = []            
        self._mv_impute_method = "drop"         

        if not HAS_MATPLOTLIB:
            messagebox.showwarning("Missing dependency",
                "matplotlib is required for EDA charts.\n\nInstall it with:\n  pip install matplotlib", parent=self)

        self._build_ui()
        self._refresh_patients()
        
        # Bind mousewheel to the left canvas and all its children after UI is built
        self._bind_left_panel_scroll()
        
        self.protocol("WM_DELETE_WINDOW", self.destroy)

    # ------------------------------------------------------------------ data
    def _collect_records(self):
        rows = []
        for pid, p in self.store.data.get("patients", {}).items():
            name = p.get("name", "")
            for v_idx, visit in enumerate(p.get("visits", [])):
                vdate = visit.get("visit_date") or ""
                for f_idx, f in enumerate(visit.get("fields", [])):
                    try:
                        val = float(str(f.get("value", "")).replace(",", "."))
                        rows.append({
                            "patient_id": pid, "patient_name": name, "visit_date": vdate,
                            "label": f.get("label", ""), "value": val,
                            "unit": f.get("unit", ""), "confidence": f.get("confidence", ""),
                            "_v_idx": v_idx, "_f_idx": f_idx
                        })
                    except (ValueError, TypeError): continue
        return rows

    def _unique_labels(self, rows): return sorted({r["label"] for r in rows if r["label"]})
    def _unique_patients(self, rows): return sorted({r["patient_id"] for r in rows})

    def _filter_rows(self, rows):
        sel_pids = self._selected_patient_ids()
        sel_labels = self._selected_labels()
        filtered = rows
        if sel_pids: filtered = [r for r in filtered if r["patient_id"] in sel_pids]
        if sel_labels: filtered = [r for r in filtered if r["label"] in sel_labels]
        return filtered

    def _build_feature_matrix(self, sel_labels, impute_method="drop"):
        if len(sel_labels) < 2:
            return None, "Select at least 2 variables for multivariate analysis."
        
        sel_pids = self._selected_patient_ids()
        raw_matrix, pids = [], []
        
        for pid, p in self.store.data.get("patients", {}).items():
            if sel_pids and pid not in sel_pids: continue
            for v_idx, visit in enumerate(p.get("visits", [])):
                row = {lbl: np.nan for lbl in sel_labels}
                has_any = False
                for f in visit.get("fields", []):
                    lbl = f.get("label", "")
                    if lbl in row:
                        try: val = float(str(f.get("value", "")).replace(",", "."))
                        except: val = np.nan
                        if not np.isnan(val): has_any = True
                        row[lbl] = val
                if has_any:
                    raw_matrix.append(row)
                    pids.append(pid)
                    
        if len(raw_matrix) < 2:
            return None, "Need at least 2 visits with data."
            
        X = np.array([[row[lbl] for lbl in sel_labels] for row in raw_matrix])
        handled_count = 0

        if impute_method == "drop":
            # Original behavior: completely remove visits with any missing data
            valid_mask = ~np.isnan(X).any(axis=1)
            X, valid_pids = X[valid_mask], [p for p, v in zip(pids, valid_mask) if v]
            handled_count = len(raw_matrix) - len(X)
        else:
            # Apply Imputation methods
            valid_pids = pids
            raw_X_copy = X.copy()
            
            if impute_method == "median":
                imputer = SimpleImputer(strategy='median')
            elif impute_method == "mean":
                imputer = SimpleImputer(strategy='mean')
            elif impute_method == "knn":
                imputer = KNNImputer(n_neighbors=5) # Default to 5 nearest neighbors
            elif impute_method == "mice":
                imputer = IterativeImputer(max_iter=10, random_state=42)
            else:
                imputer = SimpleImputer(strategy='median')
                
            try:
                X = imputer.fit_transform(X)
                # Count how many individual values were imputed (were NaN)
                handled_count = int(np.isnan(raw_X_copy).sum())
            except Exception as e:
                return None, f"Imputation failed: {str(e)}"

        if len(X) < 2:
            return None, f"Not enough data left after processing."
            
        return X, sel_labels, valid_pids, handled_count

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        # === LEFT OUTER FRAME (Fixed width container) ===
        left_outer = ttk.Frame(self, width=320)
        left_outer.pack(side=tk.LEFT, fill=tk.Y, padx=6, pady=6)
        left_outer.pack_propagate(False)

        # --- PINNED BOTTOM BUTTONS ---
        bottom_frame = ttk.Frame(left_outer)
        bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=(10, 0))
        
        # Use a custom style for the main action button to make it stand out
        style = ttk.Style()
        style.configure("Run.TButton", font=("", 10, "bold"))
        
        ttk.Button(bottom_frame, text="▶  Run Analysis", command=self.run_analysis, style="Run.TButton").pack(fill=tk.X, pady=(0, 4))
        ttk.Button(bottom_frame, text="Export filtered CSV…", command=self._export_filtered_csv).pack(fill=tk.X)

        # --- SCROLLABLE AREA ---
        canvas_container = ttk.Frame(left_outer)
        canvas_container.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self._left_canvas = tk.Canvas(canvas_container, highlightthickness=0, borderwidth=0)
        left_vsb = ttk.Scrollbar(canvas_container, orient="vertical", command=self._left_canvas.yview)
        self._left_canvas.configure(yscrollcommand=left_vsb.set)
        
        self._left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        left_vsb.pack(side=tk.RIGHT, fill=tk.Y)

        # Inner frame that actually holds the widgets
        left = ttk.Frame(self._left_canvas)
        self._left_window = self._left_canvas.create_window((0, 0), window=left, anchor="nw")

        def _on_left_configure(e):
            self._left_canvas.configure(scrollregion=self._left_canvas.bbox("all"))
        left.bind("<Configure>", _on_left_configure)

        def _on_canvas_configure(e):
            # Force the inner frame width to match the canvas width
            self._left_canvas.itemconfig(self._left_window, width=e.width)
        self._left_canvas.bind("<Configure>", _on_canvas_configure)

        # --- WIDGETS INSIDE SCROLLABLE FRAME ---

        # Patient selection
        pf = ttk.LabelFrame(left, text="Patients")
        pf.pack(fill=tk.X, pady=(0, 6), padx=2)
        btn_row = ttk.Frame(pf); btn_row.pack(fill=tk.X, padx=4, pady=2)
        ttk.Button(btn_row, text="Select All", command=self._select_all_patients).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_row, text="Deselect All", command=self._deselect_all_patients).pack(side=tk.LEFT, padx=2)
        self.patient_listbox = tk.Listbox(pf, selectmode=tk.EXTENDED, height=5, exportselection=False)
        self.patient_listbox.pack(fill=tk.X, padx=4, pady=4)

        # Field selection
        lf = ttk.LabelFrame(left, text="Fields (labels)")
        lf.pack(fill=tk.X, pady=(0, 6), padx=2)
        btn_row2 = ttk.Frame(lf); btn_row2.pack(fill=tk.X, padx=4, pady=2)
        ttk.Button(btn_row2, text="Select All", command=self._select_all_labels).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_row2, text="Deselect All", command=self._deselect_all_labels).pack(side=tk.LEFT, padx=2)
        self.label_listbox = tk.Listbox(lf, selectmode=tk.EXTENDED, height=6, exportselection=False)
        self.label_listbox.pack(fill=tk.X, padx=4, pady=4)

        # Plot type
        ptf = ttk.LabelFrame(left, text="Plot type")
        ptf.pack(fill=tk.X, pady=(0, 6), padx=2)
        self.plot_type_var = tk.StringVar(value=self.PLOT_TYPES[0][1])
        for display, val in self.PLOT_TYPES:
            if val is None: 
                ttk.Separator(ptf, orient="horizontal").pack(fill=tk.X, pady=6, padx=10)
                continue
            ttk.Radiobutton(ptf, text=display, variable=self.plot_type_var, value=val,
                            command=self.run_analysis).pack(anchor="w", padx=8, pady=1)

        # Options
        of = ttk.LabelFrame(left, text="Options (Univariate)")
        of.pack(fill=tk.X, pady=(0, 6), padx=2)
        self.show_grid_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(of, text="Show grid", variable=self.show_grid_var, command=self._render_current_plot).pack(anchor="w", padx=8)
        self.show_markers_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(of, text="Show markers", variable=self.show_markers_var, command=self._render_current_plot).pack(anchor="w", padx=8)
        self.combine_patients_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(of, text="Combine patients", variable=self.combine_patients_var, command=self._render_current_plot).pack(anchor="w", padx=8)

        # === RIGHT PANEL ===
        right = ttk.Frame(self)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=6, pady=6)
        self.notebook = ttk.Notebook(right)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        plot_tab = ttk.Frame(self.notebook)
        self.notebook.add(plot_tab, text="Plot")

        # Carousel Navigation
        nav_frame = ttk.Frame(plot_tab)
        nav_frame.pack(fill=tk.X, pady=(0, 5))
        self.prev_btn = ttk.Button(nav_frame, text="◀ Prev", command=self._prev_plot, width=10)
        self.prev_btn.pack(side=tk.LEFT, padx=4)
        self.next_btn = ttk.Button(nav_frame, text="Next ▶", command=self._next_plot, width=10)
        self.next_btn.pack(side=tk.LEFT, padx=4)
        
        self.nav_combo_var = tk.StringVar()
        self.nav_combo = ttk.Combobox(nav_frame, textvariable=self.nav_combo_var, state="readonly", width=50)
        self.nav_combo.pack(side=tk.LEFT, padx=10, fill=tk.X, expand=True)
        self.nav_combo.bind("<<ComboboxSelected>>", lambda e: self._on_nav_combo_change())
        
        self.nav_label_var = tk.StringVar(value="")
        ttk.Label(nav_frame, textvariable=self.nav_label_var, width=16, anchor="center", font=("", 9, "bold")).pack(side=tk.RIGHT, padx=4)

        # Matplotlib Canvas
        self.fig = Figure(figsize=(8, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=plot_tab)
        self.canvas_plot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas_plot.get_tk_widget().bind("<Left>", lambda e: self._prev_plot())
        self.canvas_plot.get_tk_widget().bind("<Right>", lambda e: self._next_plot())

        toolbar_frame = ttk.Frame(plot_tab)
        toolbar_frame.pack(fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas_plot, toolbar_frame)
        self.toolbar.update()

        # Summary Tab
        summary_tab = ttk.Frame(self.notebook)
        self.notebook.add(summary_tab, text="Summary Statistics")
        cols = ("label", "count", "mean", "std", "min", "q25", "median", "q75", "max", "unit")
        self.summary_tree = ttk.Treeview(summary_tab, columns=cols, show="headings", height=20)
        for c in cols:
            w = 60 if c != "label" else 180
            self.summary_tree.heading(c, text=c.upper())
            self.summary_tree.column(c, width=w, anchor="center" if c != "label" else "w")
        sum_vsb = ttk.Scrollbar(summary_tab, orient="vertical", command=self.summary_tree.yview)
        self.summary_tree.configure(yscrollcommand=sum_vsb.set)
        self.summary_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sum_vsb.pack(side=tk.RIGHT, fill=tk.Y)

        self.info_label = ttk.Label(right, text="", foreground="#555555")
        self.info_label.pack(fill=tk.X, pady=(4, 0))

        # --- Data Coverage Tab ---
        coverage_tab = ttk.Frame(self.notebook)
        self.notebook.add(coverage_tab, text="Data Coverage")

        cov_label = ttk.Label(coverage_tab, 
            text="Variable Presence Summary (Applies current Patient/Label/Date filters):", 
            font=("", 9, "bold"))
        cov_label.pack(anchor="w", padx=5, pady=(10, 2))

        cov_tree_frame = ttk.Frame(coverage_tab)
        cov_tree_frame.pack(fill=tk.BOTH, expand=True, padx=5)

        cov_cols = ("variable", "pat_total", "pat_has", "pat_miss", "pat_pct", 
                    "vis_total", "vis_has", "vis_miss", "vis_pct")
        self.cov_tree = ttk.Treeview(cov_tree_frame, columns=cov_cols, show="headings", height=8)
        
        cov_headers = {"variable": ("Variable", 160), "pat_total": ("Total Pts", 65), 
                       "pat_has": ("Pts With", 65), "pat_miss": ("Pts Miss", 65), 
                       "pat_pct": ("Pt Cov %", 70), "vis_total": ("Total Vis", 65),
                       "vis_has": ("Vis With", 65), "vis_miss": ("Vis Miss", 65), 
                       "vis_pct": ("Vis Cov %", 70)}
        
        for col, (text, width) in cov_headers.items():
            anchor = "w" if col == "variable" else "center"
            self.cov_tree.heading(col, text=text)
            self.cov_tree.column(col, width=width, anchor=anchor, stretch=False)

        cov_vsb = ttk.Scrollbar(cov_tree_frame, orient="vertical", command=self.cov_tree.yview)
        self.cov_tree.configure(yscrollcommand=cov_vsb.set)
        self.cov_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        cov_vsb.pack(side=tk.RIGHT, fill=tk.Y)

        # Matrix Section
        mat_label = ttk.Label(coverage_tab, 
            text="Patient-Variable Matrix (Variables = Rows, Patients = Columns. ✔ = Has, ✖ = Missing):", 
            font=("", 9, "bold"))
        mat_label.pack(anchor="w", padx=5, pady=(20, 2))

        mat_tree_frame = ttk.Frame(coverage_tab)
        mat_tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 10))
        mat_tree_frame.rowconfigure(0, weight=1)
        mat_tree_frame.columnconfigure(0, weight=1)

        self.matrix_tree = ttk.Treeview(mat_tree_frame, show="headings", height=10)
        mat_vsb = ttk.Scrollbar(mat_tree_frame, orient="vertical", command=self.matrix_tree.yview)
        mat_hsb = ttk.Scrollbar(mat_tree_frame, orient="horizontal", command=self.matrix_tree.xview)
        self.matrix_tree.configure(yscrollcommand=mat_vsb.set, xscrollcommand=mat_hsb.set)
        
        self.matrix_tree.grid(row=0, column=0, sticky="nsew")
        mat_vsb.grid(row=0, column=1, sticky="ns")
        mat_hsb.grid(row=1, column=0, sticky="ew")

    def _bind_left_panel_scroll(self):
        """Recursively bind mousewheel events to the left canvas and all its children."""
        def _on_mousewheel(event):
            if event.num == 4: delta = -1
            elif event.num == 5: delta = 1
            else: delta = int(event.delta / 120) * -1
            self._left_canvas.yview_scroll(delta, "units")

        def _bind_to_widget(widget):
            widget.bind("<MouseWheel>", _on_mousewheel, add="+")
            widget.bind("<Button-4>", _on_mousewheel, add="+")
            widget.bind("<Button-5>", _on_mousewheel, add="+")
            for child in widget.winfo_children():
                _bind_to_widget(child)
                
        _bind_to_widget(self._left_canvas)

    # ------------------------------------------------------------------ helpers
    def _refresh_patients(self):
        self.patient_listbox.delete(0, tk.END)
        rows = self._collect_records()
        for pid in self._unique_patients(rows):
            name = next((r["patient_name"] for r in rows if r["patient_id"] == pid), "")
            self.patient_listbox.insert(tk.END, f"{pid}  –  {name}")
        self.patient_listbox.select_set(0, tk.END)
        self._refresh_labels(rows)

    def _refresh_labels(self, rows=None):
        if rows is None: rows = self._collect_records()
        self.label_listbox.delete(0, tk.END)
        for lbl in self._unique_labels(rows): self.label_listbox.insert(tk.END, lbl)
        self.label_listbox.select_set(0, tk.END)

    def _selected_patient_ids(self):
        pids = self._unique_patients(self._collect_records())
        sel = self.patient_listbox.curselection()
        return [] if (not sel or len(sel) == len(pids)) else [pids[i] for i in sel if i < len(pids)]

    def _selected_labels(self):
        labels = [self.label_listbox.get(i) for i in range(self.label_listbox.size())]
        sel = self.label_listbox.curselection()
        return [] if (not sel or len(sel) == len(labels)) else [labels[i] for i in sel if i < len(labels)]

    def _select_all_patients(self): self.patient_listbox.select_set(0, tk.END)
    def _deselect_all_patients(self): self.patient_listbox.select_clear(0, tk.END)
    def _select_all_labels(self): self.label_listbox.select_set(0, tk.END)
    def _deselect_all_labels(self): self.label_listbox.select_clear(0, tk.END)

    def _compute_summary(self, rows):
        from collections import defaultdict
        groups = defaultdict(list)
        for r in rows: groups[r["label"]].append(r["value"])
        stats = []
        for label in sorted(groups):
            vals = sorted(groups[label]); n = len(vals)
            if n == 0: continue
            arr = np.array(vals)
            unit = next((r["unit"] for r in rows if r["label"] == label and r.get("unit")), "")
            stats.append({"label": label, "count": n, "mean": float(np.mean(arr)),
                "std": float(np.std(arr, ddof=1)) if n > 1 else 0.0,
                "min": float(np.min(arr)), "q25": float(np.percentile(arr, 25)),
                "median": float(np.median(arr)), "q75": float(np.percentile(arr, 75)),
                "max": float(np.max(arr)), "unit": unit})
        return stats

    def _populate_summary_tree(self, stats):
        self.summary_tree.delete(*self.summary_tree.get_children())
        for s in stats:
            self.summary_tree.insert("", tk.END, values=(
                s["label"], s["count"], f"{s['mean']:.3f}", f"{s['std']:.3f}",
                f"{s['min']:.3f}", f"{s['q25']:.3f}", f"{s['median']:.3f}",
                f"{s['q75']:.3f}", f"{s['max']:.3f}", s["unit"]))

    def update_coverage_tab(self):
        """Calculate missing values and populate the coverage summary & matrix."""
        from collections import defaultdict
        
        sel_pids = self._selected_patient_ids()
        sel_labels = self._selected_labels()

        active_pids = set()
        active_labels = set()
        total_visits = 0
        
        var_pat_set = defaultdict(set)  # {label: set(pids)}
        var_vis_set = defaultdict(set)  # {label: set(visit_ids)}

        for pid, p in self.store.data.get("patients", {}).items():
            if sel_pids and pid not in sel_pids: continue
            active_pids.add(pid)
            
            for v_idx, visit in enumerate(p.get("visits", [])):              
                total_visits += 1
                vis_id = f"{pid}_v{v_idx}"
                
                for f in visit.get("fields", []):
                    lbl = f.get("label", "")
                    if not lbl: continue
                    if sel_labels and lbl not in sel_labels: continue
                    
                    active_labels.add(lbl)
                    var_pat_set[lbl].add(pid)
                    var_vis_set[lbl].add(vis_id)

        # 1. Populate Summary Treeview
        self.cov_tree.delete(*self.cov_tree.get_children())
        for lbl in sorted(active_labels):
            p_has = len(var_pat_set.get(lbl, set()))
            p_miss = len(active_pids) - p_has
            p_pct = (p_has / len(active_pids) * 100) if active_pids else 0
            
            v_has = len(var_vis_set.get(lbl, set()))
            v_miss = total_visits - v_has
            v_pct = (v_has / total_visits * 100) if total_visits else 0
            
            self.cov_tree.insert("", tk.END, values=(
                lbl, len(active_pids), p_has, p_miss, f"{p_pct:.1f}%",
                total_visits, v_has, v_miss, f"{v_pct:.1f}%"
            ))

        # 2. Populate Dynamic Patient-Variable Matrix
        self.matrix_tree.delete(*self.matrix_tree.get_children())
        
        pids_list = sorted(active_pids)
        labels_list = sorted(active_labels)
        
        # Dynamically set columns (Variable + 1 column per patient)
        cols = ("variable",) + tuple(pids_list)
        self.matrix_tree["columns"] = cols
        
        self.matrix_tree.heading("variable", text="Variable")
        self.matrix_tree.column("variable", width=150, anchor="w", stretch=False)
        
        for pid in pids_list:
            self.matrix_tree.heading(pid, text=pid)
            self.matrix_tree.column(pid, width=60, anchor="center", stretch=False)
            
        # Insert rows
        for lbl in labels_list:
            row_vals = [lbl]
            for pid in pids_list:
                # Use checkmarks for better visual distinction
                row_vals.append("✔" if pid in var_pat_set.get(lbl, set()) else "✖")
            self.matrix_tree.insert("", tk.END, values=row_vals)

    # ------------------------------------------------------------------ carousel controls
    def _prev_plot(self):
        if self._is_multivariate or not self._filtered_labels: return
        self._current_idx = max(0, self._current_idx - 1)
        self.nav_combo.current(self._current_idx)
        self._render_current_plot()

    def _next_plot(self):
        if self._is_multivariate or not self._filtered_labels: return
        self._current_idx = min(len(self._filtered_labels) - 1, self._current_idx + 1)
        self.nav_combo.current(self._current_idx)
        self._render_current_plot()

    def _on_nav_combo_change(self):
        if self._is_multivariate: return
        idx = self.nav_combo.current()
        if idx >= 0: self._current_idx = idx; self._render_current_plot()

    # ------------------------------------------------------------------ plotting logic
    def run_analysis(self):
        rows = self._collect_records()
        self._filtered_rows = self._filter_rows(rows)
        plot_type = self.plot_type_var.get()

        self.info_label.config(text=f"Showing {len(self._filtered_rows)} record(s) out of {len(rows)} total.")
        stats = self._compute_summary(self._filtered_rows)
        self._populate_summary_tree(stats)

        if plot_type == "table": self.notebook.select(1); return
        if not HAS_MATPLOTLIB: self.notebook.select(1); return

        multivariate_types = {"correlation", "pca", "tsne", "umap"}
        self._is_multivariate = plot_type in multivariate_types

        if self._is_multivariate:
            self._mv_settings_applied = False 
            self.nav_combo["values"] = ["Global View"]
            self.nav_combo.current(0)
            self.nav_label_var.set("Global View")
            self.prev_btn.config(state="disabled")
            self.next_btn.config(state="disabled")
            self._render_current_plot()
        else:
            self.prev_btn.config(state="normal")
            self.next_btn.config(state="normal")
            self._filtered_labels = sorted({r["label"] for r in self._filtered_rows}) if self._filtered_rows else []
            
            if not self._filtered_labels:
                self.fig.clear(); self.ax = self.fig.add_subplot(111)
                self.ax.set_title("No data matching filters"); self.canvas_plot.draw()
                self.notebook.select(0); return

            self.nav_combo["values"] = self._filtered_labels
            if self._current_idx >= len(self._filtered_labels): self._current_idx = 0
            self.nav_combo.current(self._current_idx)
            self._render_current_plot()
            
        self.notebook.select(0)
        self.update_coverage_tab()

    def _render_current_plot(self):
        if not HAS_MATPLOTLIB: return
        if self._is_multivariate:
            self._render_multivariate()
        else:
            if not self._filtered_labels: return
            label = self._filtered_labels[self._current_idx]
            plot_type = self.plot_type_var.get()
            self.fig.clear(); self.ax = self.fig.add_subplot(111)
            label_rows = [r for r in self._filtered_rows if r["label"] == label]
            
            if plot_type == "timeseries": self._plot_timeseries(label_rows, label)
            elif plot_type == "histogram": self._plot_histogram(label_rows, label)
            elif plot_type == "boxplot": self._plot_boxplot(label_rows, label)
            elif plot_type == "bar": self._plot_bar(label_rows, label)
            
            self.fig.tight_layout(); self.canvas_plot.draw()
            self.nav_label_var.set(f"{self._current_idx + 1} / {len(self._filtered_labels)}")

    # ------------------------------------------------------------------ MULTIVARIATE PLOTS
    def _render_multivariate(self):
        # If settings aren't applied, pause rendering and open the setup dialog
        if not self._mv_settings_applied:
            self._open_mv_setup()
            return

        # Actual calculation and plotting
        plot_type = self.plot_type_var.get()
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.info_label.config(text="Processing multivariate analysis...")
        self.update_idletasks()

        if plot_type == "distribution":
            self._plot_distribution()
        else:
            # Pass the user-selected variables and imputation method to the matrix builder
            res = self._build_feature_matrix(self._mv_selected_vars, self._mv_impute_method)
            if res[0] is None:
                self.ax.set_title(res[1], fontsize=12, color='red')
                self.canvas_plot.draw()
                self.info_label.config(text=res[1])
                return
            
            X, feature_names, pids, dropped = res
            action_word = "Dropped" if self._mv_impute_method == "drop" else "Imputed"
            if dropped > 0:
                self.info_label.config(text=f"{action_word} {dropped} missing value(s) using selected method.")
            
            if plot_type == "correlation":
                self._plot_correlation(X, feature_names)
            elif plot_type == "pca":
                if not HAS_SKLEARN:
                    self.ax.set_title("scikit-learn is required for PCA", color='red')
                    self.canvas_plot.draw()
                    return
                self._plot_pca(X, feature_names, pids)
            elif plot_type == "tsne":
                if not HAS_SKLEARN:
                    self.ax.set_title("scikit-learn is required for t-SNE", color='red')
                    self.canvas_plot.draw()
                    return
                self._plot_tsne(X, feature_names, pids)
            elif plot_type == "umap":
                if not HAS_UMAP:
                    self.ax.set_title("umap-learn is required for UMAP (pip install umap-learn)", color='red')
                    self.canvas_plot.draw()
                    return
                self._plot_umap(X, feature_names, pids)

        self.fig.tight_layout()
        self.canvas_plot.draw()

    def _open_mv_setup(self):
        """Open the setup dialog and process the result."""
        from collections import Counter
        counts = Counter(r["label"] for r in self._filtered_rows if r["label"])
        vars_with_counts = sorted(counts.items(), key=lambda x: (-x[1], x[0]))

        # Pre-calculate which variables exist in which visit for the smart requirement checker
        visit_coverage = {}
        for r in self._filtered_rows:
            key = (r["patient_id"], r["_v_idx"])
            if key not in visit_coverage:
                visit_coverage[key] = set()
            visit_coverage[key].add(r["label"])

        dialog = MultivariateSetupDialog(self, vars_with_counts, visit_coverage)
        self.wait_window(dialog)  # Block until user clicks Run or Cancel

        if dialog.result:
            self._mv_selected_vars, method_str = dialog.result
            
            # Map dialog string to internal code
            method_map = {
                "Drop missing visits (Complete Cases)": "drop",
                "Impute with Median": "median",
                "Impute with Mean": "mean",
                "Impute with Zero": "zero"
            }
            self._mv_impute_method = method_map.get(method_str, "drop")
            self._mv_settings_applied = True
            
            # Render again now that settings are applied
            self._render_multivariate()
        else:
            # User cancelled, revert to a safe univariate state
            self.plot_type_var.set("timeseries")
            self._is_multivariate = False
            self.run_analysis()

    def _plot_correlation(self, X, feature_names):
        corr = np.corrcoef(X.T)
        im = self.ax.imshow(corr, cmap='coolwarm', vmin=-1, vmax=1)
        self.fig.colorbar(im, ax=self.ax)
        self.ax.set_xticks(range(len(feature_names)))
        self.ax.set_yticks(range(len(feature_names)))
        self.ax.set_xticklabels(feature_names, rotation=45, ha='right', fontsize=8)
        self.ax.set_yticklabels(feature_names, fontsize=8)
        for i in range(len(feature_names)):
            for j in range(len(feature_names)):
                text_color = "white" if abs(corr[i, j]) > 0.7 else "black"
                self.ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", color=text_color, fontsize=8)
        self.ax.set_title("Feature Correlation Matrix", fontsize=12, fontweight='bold')

    def _plot_pca(self, X, feature_names, pids):
        if X.shape[0] < 2:
            self.ax.set_title("PCA requires at least 2 visits", fontsize=12, color='red')
            return
            
        X_scaled = StandardScaler().fit_transform(X)
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X_scaled)  # This is the line that was missing!
        self._scatter_multivariate(X_pca, pids, "PCA (2D)", f"Explained Var: {sum(pca.explained_variance_ratio_)*100:.1f}%")

    def _plot_tsne(self, X, feature_names, pids):
        n_samples = X.shape[0]
        
        # t-SNE requires perplexity < n_samples. We dynamically lower it if needed.
        default_perplexity = 30.0
        max_possible_perplexity = n_samples - 1
        
        if max_possible_perplexity < 5:
            # Not enough data to run t-SNE meaningfully
            self.ax.set_title(f"t-SNE requires at least 5 visits (you have {n_samples})", fontsize=12, color='red')
            return

        safe_perplexity = min(default_perplexity, max_possible_perplexity)
        
        X_scaled = StandardScaler().fit_transform(X)
        tsne = TSNE(n_components=2, perplexity=safe_perplexity, random_state=42, learning_rate='auto', init='pca')
        X_tsne = tsne.fit_transform(X_scaled)
        self._scatter_multivariate(X_tsne, pids, "t-SNE (2D)", f"Perplexity adjusted to {safe_perplexity:.1f}")

    def _plot_umap(self, X, feature_names, pids):
        # UMAP requires a reasonable number of samples to build its neighbor graph.
        # If we have too few, it crashes on internal sparse matrix math.
        if X.shape[0] < 5:
            self.ax.set_title(
                f"UMAP requires at least 5 visits to build a reliable graph (you have {X.shape[0]})", 
                fontsize=12, color='red'
            )
            return

        X_scaled = StandardScaler().fit_transform(X)
        
        # Dynamically lower n_neighbors if the dataset is small
        # UMAP default is 15, but it must be < n_samples
        safe_n_neighbors = min(15, X.shape[0] - 1)
        
        reducer = UMAP(n_components=2, n_neighbors=safe_n_neighbors, random_state=42)
        X_umap = reducer.fit_transform(X_scaled)
        self._scatter_multivariate(X_umap, pids, "UMAP (2D)", f"n_neighbors adjusted to {safe_n_neighbors}")

    def _scatter_multivariate(self, X_transformed, pids, title, subtitle=""):
        unique_pids = list(sorted(set(pids)))
        pid_to_color = {pid: i for i, pid in enumerate(unique_pids)}
        colors = [plt.cm.tab10(pid_to_color[p] % 10) for p in pids]
        
        self.ax.scatter(X_transformed[:, 0], X_transformed[:, 1], c=colors, alpha=0.7, edgecolors='w', linewidth=0.5)
        
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=plt.cm.tab10(i%10), markersize=10, label=pid) for i, pid in enumerate(unique_pids)]
        self.ax.legend(handles=handles, fontsize=8, loc='best', title="Patient ID", bbox_to_anchor=(1.05, 1))
        
        full_title = f"{title}"
        if subtitle: full_title += f"\n{subtitle}"
        self.ax.set_title(full_title, fontsize=12, fontweight='bold')
        self.ax.set_xlabel("Component 1"); self.ax.set_ylabel("Component 2")
        if self.show_grid_var.get(): self.ax.grid(True, alpha=0.3)

    # ------------------------------------------------------------------ UNIVARIATE PLOTS
    def _parse_date(self, d):
        if not d: return None
        for fmt in ("%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y"):
            try: return datetime.strptime(d, fmt)
            except ValueError: pass
        return None

    def _plot_timeseries(self, label_rows, label):
        combine = self.combine_patients_var.get()
        sel_pids = self._selected_patient_ids()
        pids = sel_pids if sel_pids else sorted({r["patient_id"] for r in label_rows})
        colors = plt.cm.tab10.colors
        if combine:
            pts = sorted([(self._parse_date(r["visit_date"]), r["value"]) for r in label_rows if self._parse_date(r["visit_date"])])
            if pts:
                xs, ys = zip(*pts)
                self.ax.plot(xs, ys, marker="o" if self.show_markers_var.get() else None, linewidth=1.5, color=colors[0])
        else:
            for pi, pid in enumerate(pids):
                pts = sorted([(self._parse_date(r["visit_date"]), r["value"]) for r in label_rows if r["patient_id"] == pid and self._parse_date(r["visit_date"])])
                if pts:
                    xs, ys = zip(*pts)
                    mk = "o" if self.show_markers_var.get() else None
                    self.ax.plot(xs, ys, marker=mk, linewidth=1.5, label=pid, color=colors[pi % len(colors)])
            if len(pids) <= 12: self.ax.legend(fontsize=8, loc="best")
        unit = next((r["unit"] for r in label_rows if r.get("unit")), "")
        self.ax.set_title(label, fontsize=12, fontweight="bold")
        self.ax.set_ylabel(unit if unit else "value", fontsize=9)
        if self.show_grid_var.get(): self.ax.grid(True, alpha=0.3)
        self.ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        self.ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        for lbl_obj in self.ax.get_xticklabels(): lbl_obj.set_rotation(30); lbl_obj.set_ha("right")

    def _plot_histogram(self, label_rows, label):
        vals = [r["value"] for r in label_rows]
        if not vals: self.ax.set_title(f"{label} (no data)"); return
        unit = next((r["unit"] for r in label_rows if r.get("unit")), "")
        self.ax.hist(vals, bins=max(5, min(30, int(np.sqrt(len(vals))))), color=plt.cm.Set3.colors[0], edgecolor="white", alpha=0.85)
        self.ax.axvline(np.mean(vals), color="red", linestyle="--", linewidth=1.5, label=f"mean={np.mean(vals):.2f}")
        self.ax.axvline(np.median(vals), color="blue", linestyle=":", linewidth=1.5, label=f"median={np.median(vals):.2f}")
        self.ax.legend(fontsize=9); self.ax.set_title(label, fontsize=12, fontweight="bold")
        self.ax.set_xlabel(unit if unit else "value", fontsize=9)
        if self.show_grid_var.get(): self.ax.grid(True, alpha=0.3)

    def _plot_boxplot(self, label_rows, label):
        sel_pids = self._selected_patient_ids()
        pids = sel_pids if sel_pids else sorted({r["patient_id"] for r in label_rows})
        data, tick_labels = [], []
        for pid in pids:
            vals = [r["value"] for r in label_rows if r["patient_id"] == pid]
            if vals: 
                data.append(vals)
                tick_labels.append(pid)
                
        if not data: 
            self.ax.set_title(f"{label} (no data)", fontsize=12, fontweight="bold")
            return
            
        # Try the new Matplotlib 3.9+ argument first, fallback to the old one
        try:
            bp = self.ax.boxplot(data, patch_artist=True, tick_labels=tick_labels)
        except TypeError:
            bp = self.ax.boxplot(data, patch_artist=True, labels=tick_labels)
            
        # FIX: Use plt.get_cmap() which works safely in all modern Matplotlib versions
        cmap = plt.get_cmap('Pastel1')
        colors = [cmap(i % cmap.N) for i in range(len(data))]
        
        if "boxes" in bp:
            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)
                
        unit = next((r["unit"] for r in label_rows if r.get("unit")), "")
        self.ax.set_title(label, fontsize=12, fontweight="bold")
        self.ax.set_ylabel(unit if unit else "value", fontsize=9)
        self.ax.tick_params(axis="x", labelsize=8)
        for lbl in self.ax.get_xticklabels():
            lbl.set_rotation(30)
            lbl.set_ha("right")
        if self.show_grid_var.get(): 
            self.ax.grid(True, axis="y", alpha=0.3)

    def _plot_bar(self, label_rows, label):
        sel_pids = self._selected_patient_ids()
        pids = sel_pids if sel_pids else sorted({r["patient_id"] for r in label_rows})
        colors = plt.cm.tab10.colors; means, pid_labels = [], []
        for pid in pids:
            vals = [r["value"] for r in label_rows if r["patient_id"] == pid]
            if vals: means.append(np.mean(vals)); pid_labels.append(pid)
        if not means: self.ax.set_title(f"{label} (no data)"); return
        x = range(len(pid_labels))
        self.ax.bar(x, means, color=[colors[i % 10] for i in range(len(means))], edgecolor="white", alpha=0.85)
        self.ax.set_xticks(list(x)); self.ax.set_xticklabels(pid_labels, fontsize=8, rotation=30, ha="right")
        unit = next((r["unit"] for r in label_rows if r.get("unit")), "")
        self.ax.set_title(label, fontsize=12, fontweight="bold"); self.ax.set_ylabel(unit if unit else "value", fontsize=9)
        if self.show_grid_var.get(): self.ax.grid(True, axis="y", alpha=0.3)

    # ------------------------------------------------------------------ export
    def _export_filtered_csv(self):
        rows = self._collect_records(); filtered = self._filter_rows(rows)
        if not filtered: messagebox.showinfo("Export", "No records match filters.", parent=self); return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")], parent=self)
        if not path: return
        import csv as csv_mod
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv_mod.writer(f)
            w.writerow(["patient_id", "patient_name", "visit_date", "label", "value", "unit", "confidence"])
            for r in filtered: w.writerow([r["patient_id"], r["patient_name"], r["visit_date"], r["label"], r["value"], r["unit"], r["confidence"]])
        messagebox.showinfo("Export", f"Exported {len(filtered)} records to:\n{path}", parent=self)

class FixVariableNamesWindow(tk.Toplevel):
    """Window to detect and fix duplicate/variant variable names across all records."""

    def __init__(self, parent, store: DataStore):
        super().__init__(parent)
        self.title("Fix / Consolidate Variable Names")
        self.geometry("980x720")
        self.minsize(700, 500)
        self.store = store
        self.duplicates = []        # [(raw, canonical, count, method), ...]
        self.unmapped_labels = []   # [(raw, count), ...]
        self._build_ui()
        self._scan()

    # --------------------------------------------------------------- resolve
    def _resolve(self, raw_label):
        """Return (canonical_name_or_None, method_string)."""
        norm = normalize_label(raw_label)
        if not norm:
            return None, "empty"
        if norm in _LABEL_LOOKUP_NORM:
            return _LABEL_LOOKUP_NORM[norm], "exact"
        if norm in _BUILTIN_ALIAS_LOOKUP_NORM:
            return _BUILTIN_ALIAS_LOOKUP_NORM[norm], "alias"
        if norm in self.store.aliases:
            return self.store.aliases[norm]["canonical"], "user_alias"
        best = difflib.get_close_matches(
            norm, _LABEL_LOOKUP_NORM.keys(), n=1, cutoff=0.70
        )
        if best:
            return _LABEL_LOOKUP_NORM[best[0]], "fuzzy"
        return None, "unmapped"

    # --------------------------------------------------------------- scan
    def _scan(self):
        """Walk every saved field, count raw labels, and classify them."""
        label_counts = {}  # raw_label -> occurrence count
        for pid, p in self.store.data.get("patients", {}).items():
            for visit in p.get("visits", []):
                for f in visit.get("fields", []):
                    lbl = f.get("label", "").strip()
                    if lbl:
                        label_counts[lbl] = label_counts.get(lbl, 0) + 1

        self.duplicates = []
        self.unmapped_labels = []

        for raw, count in sorted(label_counts.items()):
            canonical, method = self._resolve(raw)
            if canonical is None:
                self.unmapped_labels.append((raw, count))
            elif method != "exact":
                self.duplicates.append((raw, canonical, count, method))

        self._populate()

    # --------------------------------------------------------------- UI
    def _build_ui(self):
        main = ttk.Frame(self, padding=10)
        main.pack(fill=tk.BOTH, expand=True)

        # ---- explanatory text ----
        ttk.Label(
            main,
            text=(
                "This tool finds labels that refer to the same measurement but are "
                "stored under different names (e.g. \"CCMH\" and \"MCHC\" both mean "
                "C.C.M.H.).\n"
                "The correct canonical names come from the built-in label list and "
                "your synonym dictionary."
            ),
            wraplength=900,
            foreground="#444444",
        ).pack(anchor="w", pady=(0, 8))

        # ============ TOP: auto-detectable duplicates ============
        dup_lf = ttk.LabelFrame(main, text="Auto-detectable duplicates")
        dup_lf.pack(fill=tk.BOTH, expand=True, pady=(0, 8))

        self.dup_count_label = ttk.Label(dup_lf, text="", font=("", 9, "bold"))
        self.dup_count_label.pack(anchor="w", padx=6, pady=(4, 2))

        tree_frame = ttk.Frame(dup_lf)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=6, pady=4)

        cols = ("variant", "arrow", "canonical", "count", "method")
        self.dup_tree = ttk.Treeview(tree_frame, columns=cols, show="headings",
                                      selectmode="extended", height=8)
        self.dup_tree.heading("variant",   text="Variant found in data")
        self.dup_tree.heading("arrow",     text="")
        self.dup_tree.heading("canonical", text="Should be renamed to →")
        self.dup_tree.heading("count",     text="Records")
        self.dup_tree.heading("method",    text="Match type")
        self.dup_tree.column("variant",   width=240, anchor="w")
        self.dup_tree.column("arrow",     width=30,  anchor="center", stretch=False)
        self.dup_tree.column("canonical", width=260, anchor="w")
        self.dup_tree.column("count",     width=70,  anchor="center", stretch=False)
        self.dup_tree.column("method",    width=90,  anchor="center", stretch=False)

        dup_vsb = ttk.Scrollbar(tree_frame, orient="vertical",
                                 command=self.dup_tree.yview)
        self.dup_tree.configure(yscrollcommand=dup_vsb.set)
        self.dup_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        dup_vsb.pack(side=tk.LEFT, fill=tk.Y)

        # Tag colours for match method
        self.dup_tree.tag_configure("alias",      foreground="#2980b9")
        self.dup_tree.tag_configure("user_alias", foreground="#8e44ad")
        self.dup_tree.tag_configure("fuzzy",      foreground="#d35400")

        btn_col = ttk.Frame(dup_lf)
        btn_col.pack(side=tk.RIGHT, fill=tk.Y, padx=(0, 6), pady=4)
        ttk.Button(btn_col, text="Select All",
                    command=self._sel_all_dup, width=16).pack(pady=2)
        ttk.Button(btn_col, text="Deselect All",
                    command=self._desel_all_dup, width=16).pack(pady=2)
        ttk.Separator(btn_col, orient="horizontal").pack(fill=tk.X, pady=8)
        ttk.Button(btn_col, text="Fix Selected",
                    command=self._fix_selected_dup, width=16).pack(pady=2)
        ttk.Button(btn_col, text="Fix All",
                    command=self._fix_all_dup, width=16).pack(pady=2)

        # ============ BOTTOM: unmapped labels ============
        unm_lf = ttk.LabelFrame(
            main,
            text="Unmapped labels (no automatic match — assign manually)"
        )
        unm_lf.pack(fill=tk.BOTH, expand=True, pady=(0, 8))

        self.unm_count_label = ttk.Label(unm_lf, text="", font=("", 9, "bold"))
        self.unm_count_label.pack(anchor="w", padx=6, pady=(4, 2))

        unm_body = ttk.Frame(unm_lf)
        unm_body.pack(fill=tk.BOTH, expand=True, padx=6, pady=4)

        # Left: listbox of unmapped labels
        left_unm = ttk.Frame(unm_body)
        left_unm.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.unm_listbox = tk.Listbox(
            left_unm, height=6, width=45, exportselection=False,
            font=("TkFixedFont", 10)
        )
        self.unm_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.unm_listbox.bind("<<ListboxSelect>>", self._on_unm_select)
        unm_vsb = ttk.Scrollbar(left_unm, orient="vertical",
                                 command=self.unm_listbox.yview)
        self.unm_listbox.configure(yscrollcommand=unm_vsb.set)
        unm_vsb.pack(side=tk.LEFT, fill=tk.Y)

        # Right: assignment controls
        right_unm = ttk.Frame(unm_body)
        right_unm.pack(side=tk.LEFT, fill=tk.Y, padx=(16, 0))

        ttk.Label(right_unm, text="Assign selected label to:").pack(anchor="w")
        self.unm_combo = ttk.Combobox(
            right_unm, values=LABEL_SUGGESTIONS, width=38, font=("", 10)
        )
        self.unm_combo.pack(pady=4)

        # Fuzzy suggestion label
        self.fuzzy_hint = ttk.Label(
            right_unm, text="", foreground="#888888",
            wraplength=260, font=("", 9, "italic")
        )
        self.fuzzy_hint.pack(anchor="w", pady=(0, 4))

        ttk.Button(right_unm, text="Assign & Fix",
                    command=self._assign_unmapped, width=18).pack(pady=4)
        ttk.Button(right_unm, text="Assign & Learn Alias",
                    command=self._assign_and_learn, width=18).pack(pady=2)

        # ---- status bar ----
        sep = ttk.Separator(main, orient="horizontal")
        sep.pack(fill=tk.X, pady=6)

        status_row = ttk.Frame(main)
        status_row.pack(fill=tk.X)
        self.status_var = tk.StringVar(value="Ready.")
        self.status_label = ttk.Label(
            status_row, textvariable=self.status_var,
            foreground="#006600", font=("", 10)
        )
        self.status_label.pack(side=tk.LEFT, padx=4)
        ttk.Button(status_row, text="Close", command=self.destroy,
                    width=10).pack(side=tk.RIGHT)

    # --------------------------------------------------------------- populate
    def _populate(self):
        # -- duplicates tree --
        self.dup_tree.delete(*self.dup_tree.get_children())
        for raw, canonical, count, method in self.duplicates:
            tag = method if method in ("alias", "user_alias", "fuzzy") else ""
            self.dup_tree.insert(
                "", tk.END,
                values=(raw, "→", canonical, count, method),
                tags=(tag,)
            )
        self.dup_count_label.config(
            text=f"{len(self.duplicates)} variant(s) that map to a canonical name:"
        )

        # -- unmapped listbox --
        self.unm_listbox.delete(0, tk.END)
        for raw, count in self.unmapped_labels:
            self.unm_listbox.insert(tk.END, f"{raw}  ({count} record(s))")
        self.unm_count_label.config(
            text=f"{len(self.unmapped_labels)} label(s) with no automatic match:"
        )
        self.fuzzy_hint.config(text="")
        self.unm_combo.set("")

    # --------------------------------------------------------------- dup helpers
    def _sel_all_dup(self):
        self.dup_tree.selection_set(self.dup_tree.get_children())

    def _desel_all_dup(self):
        self.dup_tree.selection_remove(self.dup_tree.get_children())

    def _gather_dup_renames(self, selection_only=False):
        """Return list of (old_label, new_label) from the tree."""
        children = (self.dup_tree.selection() if selection_only
                    else self.dup_tree.get_children())
        renames = []
        for item in children:
            vals = self.dup_tree.item(item, "values")
            if vals and vals[0] and vals[2]:
                renames.append((str(vals[0]), str(vals[2])))
        return renames

    # --------------------------------------------------------------- apply
    def _apply_renames(self, renames):
        """Walk all records and rename matching labels. Return rename count."""
        if not renames:
            return 0
        rename_map = {old: new for old, new in renames}
        count = 0
        for pid, p in self.store.data.get("patients", {}).items():
            for visit in p.get("visits", []):
                for f in visit.get("fields", []):
                    lbl = f.get("label", "").strip()
                    if lbl in rename_map:
                        f["label"] = rename_map[lbl]
                        count += 1
        if count > 0:
            self.store.save()
        return count

    # --------------------------------------------------------------- dup actions
    def _fix_selected_dup(self):
        renames = self._gather_dup_renames(selection_only=True)
        if not renames:
            messagebox.showinfo("Fix", "No rows selected.", parent=self)
            return
        self._do_fix(renames)

    def _fix_all_dup(self):
        renames = self._gather_dup_renames(selection_only=False)
        if not renames:
            messagebox.showinfo("Fix", "No duplicates to fix.", parent=self)
            return
        self._do_fix(renames)

    def _do_fix(self, renames):
        summary = "\n".join(f"  {old}  →  {new}" for old, new in renames)
        confirm = messagebox.askyesno(
            "Confirm rename",
            f"About to rename the following variant(s) in all saved records:\n\n"
            f"{summary}\n\n"
            f"Proceed?",
            parent=self,
        )
        if not confirm:
            return

        # Learn aliases before renaming so future OCR hits also benefit
        for old, new in renames:
            self.store.learn_alias(old, new)
        self.store.save_aliases()

        count = self._apply_renames(renames)
        self.status_var.set(
            f"✓ Renamed {count} field(s) across {len(renames)} variant(s). Data saved."
        )
        self._scan()

    # --------------------------------------------------------------- unmapped
    def _on_unm_select(self, event):
        sel = self.unm_listbox.curselection()
        self.fuzzy_hint.config(text="")
        self.unm_combo.set("")
        if not sel or sel[0] >= len(self.unmapped_labels):
            return
        raw, _ = self.unmapped_labels[sel[0]]

        # Try fuzzy match with a lower cutoff than the main resolver
        norm = normalize_label(raw)
        candidates = difflib.get_close_matches(
            norm,
            [normalize_label(s) for s in LABEL_SUGGESTIONS],
            n=3, cutoff=0.45,
        )
        if candidates:
            # Map back to original spelling
            originals = []
            for c in candidates:
                for s in LABEL_SUGGESTIONS:
                    if normalize_label(s) == c:
                        originals.append(s)
                        break
            self.unm_combo.set(originals[0])
            if len(originals) > 1:
                others = ", ".join(originals[1:])
                self.fuzzy_hint.config(
                    text=f"Other guesses: {others}"
                )
            else:
                self.fuzzy_hint.config(
                    text=f"Best fuzzy guess (low confidence)"
                )
        else:
            self.fuzzy_hint.config(text="No close match found — pick manually")

    def _assign_unmapped(self):
        self._do_assign(learn_alias=False)

    def _assign_and_learn(self):
        self._do_assign(learn_alias=True)

    def _do_assign(self, learn_alias=False):
        sel = self.unm_listbox.curselection()
        if not sel:
            messagebox.showinfo("Assign", "Select a label from the list first.",
                                 parent=self)
            return
        if sel[0] >= len(self.unmapped_labels):
            return
        raw, count = self.unmapped_labels[sel[0]]

        canonical = self.unm_combo.get().strip()
        if not canonical:
            messagebox.showinfo("Assign", "Pick a target canonical name from the dropdown.",
                                 parent=self)
            return

        action = "rename and learn alias for" if learn_alias else "rename"
        confirm = messagebox.askyesno(
            "Confirm",
            f"{action.title()} \"{raw}\" ({count} records) → \"{canonical}\"?\n\nProceed?",
            parent=self,
        )
        if not confirm:
            return

        if learn_alias:
            self.store.learn_alias(raw, canonical)
            self.store.save_aliases()

        n = self._apply_renames([(raw, canonical)])
        self.status_var.set(
            f"✓ Renamed {n} field(s): '{raw}' → '{canonical}'."
            + (" Alias saved." if learn_alias else "")
        )
        self._scan()

#################################################
# Developer tool: Extract variable names        #
#################################################
class ExtractVariablesWindow(tk.Toplevel):
    """Developer tool to extract, view, and copy all unique variable names from the database."""

    def __init__(self, parent, store: DataStore):
        super().__init__(parent)
        self.title("Developer Tool: Extract Variable Names")
        self.geometry("650x550")
        self.minsize(500, 400)
        self.store = store
        self.var_list = []

        self._build_ui()
        self._load_data()

    def _build_ui(self):
        main = ttk.Frame(self, padding=10)
        main.pack(fill=tk.BOTH, expand=True)

        ttk.Label(
            main,
            text="List of all unique 'label' strings currently saved in the records database:",
            font=("", 10, "bold")
        ).pack(anchor="w", pady=(0, 8))

        # Info label
        self.info_var = tk.StringVar(value="Loading...")
        ttk.Label(main, textvariable=self.info_var, foreground="#555555").pack(anchor="w", pady=(0, 4))

        # Treeview
        tree_frame = ttk.Frame(main)
        tree_frame.pack(fill=tk.BOTH, expand=True)

        cols = ("label", "count")
        self.tree = ttk.Treeview(tree_frame, columns=cols, show="headings", height=15)
        self.tree.heading("label", text="Variable Name (Label)")
        self.tree.heading("count", text="Occurrences")
        self.tree.column("label", width=450, anchor="w")
        self.tree.column("count", width=100, anchor="center", stretch=False)
        
        tree_scroll = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=tree_scroll.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Buttons
        btn_frame = ttk.Frame(main)
        btn_frame.pack(fill=tk.X, pady=(12, 0))

        ttk.Button(btn_frame, text="Copy as Python List", command=self._copy_python).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Copy as JSON Array", command=self._copy_json).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Copy Raw (Tab-separated)", command=self._copy_raw).pack(side=tk.LEFT, padx=4)
        
        ttk.Button(btn_frame, text="Export to .txt", command=self._export_txt).pack(side=tk.RIGHT, padx=4)

    def _load_data(self):
        self.var_list = self.store.get_all_variable_names()
        self.tree.delete(*self.tree.get_children())
        
        for label, count in self.var_list:
            self.tree.insert("", tk.END, values=(label, count))
            
        total_records = sum(count for _, count in self.var_list)
        self.info_var.config(text=f"Found {len(self.var_list)} unique variable name(s) across {total_records} total field(s).")

    def _copy_python(self):
        if not self.var_list:
            return
        # Format: ["LABEL1", "LABEL2"]
        py_list = ",\n    ".join(f'"{label}"' for label, _ in self.var_list)
        text = f"[\n    {py_list}\n]"
        self.clipboard_clear()
        self.clipboard_append(text)
        messagebox.showinfo("Copied", "Variable names copied to clipboard as a Python list.", parent=self)

    def _copy_json(self):
        if not self.var_list:
            return
        import json
        json_list = [label for label, _ in self.var_list]
        text = json.dumps(json_list, ensure_ascii=False, indent=2)
        self.clipboard_clear()
        self.clipboard_append(text)
        messagebox.showinfo("Copied", "Variable names copied to clipboard as a JSON array.", parent=self)

    def _copy_raw(self):
        if not self.var_list:
            return
        text = "\n".join(f"{label}\t{count}" for label, count in self.var_list)
        self.clipboard_clear()
        self.clipboard_append(text)
        messagebox.showinfo("Copied", "Variables and counts copied as tab-separated text.", parent=self)

    def _export_txt(self):
        if not self.var_list:
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Export Variable Names",
            parent=self
        )
        if not path:
            return
        with open(path, "w", encoding="utf-8") as f:
            f.write(f"# Extracted Variable Names ({len(self.var_list)} unique)\n")
            f.write(f"# Format: Variable Name <TAB> Occurrence Count\n\n")
            for label, count in self.var_list:
                f.write(f"{label}\t{count}\n")
        messagebox.showinfo("Exported", f"Saved {len(self.var_list)} variables to:\n{path}", parent=self)

class RecordsViewerWindow(tk.Toplevel):
    """Window to browse, search, edit, and delete patient record fields."""

    def __init__(self, parent, store: DataStore):
        super().__init__(parent)
        self.title("Saved Patient Records")
        self.geometry("1100x650")
        self.minsize(800, 450)
        self.store = store
        self.all_rows = [] 
        self._editing_item = None

        self._build_ui()
        self._load_data()

    def _build_ui(self):
        # ---- Search & Action Bar ----
        search_frame = ttk.Frame(self, padding=(10, 8, 10, 4))
        search_frame.pack(fill=tk.X)

        ttk.Label(search_frame, text="Search by Variable:", font=("", 10, "bold")).pack(side=tk.LEFT, padx=(0, 8))
        
        self.search_var = tk.StringVar()
        self.search_var.trace_add("write", lambda *args: self._on_search())
        
        self.search_combo = ttk.Combobox(search_frame, textvariable=self.search_var, width=35)
        self.search_combo.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 8))
        
        ttk.Button(search_frame, text="Clear", command=self._clear_search, width=8).pack(side=tk.LEFT, padx=(0, 6))
        
        # NEW: Delete Button
        self.delete_btn = ttk.Button(search_frame, text="Delete Selected", command=self._delete_selected, width=16)
        self.delete_btn.pack(side=tk.LEFT, padx=(0, 8))

        self.status_var = tk.StringVar(value="Loading data...")
        ttk.Label(search_frame, textvariable=self.status_var, foreground="#555555").pack(side=tk.RIGHT)

        # ---- Treeview ----
        tree_frame = ttk.Frame(self, padding=(10, 4, 10, 10))
        tree_frame.pack(fill=tk.BOTH, expand=True)

        cols = ("patient_id", "patient_name", "visit_date", "label", "value", "unit", "flag")
        self.tree = ttk.Treeview(tree_frame, columns=cols, show="headings", selectmode="extended")
        
        self.tree.heading("patient_id", text="Patient ID")
        self.tree.heading("patient_name", text="Name")
        self.tree.heading("visit_date", text="Visit Date")
        self.tree.heading("label", text="Variable (Label) - Dbl-click to edit")
        self.tree.heading("value", text="Value")
        self.tree.heading("unit", text="Unit")
        self.tree.heading("flag", text="Flag")

        self.tree.column("patient_id", width=80, anchor="center")
        self.tree.column("patient_name", width=120, anchor="w")
        self.tree.column("visit_date", width=100, anchor="center")
        self.tree.column("label", width=200, anchor="w")
        self.tree.column("value", width=80, anchor="center")
        self.tree.column("unit", width=80, anchor="center")
        self.tree.column("flag", width=150, anchor="w")

        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)

        # Bindings for Edit and Delete
        self.tree.bind("<Double-1>", self._on_double_click)
        self.tree.bind("<Delete>", lambda e: self._delete_selected())
        
        # NEW: Right-click context menu for deletion
        self.context_menu = tk.Menu(self, tearoff=0)
        self.context_menu.add_command(label="Delete selected field(s)", command=self._delete_selected)
        self.tree.bind("<Button-3>", self._show_context_menu) # Windows/Linux
        self.tree.bind("<Button-2>", self._show_context_menu) # macOS

    def _load_data(self):
        """Extract all records from the data store into a flat list."""
        self.all_rows = []
        unique_labels = set()

        for pid, p in self.store.data.get("patients", {}).items():
            name = p.get("name", "")
            for v_idx, visit in enumerate(p.get("visits", [])):
                vdate = visit.get("visit_date") or ""
                for f_idx, f in enumerate(visit.get("fields", [])):
                    label = f.get("label", "")
                    if label:
                        unique_labels.add(label)
                    
                    self.all_rows.append({
                        "patient_id": pid,
                        "patient_name": name,
                        "visit_date": vdate,
                        "label": label,
                        "value": f.get("value", ""),
                        "unit": f.get("unit", ""),
                        "flag": f.get("flag", ""),
                        "_v_idx": v_idx,
                        "_f_idx": f_idx
                    })

        self.search_combo["values"] = sorted(unique_labels)
        self._populate_tree(self.all_rows)

    def _populate_tree(self, rows):
        """Clear the tree and insert the provided rows."""
        self.tree.delete(*self.tree.get_children())
        for r in rows:
            iid = f"{r['patient_id']}_{r['_v_idx']}_{r['_f_idx']}"
            self.tree.insert("", tk.END, iid=iid, values=(
                r["patient_id"], r["patient_name"], r["visit_date"],
                r["label"], r["value"], r["unit"], r["flag"]
            ))
        
        total = len(self.all_rows)
        showing = len(rows)
        self.status_var.config(text=f"Showing {showing} of {total} total records")

    def _on_search(self):
        query = self.search_var.get().strip().lower()
        if not query:
            self._populate_tree(self.all_rows)
            return
        filtered_rows = [r for r in self.all_rows if query in r["label"].lower()]
        self._populate_tree(filtered_rows)

    def _clear_search(self):
        self.search_var.set("")
        self.tree.focus_set()

    # ------------------------------------------------------------------
    # RIGHT-CLICK CONTEXT MENU
    # ------------------------------------------------------------------
    def _show_context_menu(self, event):
        item_id = self.tree.identify_row(event.y)
        if item_id:
            # If right-clicked on an unselected row, select just that row first
            if item_id not in self.tree.selection():
                self.tree.selection_set(item_id)
            self.context_menu.tk_popup(event.x_root, event.y_root)

    # ------------------------------------------------------------------
    # DELETE LOGIC
    # ------------------------------------------------------------------
    def _delete_selected(self):
        selected_items = self.tree.selection()
        if not selected_items:
            messagebox.showinfo("Delete", "Please select one or more rows to delete.", parent=self)
            return

        # Prepare confirmation message
        if len(selected_items) == 1:
            vals = self.tree.item(selected_items[0], "values")
            msg = (f"Are you sure you want to delete this field?\n\n"
                   f"Patient: {vals[1]}\nDate: {vals[2]}\n"
                   f"Field: {vals[3]} ({vals[4]} {vals[5]})")
        else:
            msg = f"Are you sure you want to delete {len(selected_items)} selected field(s)?"

        if not messagebox.askyesno("Confirm Deletion", msg, parent=self):
            return

        # Parse IDs to locate them in the JSON structure
        deletes = []
        for item_id in selected_items:
            parts = item_id.split("_")
            pid = parts[0]
            v_idx = int(parts[-2])
            f_idx = int(parts[-1])
            deletes.append((pid, v_idx, f_idx, item_id))

        # CRITICAL: Sort by field index in DESCENDING order. 
        # If we delete index 0, index 1 becomes 0. Sorting descending 
        # ensures we delete from the end of the list first, preventing index shifting bugs.
        deletes.sort(key=lambda x: x[2], reverse=True)

        for pid, v_idx, f_idx, item_id in deletes:
            try:
                # Remove from database dictionary
                self.store.data["patients"][pid]["visits"][v_idx]["fields"].pop(f_idx)
                
                # Remove from internal flat list
                self.all_rows = [r for r in self.all_rows if not (
                    r["patient_id"] == pid and r["_v_idx"] == v_idx and r["_f_idx"] == f_idx
                )]
                
                # Remove from UI Treeview
                self.tree.delete(item_id)
            except (IndexError, KeyError):
                continue  # Skip if already deleted or data is out of sync

        # Save database once after all deletions
        self.store.save()
        
        # Update status bar
        showing = len(self.tree.get_children())
        total = len(self.all_rows)
        self.status_var.config(text=f"Showing {showing} of {total} total records")

    # ------------------------------------------------------------------
    # INLINE EDITING LOGIC
    # ------------------------------------------------------------------
    def _on_double_click(self, event):
        col = self.tree.identify_column(event.x)
        if col != "#4":  # Only edit the Label column
            return

        item_id = self.tree.identify_row(event.y)
        if not item_id:
            return

        bbox = self.tree.bbox(item_id, "label")
        if not bbox:
            return

        x, y, w, h = bbox
        current_val = self.tree.item(item_id, "values")[3]

        self.edit_entry = tk.Entry(self.tree, width=30)
        self.edit_entry.place(x=x, y=y, width=w, height=h + 2)
        self.edit_entry.insert(0, current_val)
        self.edit_entry.select_range(0, tk.END)
        self.edit_entry.focus_set()

        self._editing_item = item_id

        self.edit_entry.bind("<Return>", self._save_edit)
        self.edit_entry.bind("<Escape>", self._cancel_edit)
        self.edit_entry.bind("<FocusOut>", self._save_edit)

    def _save_edit(self, event=None):
        if not self._editing_item or not hasattr(self, 'edit_entry'):
            return

        new_val = self.edit_entry.get().strip()
        item_id = self._editing_item

        self.edit_entry.destroy()
        del self.edit_entry
        self._editing_item = None

        if not new_val:
            return

        try:
            parts = item_id.split("_")
            pid = parts[0]
            v_idx = int(parts[-2])
            f_idx = int(parts[-1])
        except (ValueError, IndexError):
            return

        patient = self.store.data.get("patients", {}).get(pid)
        if not patient:
            return

        visit = patient["visits"][v_idx]
        old_label = visit["fields"][f_idx].get("label", "")

        if old_label != new_val:
            visit["fields"][f_idx]["label"] = new_val
            self.store.save()

            for row in self.all_rows:
                if (row["patient_id"] == pid and 
                    row["_v_idx"] == v_idx and 
                    row["_f_idx"] == f_idx):
                    row["label"] = new_val
                    break

            self._on_search() # Refresh view

    def _cancel_edit(self, event=None):
        if hasattr(self, 'edit_entry'):
            self.edit_entry.destroy()
            del self.edit_entry
        self._editing_item = None

class App(tk.Tk):
    CONF_GOOD = 0.85
    CONF_OK = 0.60

    def __init__(self):
        super().__init__()
        self.title("Medical Report Scanner")
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        
        window_width = 800
        window_height = 600

        x = (screen_width // 2) - (window_width // 2)
        y = (screen_height // 2) - (window_height // 2)
        
        self.geometry(f"{window_width}x{window_height}+{x}+{y}")

        self.ocr_engine = OCREngine()
        self.store = None
        self.image_paths = []
        self.current_index = -1
        self.pipeline = None
        self.display_scale = 1.0
        self.fit_scale = 1.0
        self.zoom_factor = 1.0
        self.disp_w = 1
        self.disp_h = 1
        self.tk_image = None
        self.ocr_lines = []
        self.grouped_lines = []
        self.field_rows = []
        self.crop_start = None
        self.crop_rect_id = None
        self.ocr_queue = queue.Queue()
        self.current_visit_index = None
        self.record_index_map = []
        self.pending_preselect = None

        self.angle_var = tk.DoubleVar(value=0.0)
        self.denoise_var = tk.BooleanVar(value=False)
        self.contrast_var = tk.BooleanVar(value=False)
        self.shadow_var = tk.BooleanVar(value=False)
        self.grid_var = tk.BooleanVar(value=False)

        self._build_menu_bar()
        self._build_ui()

    def _build_menu_bar(self):
        menubar = tk.Menu(self, tearoff=0)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Select Folder...", command=self.select_folder,
                               accelerator="Ctrl+O")
        file_menu.add_command(label="View Saved Records...", command=self.open_records_viewer,
                               accelerator="Ctrl+H")
        file_menu.add_command(label="Manage Synonyms...", command=self.open_alias_manager)
        file_menu.add_separator()
        file_menu.add_command(label="Fix Variable Names…", command=self.open_fix_variables,
                               accelerator="Ctrl+Shift+F")        
        file_menu.add_command(label="Exploratory Data Analysis…", command=self.open_eda,
                               accelerator="Ctrl+E")
        file_menu.add_command(label="Exit", command=self.destroy, accelerator="Ctrl+Q")
        file_menu.add_separator()
        file_menu.add_command(label="Developer: Extract Variable Names…", command=self.open_extract_variables)
        menubar.add_cascade(label="File", menu=file_menu)

        image_menu = tk.Menu(menubar, tearoff=0)
        image_menu.add_command(label="Scan", command=self.auto_scan_current,
                                accelerator="Ctrl+K")
        image_menu.add_separator()
        image_menu.add_command(label="Rotate 90° Clockwise", command=lambda: self.rotate(90))
        image_menu.add_command(label="Rotate 90° Counter-clockwise", command=lambda: self.rotate(-90))
        image_menu.add_separator()
        image_menu.add_checkbutton(label="Denoise", variable=self.denoise_var,
                                    command=self.refresh_display)
        image_menu.add_checkbutton(label="Contrast (CLAHE)", variable=self.contrast_var,
                                    command=self.refresh_display)
        image_menu.add_checkbutton(label="Remove Shadows", variable=self.shadow_var,
                                    command=self.refresh_display)
        image_menu.add_checkbutton(label="Show Grid", variable=self.grid_var,
                                    command=self.refresh_display)
        image_menu.add_separator()
        image_menu.add_command(label="Clear Crop", command=self.clear_crop)
        image_menu.add_command(label="Reset All Adjustments", command=self.reset_pipeline)
        menubar.add_cascade(label="Image", menu=image_menu)

        view_menu = tk.Menu(menubar, tearoff=0)
        view_menu.add_command(label="Zoom In", command=lambda: self.zoom_step(1),
                               accelerator="Ctrl++")
        view_menu.add_command(label="Zoom Out", command=lambda: self.zoom_step(-1),
                               accelerator="Ctrl+-")
        view_menu.add_command(label="Fit to Window", command=self.zoom_fit,
                               accelerator="Ctrl+0")
        menubar.add_cascade(label="View", menu=view_menu)

        ocr_menu = tk.Menu(menubar, tearoff=0)
        ocr_menu.add_command(label="Run OCR", command=self.run_ocr, accelerator="Ctrl+R")
        ocr_menu.add_command(label="Force Reprocess Current Image", command=self.force_reprocess)
        menubar.add_cascade(label="OCR", menu=ocr_menu)

        visit_menu = tk.Menu(menubar, tearoff=0)
        visit_menu.add_command(label="Add Blank Field", command=self.add_blank_field)
        visit_menu.add_command(label="Save Visit && Next Image", command=self.save_and_next,
                                accelerator="Ctrl+S")
        visit_menu.add_command(label="Skip Image (No Save)", command=self.next_image,
                                accelerator="Ctrl+Right")
        menubar.add_cascade(label="Visit", menu=visit_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self.show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.config(menu=menubar)
        self.menubar = menubar

        self.bind_all("<Control-o>", lambda e: self.select_folder())
        self.bind_all("<Control-h>", lambda e: self.open_records_viewer())
        self.bind_all("<Control-q>", lambda e: self.destroy())
        self.bind_all("<Control-k>", lambda e: self.auto_scan_current())
        self.bind_all("<Control-r>", lambda e: self.run_ocr())
        self.bind_all("<Control-s>", lambda e: self.save_and_next())
        self.bind_all("<Control-Right>", lambda e: self.next_image())
        self.bind_all("<Control-plus>", lambda e: self.zoom_step(1))
        self.bind_all("<Control-KP_Add>", lambda e: self.zoom_step(1))
        self.bind_all("<Control-minus>", lambda e: self.zoom_step(-1))
        self.bind_all("<Control-KP_Subtract>", lambda e: self.zoom_step(-1))
        self.bind_all("<Control-0>", lambda e: self.zoom_fit())
        self.bind_all("<Control-e>", lambda e: self.open_eda())
        self.bind_all("<Control-Shift-F>", lambda e: self.open_fix_variables())        

    def show_about(self):
        messagebox.showinfo(
            "About",
            "Diabetes Records Digitizer\n\n"
            "Scan, OCR, and review FNS follow-up forms into structured "
            "per-patient records.\n\n"
            "All actions are available from the menu bar above: File, "
            "Image, View, OCR, and Visit."
        )

    def _build_ui(self):
        top = ttk.Frame(self)
        top.pack(side=tk.TOP, fill=tk.X, padx=6, pady=4)
        self.progress_label = ttk.Label(top, text="No folder selected")
        self.progress_label.pack(side=tk.LEFT, padx=12)
        self.status_label = ttk.Label(top, text="", foreground="#555555")
        self.status_label.pack(side=tk.RIGHT)

        main = ttk.Frame(self)
        main.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=6, pady=4)

        right_outer = ttk.Frame(main, width=540)
        right_outer.pack(side=tk.RIGHT, fill=tk.BOTH)
        right_outer.pack_propagate(False)

        right_vsb = ttk.Scrollbar(right_outer, orient="vertical")
        right_vsb.pack(side=tk.RIGHT, fill=tk.Y)

        self.right_canvas = tk.Canvas(right_outer, borderwidth=0, highlightthickness=0)
        self.right_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.right_canvas.configure(yscrollcommand=right_vsb.set)
        right_vsb.config(command=self.right_canvas.yview)

        right = ttk.Frame(self.right_canvas)
        self.right_canvas_window = self.right_canvas.create_window((0, 0), window=right, anchor="nw")

        def _on_right_configure(e):
            self.right_canvas.configure(scrollregion=self.right_canvas.bbox("all"))
        right.bind("<Configure>", _on_right_configure)

        def _on_right_canvas_configure(e):
            self.right_canvas.itemconfig(self.right_canvas_window, width=e.width)
        self.right_canvas.bind("<Configure>", _on_right_canvas_configure)

        left = ttk.Frame(main)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        controls = ttk.Frame(left)
        controls.pack(side=tk.TOP, fill=tk.X, pady=4)
        ttk.Label(controls, text="Fine angle:").pack(side=tk.LEFT, padx=(0, 2))
        angle_scale = ttk.Scale(controls, from_=-15, to=15, variable=self.angle_var,
                                 orient=tk.HORIZONTAL, length=150,
                                 command=lambda e: self.on_angle_change())
        angle_scale.pack(side=tk.LEFT)

        ttk.Label(controls, text="Zoom:").pack(side=tk.LEFT, padx=(12, 2))
        self.zoom_label = ttk.Label(controls, text="100%", width=5, anchor="center")
        self.zoom_label.pack(side=tk.LEFT)

        ttk.Label(
            controls,
            text="(drag = crop, wheel = zoom, middle-drag = pan — Scan, Rotate, Denoise, "
                 "Contrast, Shadows, Grid, and Zoom are all in the Image/View menus)"
        ).pack(side=tk.LEFT, padx=8)

        canvas_frame = ttk.Frame(left)
        canvas_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        canvas_frame.rowconfigure(0, weight=1)
        canvas_frame.columnconfigure(0, weight=1)

        self.canvas = tk.Canvas(canvas_frame, bg="#222222", cursor="cross")
        self.canvas.grid(row=0, column=0, sticky="nsew")
        vbar = ttk.Scrollbar(canvas_frame, orient="vertical", command=self.canvas.yview)
        vbar.grid(row=0, column=1, sticky="ns")
        hbar = ttk.Scrollbar(canvas_frame, orient="horizontal", command=self.canvas.xview)
        hbar.grid(row=1, column=0, sticky="ew")
        self.canvas.configure(yscrollcommand=vbar.set, xscrollcommand=hbar.set)

        self.canvas.bind("<ButtonPress-1>", self.on_crop_start)
        self.canvas.bind("<B1-Motion>", self.on_crop_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_crop_end)

        self.canvas.bind("<MouseWheel>", self.on_mousewheel_zoom)
        self.canvas.bind("<Button-4>", self.on_mousewheel_zoom)
        self.canvas.bind("<Button-5>", self.on_mousewheel_zoom)

        self.canvas.bind("<ButtonPress-2>", lambda e: self.canvas.scan_mark(e.x, e.y))
        self.canvas.bind("<B2-Motion>", lambda e: self.canvas.scan_dragto(e.x, e.y, gain=1))

        self.bind_all("<MouseWheel>", self.on_global_mousewheel)
        self.bind_all("<Button-4>", self.on_global_mousewheel)
        self.bind_all("<Button-5>", self.on_global_mousewheel)

        patient_box = ttk.LabelFrame(right, text="Patient & visit")
        patient_box.pack(side=tk.TOP, fill=tk.X, padx=4, pady=4)

        ttk.Label(patient_box, text="Patient ID:").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self.patient_id_var = tk.StringVar()
        self.patient_combo = ttk.Combobox(patient_box, textvariable=self.patient_id_var, values=[])
        self.patient_combo.grid(row=0, column=1, sticky="we", padx=4, pady=2)
        self.patient_combo.bind("<<ComboboxSelected>>", self.on_patient_selected)
        self.patient_combo.bind("<FocusOut>", self.on_patient_selected)
        self.patient_combo.bind("<Return>", self.on_patient_selected)

        ttk.Label(patient_box, text="Patient name:").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.patient_name_var = tk.StringVar()
        ttk.Entry(patient_box, textvariable=self.patient_name_var).grid(row=1, column=1, sticky="we", padx=4, pady=2)

        ttk.Label(patient_box, text="Record (visit):").grid(row=2, column=0, sticky="w", padx=4, pady=2)
        self.record_var = tk.StringVar()
        self.record_combo = ttk.Combobox(patient_box, textvariable=self.record_var,
                                          values=["+ New record"], state="readonly")
        self.record_combo.set("+ New record")
        self.record_combo.grid(row=2, column=1, sticky="we", padx=4, pady=2)
        self.record_combo.bind("<<ComboboxSelected>>", self.on_record_selected)
        ttk.Label(patient_box, text="(pick an existing record to add this image to it)",
                  foreground="#777777").grid(row=3, column=0, columnspan=2, sticky="w", padx=4)

        ttk.Label(patient_box, text="Visit date:").grid(row=4, column=0, sticky="w", padx=4, pady=2)
        self.visit_date_var = tk.StringVar()
        ttk.Entry(patient_box, textvariable=self.visit_date_var).grid(row=4, column=1, sticky="we", padx=4, pady=2)
        ttk.Label(patient_box, text="(YYYY-MM-DD, auto-suggested from OCR)",
                  foreground="#777777").grid(row=5, column=0, columnspan=2, sticky="w", padx=4)
        patient_box.columnconfigure(1, weight=1)

        review_label = ttk.Label(right, text="Detected fields (edit / label / keep):")
        review_label.pack(side=tk.TOP, anchor="w", padx=4, pady=(8, 0))

        review_container = ttk.Frame(right)
        review_container.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.review_canvas = tk.Canvas(review_container, borderwidth=0, height=400)
        vsb = ttk.Scrollbar(review_container, orient="vertical", command=self.review_canvas.yview)
        self.review_frame = ttk.Frame(self.review_canvas)
        self.review_frame.bind(
            "<Configure>",
            lambda e: self.review_canvas.configure(scrollregion=self.review_canvas.bbox("all"))
        )
        self.review_canvas.create_window((0, 0), window=self.review_frame, anchor="nw")
        self.review_canvas.configure(yscrollcommand=vsb.set)
        self.review_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

        row_header = ttk.Frame(self.review_frame)
        row_header.pack(side=tk.TOP, fill=tk.X)
        for text, w in [("Keep", 4), ("Label", 22), ("Value", 12), ("Unit", 12), ("Conf", 5)]:
            ttk.Label(row_header, text=text, width=w, font=("TkDefaultFont", 9, "bold")).pack(side=tk.LEFT, padx=1)

        bottom_right = ttk.Frame(right)
        bottom_right.pack(side=tk.BOTTOM, fill=tk.X, padx=4, pady=6)
        ttk.Label(
            bottom_right,
            text="Use the Visit menu to add a blank field, save this visit, or skip this image.",
            foreground="#777777"
        ).pack(side=tk.LEFT)

    def on_global_mousewheel(self, event):
        try:
            widget = self.winfo_containing(event.x_root, event.y_root)
        except KeyError:
            return
        if not widget:
            return

        if event.num == 4:
            delta = 1
        elif event.num == 5:
            delta = -1
        else:
            delta = int(event.delta / 120)

        widget_str = str(widget)

        if widget_str.startswith(str(self.canvas)):
            cx = event.x_root - self.canvas.winfo_rootx()
            cy = event.y_root - self.canvas.winfo_rooty()
            direction = 1 if delta > 0 else -1
            self._zoom_at(cx, cy, direction)
        elif widget_str.startswith(str(self.review_canvas)):
            self.review_canvas.yview_scroll(-1 * delta, "units")
        elif widget_str.startswith(str(self.right_canvas)) or any(widget_str.startswith(str(child)) for child in self.right_canvas.children.values()):
            self.right_canvas.yview_scroll(-1 * delta, "units")

    def select_folder(self):
        folder = filedialog.askdirectory(title="Select folder with patient photos")
        if not folder:
            return
        self.store = DataStore(folder)
        all_images = []
        for root, _dirs, files in os.walk(folder):
            if os.path.basename(root) == "patients":
                continue
            for fn in files:
                if os.path.splitext(fn)[1].lower() in IMAGE_EXTENSIONS:
                    all_images.append(os.path.join(root, fn))
        all_images.sort()
        self.image_paths = all_images
        self.patient_combo["values"] = self.store.patient_ids()
        self.current_index = -1
        if not all_images:
            messagebox.showinfo("No images found", "No image files were found in that folder (or subfolders).")
            return
        self.next_image()

    def next_image(self):
        self.current_index += 1
        while (self.current_index < len(self.image_paths)
               and self.store.is_processed(self.image_paths[self.current_index])):
            self.current_index += 1
        if self.current_index >= len(self.image_paths):
            self.progress_label.config(text=f"All {len(self.image_paths)} images processed.")
            messagebox.showinfo("Done", "All images in this folder have been processed.")
            return
        self.load_current_image()

    def force_reprocess(self):
        if self.current_index < 0 or self.current_index >= len(self.image_paths):
            return
        path = self.image_paths[self.current_index]
        h = file_hash(path)
        entry = self.store.manifest["processed"].get(h)
        if entry:
            self.pending_preselect = (entry["patient_id"], entry["visit_index"])
        self.store.manifest["processed"].pop(h, None)
        self.load_current_image()

    def load_current_image(self):
        path = self.image_paths[self.current_index]
        try:
            self.pipeline = ImagePipeline(path)
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            self.next_image()
            return
        self.angle_var.set(0.0)
        self.denoise_var.set(False)
        self.contrast_var.set(False)
        self.grid_var.set(False)
        self.shadow_var.set(False)
        self.zoom_factor = 1.0
        self.ocr_lines = []
        self.grouped_lines = []
        self._clear_review_rows()
        self.current_visit_index = None
        self.patient_id_var.set("")
        self.patient_name_var.set("")
        self.visit_date_var.set("")
        self.record_combo["values"] = ["+ New record"]
        self.record_combo.set("+ New record")
        total = len(self.image_paths)
        done = len(self.store.manifest["processed"])
        self.progress_label.config(
            text=f"Image {self.current_index + 1} / {total}  |  {done} already saved  |  {os.path.basename(path)}"
        )
        self.refresh_display()

        if self.pending_preselect:
            pid, visit_index = self.pending_preselect
            self.pending_preselect = None
            self.patient_id_var.set(pid)
            self.patient_name_var.set(self.store.patient_name(pid))
            self.refresh_record_options()
            for combo_idx, v_idx in enumerate(self.record_index_map):
                if v_idx == visit_index:
                    self.record_combo.current(combo_idx + 1)
                    self.on_record_selected()
                    break

    def refresh_record_options(self):
        pid = self.patient_id_var.get().strip()
        options = self.store.record_options(pid) if (self.store and pid) else []
        self.record_index_map = [idx for idx, _label in options]
        labels = ["+ New record"] + [label for _idx, label in options]
        self.record_combo["values"] = labels
        self.record_combo.set("+ New record")
        self.current_visit_index = None

    def rotate(self, delta):
        if not self.pipeline:
            return
        self.pipeline.rotation_90 = (self.pipeline.rotation_90 + delta) % 360
        self.refresh_display()

    def on_angle_change(self):
        if not self.pipeline:
            return
        self.pipeline.fine_angle = self.angle_var.get()
        self.refresh_display()

    def clear_crop(self):
        if not self.pipeline:
            return
        self.pipeline.crop_rect = None
        self.refresh_display()

    def reset_pipeline(self):
        if not self.pipeline:
            return
        self.pipeline.reset()
        self.angle_var.set(0.0)
        self.denoise_var.set(False)
        self.contrast_var.set(False)
        self.grid_var.set(False)
        self.shadow_var.set(False)
        self.refresh_display()

    def auto_scan_current(self):
        if not self.pipeline:
            return
        self.status_label.config(text="Scanning...")
        self.update_idletasks()

        found = self.pipeline.auto_scan()

        self.angle_var.set(0.0)
        self.contrast_var.set(True)
        self.shadow_var.set(True)
        self.zoom_factor = 1.0

        self.ocr_lines = []
        self.grouped_lines = []

        self.refresh_display()  
        self.zoom_fit()

        if found:
            self.status_label.config(
                text="Scan complete: document detected, perspective-corrected, "
                     "auto-rotated, deskewed, and shadows removed."
            )
        else:
            self.status_label.config(
                text="Scan complete: no clear document edges found, so the full "
                     "frame was used (still auto-rotated, deskewed, and shadow-corrected). "
                     "Use Crop/Rotate manually if needed."
            )

    def refresh_display(self):
        if not self.pipeline:
            return
        self.pipeline.denoise = self.denoise_var.get()
        self.pipeline.contrast = self.contrast_var.get()
        self.pipeline.shadow_removal = self.shadow_var.get()
        img = self.pipeline.render()
        self._show_image(img, boxes=self.ocr_lines, row_boxes=self.grouped_lines)

    def _show_image(self, bgr_img, boxes=None, row_boxes=None, keep_view=False):
        self.canvas.update_idletasks()
        cw = max(self.canvas.winfo_width(), 400)
        ch = max(self.canvas.winfo_height(), 400)
        h, w = bgr_img.shape[:2]
        self.fit_scale = max(min(cw / w, ch / h), 0.02)
        scale = max(self.fit_scale * self.zoom_factor, 0.02)
        self.display_scale = scale

        prev_xview = self.canvas.xview() if keep_view else None
        prev_yview = self.canvas.yview() if keep_view else None

        disp_w, disp_h = max(1, int(w * scale)), max(1, int(h * scale))
        self.disp_w, self.disp_h = disp_w, disp_h
        interp = cv2.INTER_AREA if scale <= 1.0 else cv2.INTER_LINEAR
        resized = cv2.resize(bgr_img, (disp_w, disp_h), interpolation=interp)
        rgb = cv2.cvtColor(resized, cv2.COLOR_BGR2RGB)
        pil_img = Image.fromarray(rgb)
        self.tk_image = ImageTk.PhotoImage(pil_img)
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor="nw", image=self.tk_image, tags="img")
        self.canvas.config(scrollregion=(0, 0, disp_w, disp_h))

        if self.grid_var.get():
            step = 40
            for x in range(0, disp_w, step):
                self.canvas.create_line(x, 0, x, disp_h, fill="#444444", width=1, tags="grid")
            for y in range(0, disp_h, step):
                self.canvas.create_line(0, y, disp_w, y, fill="#444444", width=1, tags="grid")

        if boxes:
            for i, line in enumerate(boxes):
                conf = line["confidence"]
                color = "#2ecc71" if conf >= self.CONF_GOOD else ("#f1c40f" if conf >= self.CONF_OK else "#e74c3c")
                pts = []
                for x, y in line["bbox"]:
                    pts.extend([x * scale, y * scale])
                self.canvas.create_polygon(*pts, outline=color, fill="", width=2, tags=(f"box_{i}",))

        if row_boxes:
            for i, row in enumerate(row_boxes):
                pts = []
                for x, y in row["bbox"]:
                    pts.extend([x * scale, y * scale])
                self.canvas.create_polygon(*pts, outline="#3498db", fill="", width=1,
                                            dash=(4, 2), tags=(f"row_{i}",))

        if prev_xview is not None:
            self.canvas.xview_moveto(prev_xview[0])
            self.canvas.yview_moveto(prev_yview[0])

        self.zoom_label.config(text=f"{int(self.zoom_factor * 100)}%")

    def zoom_step(self, direction):
        if not self.pipeline:
            return
        cx = self.canvas.winfo_width() / 2
        cy = self.canvas.winfo_height() / 2
        self._zoom_at(cx, cy, direction)

    def zoom_fit(self):
        if not self.pipeline:
            return
        self.zoom_factor = 1.0
        self.refresh_display()
        self.canvas.xview_moveto(0)
        self.canvas.yview_moveto(0)

    def on_mousewheel_zoom(self, event):
        pass

    def _zoom_at(self, widget_x, widget_y, direction):
        old_scale = self.display_scale
        canvas_x = self.canvas.canvasx(widget_x)
        canvas_y = self.canvas.canvasy(widget_y)
        img_x = canvas_x / old_scale
        img_y = canvas_y / old_scale

        step = 1.25 if direction > 0 else (1 / 1.25)
        new_zoom = self.zoom_factor * step
        new_zoom = max(0.2, min(new_zoom, 8.0))
        if abs(new_zoom - self.zoom_factor) < 1e-3:
            return
        self.zoom_factor = new_zoom

        self.refresh_display()

        new_scale = self.display_scale
        new_canvas_x = img_x * new_scale
        new_canvas_y = img_y * new_scale
        frac_x = (new_canvas_x - widget_x) / max(self.disp_w, 1)
        frac_y = (new_canvas_y - widget_y) / max(self.disp_h, 1)
        frac_x = min(max(frac_x, 0.0), 1.0)
        frac_y = min(max(frac_y, 0.0), 1.0)
        self.canvas.xview_moveto(frac_x)
        self.canvas.yview_moveto(frac_y)

    def on_crop_start(self, event):
        self.crop_start = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
        if self.crop_rect_id:
            self.canvas.delete(self.crop_rect_id)
            self.crop_rect_id = None

    def on_crop_drag(self, event):
        if not self.crop_start:
            return
        if self.crop_rect_id:
            self.canvas.delete(self.crop_rect_id)
        x0, y0 = self.crop_start
        cx, cy = self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)
        self.crop_rect_id = self.canvas.create_rectangle(x0, y0, cx, cy, outline="#00d4ff", width=2)

    def on_crop_end(self, event):
        if not self.crop_start or not self.pipeline:
            return
        x0, y0 = self.crop_start
        x1, y1 = self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)
        self.crop_start = None
        if abs(x1 - x0) < 8 or abs(y1 - y0) < 8:
            return
        scale = self.display_scale or 1.0
        rect = (int(x0 / scale), int(y0 / scale), int(x1 / scale), int(y1 / scale))
        if self.pipeline.crop_rect:
            ox1, oy1, _ox2, _oy2 = self.pipeline.crop_rect
            rect = (rect[0] + ox1, rect[1] + oy1, rect[2] + ox1, rect[3] + oy1)
        self.pipeline.crop_rect = rect
        self.refresh_display()

    def run_ocr(self):
        if not self.pipeline:
            return
        self.status_label.config(text="Running OCR... (first run downloads models, please wait)")
        self.update_idletasks()
        img = self.pipeline.render()

        def worker():
            try:
                lines = self.ocr_engine.run(img)
                self.ocr_queue.put(("ok", lines))
            except Exception as e:
                self.ocr_queue.put(("error", str(e)))

        threading.Thread(target=worker, daemon=True).start()
        self.after(200, self._poll_ocr_queue)

    def _poll_ocr_queue(self):
        try:
            status, payload = self.ocr_queue.get_nowait()
        except queue.Empty:
            self.after(200, self._poll_ocr_queue)
            return
        if status == "error":
            self.status_label.config(text="OCR failed")
            messagebox.showerror("OCR error", payload)
            return
        self.ocr_lines = payload
        self.grouped_lines = group_ocr_lines(payload)
        self.status_label.config(
            text=f"OCR found {len(payload)} boxes, reconstructed into {len(self.grouped_lines)} line(s)"
        )
        self._show_image(self.pipeline.render(), boxes=self.ocr_lines, row_boxes=self.grouped_lines)
        self._clear_review_rows(origin="ocr")
        for row_line in self.grouped_lines:
            raw_label, value, unit = parse_label_value_unit(row_line["text"])
            extra_aliases = self.store.aliases if self.store else None
            label = match_known_label(raw_label, extra_aliases=extra_aliases)
            self._add_review_row(
                text=value or row_line["text"], confidence=row_line["confidence"],
                origin="ocr", label=label, unit=unit, raw_label=raw_label,
            )
        self._suggest_date()

    def _suggest_date(self):
        if self.visit_date_var.get():
            return
        for line in self.grouped_lines:
            d = try_parse_date(line["text"])
            if d:
                self.visit_date_var.set(d)
                return
        for line in self.ocr_lines:
            d = try_parse_date(line["text"])
            if d:
                self.visit_date_var.set(d)
                break

    def _clear_review_rows(self, origin=None):
        keep_rows = []
        for row in self.field_rows:
            if origin is None or row["origin"] == origin:
                row["frame"].destroy()
            else:
                keep_rows.append(row)
        self.field_rows = keep_rows

    def add_blank_field(self):
        self._add_review_row(text="", confidence=None, origin="manual")

    def _add_review_row(self, text, confidence, origin="ocr", label="", unit="", raw_label=""):
        frame = ttk.Frame(self.review_frame)
        frame.pack(side=tk.TOP, fill=tk.X, pady=1)

        keep_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frame, variable=keep_var, width=4).pack(side=tk.LEFT, padx=1)

        label_var = tk.StringVar(value=label)
        label_combo = ttk.Combobox(
            frame, 
            textvariable=label_var, 
            values=LABEL_SUGGESTIONS, 
            height=8, 
            width=20
        )
        label_combo.pack(side=tk.LEFT, padx=1)

        value_var = tk.StringVar(value=text)
        value_entry = ttk.Entry(frame, textvariable=value_var, width=12)
        value_entry.pack(side=tk.LEFT, padx=1)

        unit_var = tk.StringVar(value=unit)
        unit_entry = ttk.Entry(frame, textvariable=unit_var, width=12)
        unit_entry.pack(side=tk.LEFT, padx=1)

        conf_text = f"{confidence:.2f}" if confidence is not None else "--"
        conf_color = "#2ecc71" if (confidence or 0) >= self.CONF_GOOD else (
            "#f1c40f" if (confidence or 0) >= self.CONF_OK else "#e74c3c")
        conf_label = ttk.Label(frame, text=conf_text, width=5, foreground=conf_color)
        conf_label.pack(side=tk.LEFT, padx=1)

        remove_btn = ttk.Button(frame, text="x", width=2,
                                 command=lambda: self._remove_row(frame))
        remove_btn.pack(side=tk.LEFT, padx=1)

        tag_text = "(existing)" if origin == "existing" else ""
        warn_label = ttk.Label(frame, text=tag_text, foreground="#7f8c8d")
        warn_label.pack(side=tk.LEFT, padx=4)

        def on_change(*_args):
            warning = validate_value(label_var.get(), value_var.get())
            warn_label.config(text=warning or tag_text, foreground="#e67e22" if warning else "#7f8c8d")

        label_var.trace_add("write", on_change)
        value_var.trace_add("write", on_change)

        self.field_rows.append({
            "frame": frame, "keep": keep_var, "label": label_var,
            "value": value_var, "unit": unit_var, "confidence": confidence,
            "origin": origin, "raw_label": raw_label,
        })

    def _remove_row(self, frame):
        self.field_rows = [r for r in self.field_rows if r["frame"] is not frame]
        frame.destroy()

    def on_patient_selected(self, _event=None):
        pid = self.patient_id_var.get().strip()
        if pid and self.store:
            self.patient_name_var.set(self.store.patient_name(pid))
        self.refresh_record_options()

    def on_record_selected(self, _event=None):
        choice = self.record_combo.current()
        if choice <= 0:
            self.current_visit_index = None
            self._clear_review_rows(origin="existing")
            return
        visit_index = self.record_index_map[choice - 1]
        pid = self.patient_id_var.get().strip()
        visit = self.store.get_visit(pid, visit_index)
        if not visit:
            return
        self.current_visit_index = visit_index
        self.visit_date_var.set(visit.get("visit_date") or "")
        self._clear_review_rows(origin="existing")
        existing_fields = visit.get("fields", [])
        for field in existing_fields:
            self._add_review_row(
                text=field.get("value", ""),
                confidence=field.get("confidence"),
                origin="existing",
                label=field.get("label", ""),
                unit=field.get("unit", ""),
            )
        n_images = len(visit.get("images", []))
        self.status_label.config(
            text=f"Loaded record with {len(existing_fields)} existing field(s) from {n_images} image(s). "
                 "This image will be added to it."
        )

    def save_and_next(self):
        if not self.store or not self.pipeline:
            return
        pid = self.patient_id_var.get().strip()
        if not pid:
            messagebox.showwarning("Missing patient", "Please enter or select a Patient ID before saving.")
            return
        name = self.patient_name_var.get().strip()
        visit_date = self.visit_date_var.get().strip()

        fields = []
        learned_count = 0
        for row in self.field_rows:
            if not row["keep"].get():
                continue
            label = row["label"].get().strip()
            value = row["value"].get().strip()
            unit = row["unit"].get().strip()
            if not label and not value:
                continue

            raw_label = row.get("raw_label", "")
            if raw_label and normalize_label(raw_label) != normalize_label(label):
                canonical = _LABEL_LOOKUP_NORM.get(normalize_label(label))
                if canonical and self.store.learn_alias(raw_label, canonical):
                    learned_count += 1

            warning = validate_value(label, value)
            spec = KNOWN_FIELDS.get(normalize_label(label))
            fields.append({
                "label": label,
                "value": value,
                "confidence": row["confidence"],
                "unit": unit or (spec["unit"] if spec else ""),
                "flag": warning or "",
            })

        self.store.save_or_update_visit(
            pid, name, self.current_visit_index,
            self.image_paths[self.current_index], visit_date, fields
        )
        self.store.save()
        self.patient_combo["values"] = self.store.patient_ids()
        action = "Updated" if self.current_visit_index is not None else "Saved new"
        msg = f"{action} record for patient {pid} ({len(fields)} fields)"
        if learned_count:
            msg += f" - learned {learned_count} new label synonym(s)"
        self.status_label.config(text=msg)
        self.next_image()

    def open_fix_variables(self):
        if not self.store:
            messagebox.showinfo(
                "Fix Variable Names",
                "Please select a folder first (File → Select Folder)."
            )
            return
        total_fields = sum(
            len(f.get("fields", []))
            for p in self.store.data.get("patients", {}).values()
            for f in p.get("visits", [])
        )
        if total_fields == 0:
            messagebox.showinfo(
                "Fix Variable Names",
                "No extracted fields found in the current data store.\n\n"
                "Run OCR and save some visits first."
            )
            return
        FixVariableNamesWindow(self, self.store)

    ################################
    # Developer / advanced features
    ################################
    def open_extract_variables(self):
        if not self.store:
            messagebox.showinfo("Extract Variables", "Please select a folder first (File → Select Folder).")
            return
        ExtractVariablesWindow(self, self.store)

    def open_records_viewer(self):
        if not self.store:
            messagebox.showinfo(
                "Records Viewer",
                "Please select a folder first (File → Select Folder)."
            )
            return
        if not self.store.data.get("patients"):
            messagebox.showinfo(
                "Records Viewer",
                "No patient records found in the current data store."
            )
            return
        RecordsViewerWindow(self, self.store)

        win = tk.Toplevel(self)
        win.title("Saved Patient Records History")
        win.geometry("widthxheight+X+Y")

        top_frame = ttk.Frame(win, padding=6)
        top_frame.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(top_frame, text="Search Patient ID/Name:").pack(side=tk.LEFT, padx=4)
        search_var = tk.StringVar()
        search_entry = ttk.Entry(top_frame, textvariable=search_var, width=30)
        search_entry.pack(side=tk.LEFT, padx=4)

        tree_frame = ttk.Frame(win, padding=6)
        tree_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        columns = ("metric", "value", "unit", "flag")
        tree = ttk.Treeview(tree_frame, columns=columns, show="tree headings")
        tree.heading("#0", text="Patient / Visit Date")
        tree.heading("metric", text="Laboratory Metric / Field")
        tree.heading("value", text="Value")
        tree.heading("unit", text="Unit")
        tree.heading("flag", text="Anomalies / Flags")

        tree.column("#0", width=240, anchor="w")
        tree.column("metric", width=220, anchor="w")
        tree.column("value", width=90, anchor="center")
        tree.column("unit", width=80, anchor="center")
        tree.column("flag", width=150, anchor="w")

        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)

        def populate_tree(*_args):
            tree.delete(*tree.get_children())
            query = search_var.get().strip().lower()

            for pid, patient in sorted(self.store.data.get("patients", {}).items()):
                p_name = patient.get("name", "")
                if query and (query not in pid.lower() and query not in p_name.lower()):
                    continue

                p_node = tree.insert("", tk.END, text=f"{pid} - {p_name}", open=False)
                
                visits = patient.get("visits", [])
                sorted_visits = sorted(visits, key=lambda v: v.get("visit_date") or "", reverse=True)
                
                for idx, visit in enumerate(sorted_visits):
                    v_date = visit.get("visit_date") or "Unknown Date"
                    img_count = len(visit.get("images", []))
                    v_node = tree.insert(p_node, tk.END, text=f"Visit: {v_date} ({img_count} Img)", open=False)
                    
                    for field in visit.get("fields", []):
                        tree.insert(
                            v_node, 
                            tk.END, 
                            text="", 
                            values=(
                                field.get("label", ""),
                                field.get("value", ""),
                                field.get("unit", ""),
                                field.get("flag", "")
                            )
                        )

        search_var.trace_add("write", populate_tree)
        populate_tree()

        btn_frame = ttk.Frame(win, padding=6)
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Button(btn_frame, text="Close", command=win.destroy).pack(side=tk.RIGHT, padx=4)

    def open_eda(self):
        if not self.store:
            messagebox.showinfo("EDA", "Please select a folder first (File → Select Folder).")
            return
        if not self.store.data.get("patients"):
            messagebox.showinfo("EDA", "No patient records found in the current data store.")
            return
        if not HAS_MATPLOTLIB:
            proceed = messagebox.askyesno(
                "Missing matplotlib",
                "matplotlib is not installed. You can still view the Summary Statistics "
                "table, but charts will not be available.\n\n"
                "Install it with:  pip install matplotlib\n\n"
                "Open EDA anyway?",
            )
            if not proceed:
                return
        EDAWindow(self, self.store)

    def open_alias_manager(self):
        if not self.store:
            messagebox.showinfo("No folder open", "Select a folder first - synonyms are saved per folder.")
            return

        win = tk.Toplevel(self)
        win.title("Manage Label Synonyms")
        win.geometry("640x480")

        ttk.Label(
            win,
            text="Different labs abbreviate the same test differently.\n"
                 "Built-in synonyms are always available. Entries you add or that\n"
                 "get learned automatically (by correcting a label) are saved to\n"
                 "label_aliases.json in this folder and used for future auto-fill.",
            justify="left"
        ).pack(side=tk.TOP, anchor="w", padx=10, pady=8)

        columns = ("alias", "canonical", "source")
        tree = ttk.Treeview(win, columns=columns, show="headings", height=15)
        tree.heading("alias", text="Alias (as typed/seen)")
        tree.heading("canonical", text="Maps to")
        tree.heading("source", text="Source")
        tree.column("alias", width=200)
        tree.column("canonical", width=260)
        tree.column("source", width=100, anchor="center")
        tree.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=4)

        def refresh_tree():
            tree.delete(*tree.get_children())
            for canonical, synonyms in sorted(BUILTIN_ALIASES.items()):
                for syn in synonyms:
                    tree.insert("", tk.END, iid=f"builtin::{normalize_label(syn)}",
                                values=(syn, canonical, "built-in"))
            for norm, entry in sorted(self.store.aliases.items(), key=lambda kv: kv[1]["example"]):
                tree.insert("", tk.END, iid=f"learned::{norm}",
                            values=(entry.get("example", norm), entry.get("canonical", ""), "learned/added"))

        refresh_tree()

        add_frame = ttk.LabelFrame(win, text="Add a synonym")
        add_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=8)

        ttk.Label(add_frame, text="Alias:").grid(row=0, column=0, padx=4, pady=4, sticky="w")
        alias_var = tk.StringVar()
        ttk.Entry(add_frame, textvariable=alias_var, width=20).grid(row=0, column=1, padx=4, pady=4)

        ttk.Label(add_frame, text="Maps to:").grid(row=0, column=2, padx=4, pady=4, sticky="w")
        canonical_var = tk.StringVar()
        canonical_combo = ttk.Combobox(add_frame, textvariable=canonical_var,
                                        values=LABEL_SUGGESTIONS, width=28)
        canonical_combo.grid(row=0, column=3, padx=4, pady=4)

        def add_alias():
            alias = alias_var.get().strip()
            canonical = canonical_var.get().strip()
            if not alias or not canonical:
                messagebox.showwarning("Missing info", "Enter both an alias and what it maps to.", parent=win)
                return
            if canonical not in LABEL_SUGGESTIONS:
                if not messagebox.askyesno(
                    "Not a known field",
                    f'"{canonical}" isn\'t one of the existing label suggestions.\n'
                    "Add it as a synonym target anyway?", parent=win
                ):
                    return
            self.store.learn_alias(alias, canonical)
            self.store.save_aliases()
            alias_var.set("")
            canonical_var.set("")
            refresh_tree()

        ttk.Button(add_frame, text="Add", command=add_alias).grid(row=0, column=4, padx=8, pady=4)

        def remove_selected():
            sel = tree.selection()
            if not sel:
                return
            changed = False
            for iid in sel:
                if iid.startswith("learned::"):
                    self.store.remove_alias(iid.split("::", 1)[1])
                    changed = True
            if changed:
                self.store.save_aliases()
                refresh_tree()

        bottom = ttk.Frame(win)
        bottom.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=8)
        ttk.Button(bottom, text="Remove Selected (learned/added only)",
                   command=remove_selected).pack(side=tk.LEFT)
        ttk.Button(bottom, text="Close", command=win.destroy).pack(side=tk.RIGHT)


def main():
    if PaddleOCR is None:
        print("WARNING: paddleocr is not installed. The app will open, but "
              "'Run OCR' will fail until you run:\n"
              "  pip install paddleocr paddlepaddle")
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()